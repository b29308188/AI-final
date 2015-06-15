
# A simple player based on the C++ example player, gather.  Feel free
# to use this as a starting point for your own python player.  The
# player uses the home region to try to convert markers to the
# player's color.  If the home region has enough of the player's
# markers touching it, the player tries to move them elsewhere on the
# field.
#
# I'm not really a Python programmer, so feel free to send in any
# suggestions you have.
#
# ICPC Challenge
# Sturgill, NC State University

import random
import sys
import string
import math
import numpy as np
import sys

# Width and height of the world in game units.
FIELD_SIZE = 100
  
# Number of pushers per side.
PCOUNT = 3

# Number of total markers
MCOUNT = 22

# Radius of the pusher.
PUSHER_RADIUS = 1

# Mass of the pusher.
PUSHER_MASS = 1

# Maximum velocity for a pusher.
PUSHER_SPEED_LIMIT = 6.0
  
# Maximum acceleration for a pusher.
PUSHER_ACCEL_LIMIT = 2.0
  
# Radius of the marker.
MARKER_RADIUS = 2

# Mass of the marker.
MARKER_MASS = 3

# Marker velocity lost per turn
MARKER_FRICTION = 0.35

# Width and height of the home region.
HOME_SIZE = 20

# Color value for the red player
RED = 0

# Color value for the blue player
BLUE = 1

# Neutral color
GREY = 2

# Action being performed on marker
NO_INFLUENCE = 0

MIGRATE = 1

GATHER = 2

# Source of randomness
rnd = random.Random()

# Simple 2D Point/Vector representation along with common utility
# functions.
class Vector2D:
    # Make a vector with the given x, y coordinates.
    def __init__( self, vx, vy ):
        self.x = vx
        self.y = vy

    # Store a copy of a in this vector
    def set( self, a ):
        self.x = a.x
        self.y = a.y

    # Return the squared magnitude of this vector.
    def squaredMag( self ):
        return self.x * self.x + self.y * self.y

    # Return the magnitude of this vector.
    def mag( self ):
        return math.sqrt( self.x * self.x + self.y * self.y )

    # Return a unit vector pointing in the same direction as this.
    def norm( self ):
        m = self.mag()
        return Vector2D( self.x / m, self.y / m )

    # Return a CCW perpendicular to this vector.
    def perp( self ):
        return Vector2D( -self.y, self.x )

    # Return a cross product of this and b.
    def cross( self, b ):
        return self.x * b.y - self.y * b.x

    # Return a vector pointing in the same direction as this, but with
    # magnitude no greater than d.
    def limit( self, d ):
        m = self.mag()
        if m > d:
          return Vector2D( d * self.x / m, d * self.y / m )
        else:
          return Vector2D( self.x, self.y )

# Return a vector that's the sum of a and b.
def vecAdd( a, b ):
    return Vector2D( a.x + b.x, a.y + b.y )

#* Return a vector that's a minus b.
def vecSub( a, b ):
    return Vector2D( a.x - b.x, a.y - b.y )

# Return a copy of a that's a scaled by b.
def vecScale( a, b ):
  return Vector2D( a.x * b, a.y * b )

# Return the dot product of a and b.
def vecDot( a, b ):
  return a.x * b.x + a.y * b.y

# Print out a copy of a and b.
def vecPrint( a ):
    sys.stdout.write( "( %f %f )" * ( a.x, a.y ) )

# Return the vector from vectors that is nearest to a
def nearestVector(a, vectors):
    return min(vectors, key=lambda v:vecSub(a.pos, v.pos).mag())

# Simple 3D Point/Vector representation.
class Vector3D:
    # Initialize with given coordinates. */
    def __init__( self, vx, vy, vz ):
        self.x = vx
        self.y = vy
        self.z = vz

def compute_force(p): 
    force = Vector2D(0, 0)
    marker = p.target_marker
    if getBehind(p, marker, p.targetPos, force):
        runTo( p.pos, p.vel, 
               vecSub( marker.pos, vecSub( p.targetPos, marker.pos ).norm() ),
               force, 0.1 )
    return force

# return true if marker is 
def being_coerced(marker):
    return len(marker.touching_regions) >= 1 and\
           len([1 for r in marker.touching_regions if r.my_weight < 1]) == 0 and\
           len([1 for r in marker.touching_regions if r.color != RED]) == 0  

class State:
    def __init__(self, board, markers, my_pushers, 
                 opp_pushers, turn_no, scores):
        self.board = board
        self.markers = markers
        self.my_pushers = my_pushers
        self.opp_pushers = opp_pushers
        self.turn_no = turn_no
        self.scores = scores

    def update_notification(self):
        # Process home, opposing markers.
        self.my_markers = [m for m in self.markers if m.color == RED]
        self.opp_markers = [m for m in self.markers if m.color == BLUE]
        self.neut_markers = [m for m in self.markers if m.color == GREY]
        # Process home, opposing region colors.
        self.board.my_regions =\
        [r for r in self.board.regions if r.color == RED]
        self.board.opp_regions =\
        [r for r in self.board.regions if r.color == BLUE]
        self.board.neut_regions =\
        [r for r in self.board.regions if r.color == GREY]


        for r in self.board.regions:
            r.markers_touching = []
        # Process markers in regions. 
        for m in self.markers:
            m.touching_regions = []
            # Each circle is approximated as touching a 4X4 square. Iterate
            # over 16 points representing circle
            for x in xrange(int(m.pos.x-2),int(m.pos.x+2)):
                for y in xrange(int(m.pos.y-2), int(m.pos.y+2)):
                    # Add regions touched by point to regions touched by
                    # marker.
                    point_touched = self.board.regions_touched[x][y] 
                    m.touching_regions = list(set(m.touching_regions+\
                    point_touched))

            # Indicate to region that marker is touching it
            for r in m.touching_regions:
                r.markers_touching.append(m)

        color_weight = {RED:1, BLUE:-1, GREY:-1}
        for r in self.board.regions:
            r.my_weight = sum([color_weight[m.color] for m in r.markers_touching])

        #print >> sys.stderr, len(self.board.regions[0].markers_touching)

class Board:
    def __init__(self, regions, vertices):
        self.regions = regions
        self.vertices = vertices
        self.init_regions_touched()

    def init_regions_touched(self):
        # list of regions touched by every point on the map
        self.regions_touched = [[[] for y in range(101)] for x in range(101)]
        for r in self.regions:
            min_y = min(r.vertices, key=lambda v:v.pos.y).pos.y
            max_y = max(r.vertices, key=lambda v:v.pos.y).pos.y
            # Iterate over horizontal lines of polygon.
            for y in xrange(min_y, max_y+1):
                verts = r.vertices
                s = lambda a,b: (a,b) if a.x < b.x else (b,a)
                edges = [s(v[0].pos,v[1].pos) for v in zip(verts, verts[1:]+[verts[0]])]
                edges_with_y = [e for e in edges if between(y, e[0].y, e[1].y)]
                x_vals = [x_on_line(e,y) for e in edges_with_y]
                x_vals.sort()
                # mark every point on horizontal line as touching this region
                for x in xrange(x_vals[0], x_vals[1]+1):
                    if r not in self.regions_touched[x][y]:
                        self.regions_touched[x][y].append(r)
                                

class Region:
    def __init__(self):
        self.color = GREY
        self.markers_touching = []
        self.my_weight = 0
        self.vertices = []
        self.center = Vector2D(0,0)
        

class Vertex:
    def __init__(self):
        self.pos = Vector3D(0,0,0)
        # colors is a bitmap of colors incident on this vertex.
        self.colors = 0


# Simple representation for a marker.
class Marker:
    def __init__( self ):
        self.pos = Vector2D( 0, 0 )
        self.vel = Vector2D( 0, 0 )
        self.color = GREY
        self.shapely_point = None
        self.gather = False
        self.current_influence = NO_INFLUENCE
        #self.regions = []
        self.touching_regions = []

# Simple representation for a pusher.
class Pusher:
    def __init__( self ):
        # Position of the Pusher.
        self.pos = Vector2D(0, 0)

        # Pusher velocity
        self.vel = Vector2D(0, 0)

        # True if this pusher has a job.
        self.busy = False

        # How long we've been doing the current job.  If
        # this number gets to large, we'll pick a new job.
        self.jobTime = 0

        # Index of the marker each pusher is working with.
        self.mdex = 0

        # The marker that this pusher is trying to move.
        self.target_marker = None

        # Location the pusher is trying to move it's marker to.
        self.targetPos = Vector2D( 0, 0 )



class Agent:
    def get_actions(self, state):        
        raise NotImplementedError("") 
        

class GatherMigrateAgent(Agent):
    def get_actions(self, state):

        for p in state.my_pushers:
            self.check_pusher_done(p)
       
        free_pushers = [p for p in state.my_pushers if p.busy == False]
        # Candidate migrate markers are our markers that we are not 
        # currently not migrating
        cand_migrate_markers = [m for m in state.my_markers if\
                                m.current_influence == NO_INFLUENCE]
        # Candidate vertices are vertices touched by red regions, but not only
        # touched by red regions. There should always be at least 3 candidate
        # vertices because of opponent's home region
        cand_verts = [v for v in state.board.vertices\
                      if (v.colors & 1 << RED) == 1 and v.colors != 1 << RED]

        # Candidate gather markers are markers not of our color that we are
        # currently not pushing or coercing.
        cand_gather_markers = [m for m in 
                               (state.neut_markers+state.opp_markers)\
                               if m.current_influence == NO_INFLUENCE and\
                               not being_coerced(m)]
        if len(cand_gather_markers) == 0:
            cand_gather_markers = state.markers

        cand_gather_regions = [r for r in state.board.my_regions if r.my_weight >= 1]
        if len(cand_gather_regions) == 0:
            cand_gather_regions = [state.board.regions[0]]

        for p in free_pushers:
            p.jobTime = 0
            p.busy = True
            # XXX need better heuristic for picking between migrate and gather
            percent_markers = float(len(state.my_markers)) / float(MCOUNT)
            percent_regions = float(len(state.board.my_regions)) / float(len(state.board.regions))

#            if(percent_markers > percent_regions and \
#               len(cand_migrate_markers) >= 1):
            if(len(cand_migrate_markers) > len(state.board.my_regions)):
                # Migrate nearest marker to candidate vertex.
                nearest_migrate_marker = min(cand_migrate_markers, 
                                     key=lambda m:vecSub(p.pos, m.pos).mag())
                cand_migrate_markers.remove(nearest_migrate_marker)
                p.target_marker = nearest_migrate_marker
                target_vert = cand_verts.pop()
                p.targetPos = Vector2D(target_vert.pos.x, target_vert.pos.y)
                nearest_migrate_marker.current_influence = MIGRATE
            else:
                # Gather nearest marker to pusher, and push to nearest
                # already captured region.
                nearest_gather_marker = min(cand_gather_markers, 
                                     key=lambda m:vecSub(p.pos, m.pos).mag())
                cand_gather_markers.remove(nearest_gather_marker)

                nearest_region = min(cand_gather_regions, 
                                     key=lambda r:vecSub(
                                     nearest_gather_marker.pos, 
                                     r.center).mag())

                p.target_marker = nearest_gather_marker
                nearest_gather_marker.current_influence = GATHER
                p.targetPos = Vector2D(nearest_region.center.x,
                                       nearest_region.center.y)
 
        actions = [compute_force(p) for p in state.my_pushers]
        return actions
            
    def check_pusher_done(self, pusher): 
        if(pusher.busy): 
            pushed_marker = pusher.target_marker
            # Free pusher if gather marker is already being coerced.
            cond1 = pushed_marker.current_influence == GATHER and\
               being_coerced(pushed_marker)
            # Free pusher if it has been working on marker too long, or if
            # it has already pushed its marker to the destination
            cond2 = pusher.jobTime > 75 or\
               vecSub( pusher.target_marker.pos, pusher.targetPos ).mag() < 5
            if(cond1 or cond2):
                pushed_marker.current_influence = NO_INFLUENCE
                pusher.busy = False

# Return a copy of x that's constrained to be between low and high.
def clamp( x, low, high ):
  if x < low:
    x = low
  if x > high:
    x = high
  return x

def between(x, a, b): 
    return (a != b) and ((a <= x and x <= b) or (b <= x and x <= a))

# Return x value of point with y value k that is on line e. e is (a,b),
# where a and b are vertices with a.x <= b.x
def x_on_line(e, k):
    # check for vertical line
    if( e[1].x - e[0].x == 0):
        return e[0].x
    slope = float(e[1].y - e[0].y) / float(e[1].x - e[0].x)
    return int(float(k-e[0].y)/slope + float(e[0].x))
    



# Compute a force vector that can be applied to a pusher to get it
# to run through the given target location.  Pos and vel are the
# pusher's current position and velocity.  Target is the position we
# want to run through and force is a returned vector that will move
# the pusher toward the target.  The function returns true until it
# looks like the next move will take us through the target
# location. 
def runTo( pos, vel, target, force, epsilon ):
    # Get a unit vector in the direction we need to move.
    direction = vecSub( target, pos ).norm()
    
    # First, cancel out any movement that is perpendicular to the desired
    # movement direction.
    perp = direction.perp()
    force.set( vecScale( perp, -vecDot( perp, vel ) ).limit( PUSHER_ACCEL_LIMIT ) )

    # Use all the residual force to move toward the target.
    resForce = PUSHER_ACCEL_LIMIT - force.mag()
    force.set( vecAdd( force, vecScale( direction.norm(), resForce ) ) )

    # See if this move will cross close enough to the target location.
    nvel = vecAdd( vel, force ).limit( PUSHER_SPEED_LIMIT )
    t = clamp( vecDot( vecSub( target, pos ), nvel ) / vecDot( nvel, nvel ), 0, 1 )
    if vecSub( vecAdd( pos, vecScale( nvel, t ) ), target ).mag() < epsilon:
        return True
    
    return False

# Fill in the given force vector with a force intended to move the
# given pusher around behind the given marker so that it can be
# pushed toward the target destination.  Return true if the pusher
# is already behind the marker.
def getBehind( p, m, target, force ):
    # Make sure we're behind the target marker.
    mToT = vecSub( target, m.pos ).norm()
    pToM = vecSub( m.pos, p.pos ).norm()

    # See if we're already behind the marker.
    if vecDot( mToT, pToM ) > 0.7:
        return True

    # We're not, decide which way to go around.
    if pToM.cross( mToT ) > 0:
        # Try to go around to the right.
        force.set( vecScale( pToM.perp(), -PUSHER_ACCEL_LIMIT ) )
    else:
        # Try to go around to the left.
        force.set( vecScale( pToM.perp(), PUSHER_ACCEL_LIMIT ) )

    # Try to get closer to the marker if we're far away.
    maxDist = 8.0
    minDist = 6.0
    dist = vecSub( m.pos, p.pos ).mag()
    # Add a vector to help move in or out to get to the right distance.
    # from the marker.
    if dist > maxDist:
        force.set( vecAdd( force, vecScale( pToM, dist - maxDist ) ) )
    elif dist < minDist:
        force.set( vecSub( force, vecScale( pToM, minDist - dist ) ) )
    else:
        # cancel out any inward/outward velocity if the distance is good.
        inward = vecDot( p.vel, pToM )
        force.set( vecSub( force, vecScale( pToM, inward ) ) )
    return False

# Return true if the given marker is my color and is touching my home
# region.
def atHome( m ):
    return ( m.color == RED and
             m.pos.x < HOME_SIZE + MARKER_RADIUS and
             m.pos.y < HOME_SIZE + MARKER_RADIUS )

# Return a random field locatmion where we could move a marker.
def randomFieldPosition():
    return Vector2D( rnd.uniform( MARKER_RADIUS, FIELD_SIZE - MARKER_RADIUS ),
                     rnd.uniform( MARKER_RADIUS, FIELD_SIZE - MARKER_RADIUS ) )

#######################################################################

def create_state():
    # Read the static parts of the map.
    vertices = read_vertices()
    regions = read_regions(vertices)
    board = Board(regions, vertices)

    # These are updated on every turn snapshot from the game.
    pushers = tuple([ Pusher() for x in xrange( 2 * PCOUNT ) ])
    markers = tuple([Marker() for i in xrange(MCOUNT)])
    turn_no = int( sys.stdin.readline() )
    
    scores = (0,0)
    state = State(board, markers, pushers[:PCOUNT], pushers[PCOUNT:], 
                  turn_no, scores)
    return state

def read_vertices():
    n = int( sys.stdin.readline() )
    vertices = tuple([Vertex() for i in xrange(n)])
    # read list of points in the map.
    for vertex in vertices:
        vert_tokens = [int(i) for i in string.split( sys.stdin.readline() )]
        vertex.pos = Vector3D(*tuple(vert_tokens))    
    return vertices

def read_regions(vertices): 
    n = int(sys.stdin.readline())
    regions = tuple([Region() for i in xrange(n)]) 
    for region in regions:
        vert_indices =[int(i) for i in string.split(sys.stdin.readline())][1:]
        region.vertices = [vertices[i] for i in vert_indices]
        region.center.x = np.mean([v.pos.x for v in region.vertices])
        region.center.y = np.mean([v.pos.y for v in region.vertices])
    return regions

def update_state(state):
    scores = tuple([int(score) for score in string.split(sys.stdin.readline())])
    state.scores = scores
    update_region_colors(state.board.regions)
    update_pushers(state.my_pushers + state.opp_pushers)
    update_markers(state.markers)
    state.update_notification()

def update_region_colors(regions):
    # update regions' colors
    colorTokens = string.split(sys.stdin.readline())
    regionColors = tuple([int(color) for color in colorTokens[1:]]) 
    for i,region in enumerate(regions):
        region.color = regionColors[i]
    # update vertices' incident colors
    for region in regions:
        for vertex in region.vertices:
            vertex.colors |= 1 << region.color

# return tuple (myPushers, oppPushers). myPushers, oppPushers are each a tuple
# of pushers. myPushers are of my color. oppPushers are of the opponents' color
def update_pushers(pushers):
    # Read number of pushers. Value is already known to be 6
    sys.stdin.readline()
    # Update each pusher with position, velocity read from std in.
    for p in pushers:
        tokens = string.split( sys.stdin.readline() )
        p.pos.x = float(tokens[0])
        p.pos.y = float(tokens[1])
        p.vel.x = float(tokens[2])
        p.vel.y = float(tokens[3])
        p.jobTime += 1

def update_markers(markers):
    # Read number of markers. Value is already known.
    sys.stdin.readline()
    # Update each marker with position, velocity, and color read from std in.
    for m in markers:
        tokens = string.split( sys.stdin.readline() ) 
        m.pos.x = float(tokens[0])
        m.pos.y = float(tokens[1])
        m.vel.x = float(tokens[2])
        m.vel.y = float(tokens[3])
        m.color = int(tokens[4])

#######################################################################

state = create_state()
while state.turn_no >= 0:

    update_state(state)
    agent = GatherMigrateAgent()
    actions = agent.get_actions(state)
    # Write forces corresponding to pushers' next actions
    for force in actions:
        sys.stdout.write("%f %f "%(force.x, force.y ))

    sys.stdout.write( "\n" )
    sys.stdout.flush()
    turn_no = int( sys.stdin.readline() )


