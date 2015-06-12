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

class Board:
    def __init__(self, regions, vertices):
        self.regions = regions
        self.vertices = vertices
        self.my_regions = [r for r in self.regions if r.color == RED]
        self.opp_regions = [r for r in self.regions if r.color == BLUE]
        self.neut_regions = [r for r in self.regions if r.color == GREY]


class Region:
    def __init__(self):
        self.color = GREY
        #self.markers = []
        self.vertices = []

class Vertex:
    def __init__(self):
        self.pos = Vector3D(0,0,0)
        # colors is a bitmap of colors incident on this vertex.
        self.colors = 0


# Simple representation for a marker.
class Marker:
    def __init__( self ):
        # Position of the marker.
        self.pos = Vector2D( 0, 0 )

        # Marker velocity
        self.vel = Vector2D( 0, 0 )

        # Marker color
        self.color = GREY


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
        free_pushers = [p for p in state.my_pushers if self.is_free(p)] 
        moving_markers = [p.target_marker for p in state.my_pushers if not self.is_free(p)]   
        
        home_markers = [m for m in state.my_markers if atHome(m) and m not in moving_markers]
        cand_verts = [v for v in state.board.vertices\
                      if (v.colors & 1 << RED) == 1 and v.colors != 1 << RED]
            
        for h_m in home_markers[:len(free_pushers)]:
            p = free_pushers.pop()
            p.target_marker = h_m 
            # XXX doesn't check empty cand_vert
            target_vert = cand_verts.pop()
            p.targetPos = Vector2D(target_vert.pos.x, target_vert.pos.y)

        # otherwise, use pushers to gather markers 
        cand_markers = [m for m in (state.neut_markers+state.opp_markers)\
                        if m not in moving_markers]
        for p in free_pushers:
            # XXX doesn't check empty cand_markers
            target_marker = cand_markers.pop()
            p.target_marker = target_marker
            p.targetPos = Vector2D(10, 10)

        actions = [compute_force(p) for p in state.my_pushers]
        return actions
            
    def is_free(self, pusher):
        if(pusher.jobTime > 75 or pusher.target_marker == None or 
           vecSub( pusher.target_marker.pos, pusher.targetPos ).mag() < 5):
            return True
        return False

# Return a copy of x that's constrained to be between low and high.
def clamp( x, low, high ):
  if x < low:
    x = low
  if x > high:
    x = high
  return x

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


