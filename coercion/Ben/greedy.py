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
# Width and height of the world in game units.
FIELD_SIZE = 100
  
# Number of pushers per side.
PCOUNT = 3

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

# Simple representation for a pusher.
class Pusher:
	def __init__( self ):
		# Position of the Pusher.
		self.pos = Vector2D( 0, 0 )

		# Pusher velocity
		self.vel = Vector2D( 0, 0 )

		# True if this pusher has a job.
		self.busy = False

		# How long we've been doing the current job.  If
		# this number gets to large, we'll pick a new job.
		self.jobTime = 0

		# Index of the marker each pusher is working with.
		self.mdex = 0

		# Location the pusher is trying to move it's marker to.
		self.targetPos = Vector2D( 0, 0 )

# Simple representation for a marker.
class Marker:
	def __init__( self ):
		# Position of the marker.
		self.pos = Vector2D( 0, 0 )

		# Marker velocity
		self.vel = Vector2D( 0, 0 )

		# Marker color
		self.color = GREY

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
	return ( m.pos.x < HOME_SIZE + MARKER_RADIUS and m.pos.y < HOME_SIZE + MARKER_RADIUS )
def isRedHome(x, y):
	return ( x < HOME_SIZE  and y < HOME_SIZE )
def isBlueHome(x, y):
	return ( x > FIELD_SIZE- HOME_SIZE  and y > FIELD_SIZE- HOME_SIZE )

# Return a random field locatmion where we could move a marker.
def randomFieldPosition():
	return Vector2D( rnd.uniform( MARKER_RADIUS, FIELD_SIZE - MARKER_RADIUS ),
					 rnd.uniform( MARKER_RADIUS, FIELD_SIZE - MARKER_RADIUS ) )

#######################################################################

# current score for each player.
score = [ 0, 0 ]

# Read the static parts of the map.

# Read the list of vertex locations.
n = int( sys.stdin.readline() )
# List of points in the map.
vertexList = []
for i in range( n ):
	tokens = string.split( sys.stdin.readline() )
	vertexList.append( Vector3D( int( tokens[ 0 ] ), 
								 int( tokens[ 1 ] ),
								 int( tokens[ 2 ] ) ) )

# Read the list of region outlines.
n = int( sys.stdin.readline() )
# List of regions in the map
regionList = []
for i in range( n ):
	tokens = string.split( sys.stdin.readline() )
	# build a vertex list for the region.
	m = int( tokens[ 0 ] )
	regionList.append( [] )
	for j in range( m ):
		regionList[ i ].append( int( tokens[ j + 1 ] ) )

# List of current region colors, pusher and marker locations.
# These are updated on every turn snapshot from the game.
regionColors = [ GREY for x in regionList ]
pList = [ Pusher() for x in range( 2 * PCOUNT ) ]
mList = []

log_file = open("log", "w")
turnNum = int( sys.stdin.readline() )
while turnNum >= 0:
	#sys.stderr.write("======TURN  %d======\n"%(turnNum))
	tokens = string.split( sys.stdin.readline() )
	score[ RED ] = int( tokens[ 0 ] )
	score[ BLUE ] = int( tokens[ 1 ] )

	# Read all the region colors.
	tokens = string.split( sys.stdin.readline() )
	for i in range( len( regionList ) ):
		regionColors[ i ] = int( tokens[ i + 1 ] )
		
	# Read all the pusher locations.
	n = int( sys.stdin.readline() )
	for i in range( len( pList ) ):
		tokens = string.split( sys.stdin.readline() )
		pList[ i ].pos.x = float( tokens[ 0 ] )
		pList[ i ].pos.y = float( tokens[ 1 ] )
		pList[ i ].vel.x = float( tokens[ 2 ] )
		pList[ i ].vel.y = float( tokens[ 3 ] )

	# Read all the marker locations.
	n = int( sys.stdin.readline() )
	mList = [ Marker() for x in range( n ) ]
	for i in range( n ):
		tokens = string.split( sys.stdin.readline() )
		mList[ i ].pos.x = float( tokens[ 0 ] )
		mList[ i ].pos.y = float( tokens[ 1 ] )
		mList[ i ].vel.x = float( tokens[ 2 ] )
		mList[ i ].vel.y = float( tokens[ 3 ] )
		mList[ i ].color = int( tokens[ 4 ] )
	
	# Choose a next action for each pusher.
	for pdex in range( PCOUNT ):
		p = pList[ pdex ]
	  
		# See how long this pusher has been doing its job.
		if p.busy:
			#sys.stderr.write(str(pdex)+" is Busy at turn"+str(turnNum)+"\n")
			# Go to idle if we work to long on the same job.
			p.jobTime += 1
			#sys.stderr.write("JOB time :"+str(p.jobTime)+"\n")
			if p.jobTime >= 35:
				p.busy = False

			# Go back to idle if we finish our job.
			if vecSub( mList[ p.mdex ].pos, p.targetPos ).mag() < 3:
				p.busy = False
		if not p.busy:
			#sys.stderr.write(str(pdex)+" NOT Busy at turn"+str(turnNum)+"\n")
			searchSpace = [(i, j) for i in range(len(mList)) for j in range(len(regionList))]
			D = {}
			P = {}
			for ri in range(len(regionList)):
				x = np.mean([vertexList[v].x for v in regionList[ri]]) 
				y = np.mean([vertexList[v].y for v in regionList[ri]]) 
				P[ri] = (x, y)

			for (mi , ri) in searchSpace:
				(x, y) = P[ri]
				moveDistance = ( abs(mList[mi].pos.x -x) + abs(mList[mi].pos.y -y) ) / (FIELD_SIZE*2)
				numRedMarkers = 0.0
				numBlueMarkers = 0.0
				numRedRegions = 0.0
				numBlueRegions = 0.0
				bias = 0.0
				if mList[mi].color is BLUE and regionColors[ri] is RED:
					numRedMarkers += 1.0/40
					numBlueMarkers -= 1.0/40
					if not isRedHome(x, y):
						numRedRegions -= 1.0/20
				elif mList[mi].color is BLUE and regionColors[ri] is GREY:
					numBlueMarkers -= 1.0/40
					numBlueRegions += 1.0/20
				elif mList[mi].color is RED and regionColors[ri] is BLUE:
					numBlueMarkers += 1.0/40
					numRedMarkers -= 1.0/40
					if not isBlueHome(x, y):
						numBlueRegions -= 1.0/20
				elif mList[mi].color is RED and regionColors[ri] is GREY:
					numRedMarkers -= 1.0/40
					if atHome(mList[mi]): 	
						numRedRegions += 1.0/20
					bias -= moveDistance/20 
				elif mList[mi].color is GREY and regionColors[ri] is BLUE:
					numBlueMarkers += 1.0/40
					if not isBlueHome(x, y):
						numBlueRegions -= 1.0/20
				elif mList[mi].color is GREY and regionColors[ri] is RED:
					numRedMarkers += 1.0/40
					if not isRedHome(x, y):
						numRedRegions -= 1.0/20

				t = turnNum / 900.0
				D[(mi, x, y)] =  10*(1-t)*numRedMarkers - (1-t)*numBlueMarkers +20*t*numRedRegions - 20*t*numBlueRegions+bias
				
			available = False
			while available is False:
				(mdex, x, y) = max(D, key = D.get)
				available = True
				for j in range( PCOUNT ):
					if  j != pdex and pList[ j ].busy and pList[ j ].mdex == mdex :
						available = False
						D[(mdex, x, y)] = -np.inf
						break
			p.mdex = mdex
			p.targetPos = Vector2D(x, y)
			#sys.stderr.write(str(x)+"\t"+str(y)+"\n")
			p.busy = True
			p.jobTime = 0
			#sys.stderr.write(str(pdex)+":"+str(p.mdex)+"\t"+str(p.targetPos.x)+"\t"+str(p.targetPos.y)+"\n")

		# Choose a move direction in support of our current goal.
		force = Vector2D( 0, 0 )
		if p.busy:
			marker = mList[ p.mdex ]
		
			if getBehind( p, marker, p.targetPos, force ):
				runTo( p.pos, p.vel, 
					   vecSub( marker.pos, vecSub( p.targetPos, marker.pos ).norm() ),
					   force, 0.1 )

		sys.stdout.write( "%f %f " % ( force.x, force.y ) )

	sys.stdout.write( "\n" )
	sys.stdout.flush()

	turnNum = int( sys.stdin.readline() )
