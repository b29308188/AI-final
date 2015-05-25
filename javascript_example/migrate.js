// A sample player implemented in JavaScript, adapted from the Java
// Migrate player.  In this player, each pusher takes one of the red
// markers and tries to gradually move it to vertices that are on the
// boundary between red and non-red.
//
// ICPC Challenge
// Sturgill, NC State University

importPackage( java.lang );
importPackage( java.util );
 
// Width and height of the world in game units.
var FIELD_SIZE = 100;
  
// Number of pushers per side.
var PCOUNT = 3;

// Radius of the pusher. 
var PUSHER_RADIUS = 1;

// Mass of the pusher. 
var PUSHER_MASS = 1;

// Maximum velocity for a pusher.
var PUSHER_SPEED_LIMIT = 6.0;
  
// Maximum acceleration for a pusher.
var PUSHER_ACCEL_LIMIT = 2.0;
  
// Total number of markers on the field.
var MCOUNT = 22;

// Radius of the marker.
var MARKER_RADIUS = 2;

// Mass of the marker.
var MARKER_MASS = 3;

// Marker velocity lost per turn
var MARKER_FRICTION = 0.35;

// Width and height of the home region.
var HOME_SIZE = 20;

// Color value for the red player.
var RED = 0;
  
// Color value for the blue player.
var BLUE = 1;
  
// Color value for unclaimed pucks.
var GREY = 2;

// Object representing a 3d vertex on the field.
function Vertex3D( vx, vy, vz ) {
  this.x = vx;
  this.y = vy;
  this.z = vz;
}

// Object representing a 2d vector / point, along with userfule operations.
function Vector2D( vx, vy ) {
  this.x = vx;
  this.y = vy;
}

// Return the magnitude of this vector.
Vector2D.prototype.mag = function() {
  return Math.sqrt( this.x * this.x + this.y * this.y );
}

// Return a new vector containing normalized version of this vector.
Vector2D.prototype.norm = function() {
  var m = this.mag();
  return new Vector2D( this.x / m,
                       this.y / m );
}

// Return a ccw perpendicular vector for this vector
Vector2D.prototype.perp = function() {
  return new Vector2D( -this.y, this.x );
}

// Return a new vector containing a scaled by scaling factor s.
Vector2D.prototype.scale = function( s ) {
  return new Vector2D( this.x * s, this.y * s );
};


/* Return a vector pointing in the same direction as v, but with
   magnitude no greater than d. */
Vector2D.prototype.limit = function( d ) {
  var m = this.mag();
  if ( m > d )
    return new Vector2D( d * this.x / m, d * this.y / m );
  else
    return new Vector2D( this.x, this.y );
}

// Return a new vector thats the given vector rotated by r CCW.
Vector2D.prototype.rotate = function( r ) {
  var s = Math.sin( r );
  var c = Math.cos( r );
  return new Vector2D( this.x * c - this.y * s,
                             this.x * s + this.y * c );
}

// Return a new vector containing the sum of a and b.
function vecAdd( a, b ) {
  return new Vector2D( a.x + b.x, a.y + b.y );
}

// Return a new vector containing the sum of a and b.
function vecSub( a, b ) {
  return new Vector2D( a.x - b.x, a.y - b.y );
}

// Return the dot product of a and b.
function vecDot( a, b ) {
  return a.x * b.x + a.y * b.y;
}

// Return the cross product of a and b.
function cross( a, b ) {
  return a.x * b.y - a.y * b.x;
}

// Simple representation for a pusher.
function Pusher() {
  // Position of the pusher
  this.pos = new Vector2D( 0, 0 );
    
  // Pusher velocity
  this.vel = new Vector2D( 0, 0 );
    
  // True if this pusher has a job.
  this.busy = false;
    
  // How long we've been doing the current job.  If
  // this number gets to large, we'll pick a new job.
  this.jobTime = 0;
    
  // Target vertex for this pusher.
  targetVertex = 0;
};

// Simple representation for a marker.
function Marker() {
  // Position of the marker.
  this.pos = new Vector2D( 0, 0 );

  // Marker velocity
  this.vel = new Vector2D( 0, 0 );

  // Marker color
  this.color = GREY;
};

// Return the value of a, clamped to the [ b, c ] range
function clamp( a, b, c ) {
  if ( a < b )
    return b;
  if ( a > c )
    return c;
  return a;
}

/* One dimensional function to help compute acceleration
   vectors. Return an acceleration that can be applied to a pusher
   at pos and moving with velocity vel to get it to target.  The
   alim parameter puts a limit on the acceleration available.  This
   function is used by the two-dimensional moveTo function to
   compute an acceleration vector toward the target after movement
   perp to the target direction has been cancelled out.  */
function moveToLinear( pos, vel, target, alim ) {
  // Compute how far pos has to go to hit target.
  var dist = target - pos;

  // Kill velocity if we are close enough.
  if ( Math.abs( dist ) < 0.01 )
    return clamp( -vel, -alim, alim );
    
  // How many steps, at minimum, would cover the remaining distance
  // and then stop.
  var steps = Math.ceil(( -1 + Math.sqrt(1 + 8.0 * Math.abs(dist) / alim)) 
                           / 2.0);
  if ( steps < 1 )
    steps = 1;
    
  // How much acceleration would we need to apply at each step to
  // cover dist.
  var accel = 2 * dist / ( ( steps + 1 ) * steps );
    
  // Ideally, how fast would we be going now
  var ivel = accel * steps;

  // Return the best change in velocity to get vel to ivel.
  return clamp( ivel - vel, -alim, alim );
}

/** Print out a force vector that will move the given pusher to
    the given target location. */
function moveTo( p, target ) {
  // Compute a frame with axis a1 pointing at the target.
  var a1, a2;

  // Build a frame (a trivial one if we're already too close).
  var dist = vecSub( target, p.pos ).mag();
  if ( dist < 0.0001 ) {
    a1 = new Vector2D( 1.0, 0.0 );
    a2 = new Vector2D( 0.0, 1.0 );
  } else {
    a1 = vecSub( target, p.pos ).scale( 1.0 / dist );
    a2 = a1.perp();
  }
        
  // Represent the pusher velocity WRT that frame.
  var v1 = vecDot( a1, p.vel );
  var v2 = vecDot( a2, p.vel );

  // Compute a force vector in this frame, first cancel out velocity
  // perp to the target.
  var f1 = 0;
  var f2 = -v2;

  // If we have remaining force to spend, use it to move toward the target.
  if ( Math.abs( f2 ) < PUSHER_ACCEL_LIMIT ) {
    var raccel = Math.sqrt( PUSHER_ACCEL_LIMIT * PUSHER_ACCEL_LIMIT - 
                            v2 * v2 );
    f1 = moveToLinear( -dist, v1, 0.0, raccel );
  }

  // Convert force 
  var force = vecAdd( a1.scale( f1 ), a2.scale( f2 ) );
  System.out.print( force.x + " " + force.y );
}

/* Print out a force vector that will move the given pusher around
   to the side of marker m that's opposite from target.  Return true
   if we're alreay behind the marker.  */
function moveAround( p, m, target ) {
    // Compute vectors pointing from marker-to-target and Marker-to-pusher
  var mToT = vecSub( target, m.pos ).norm();
  var mToP = vecSub( p.pos, m.pos ).norm();
    
  // See if we're already close to behind the marker.
  if ( vecDot( mToT, mToP ) < -0.8 )
    return true;

  // Figure out how far around the target we need to go, we're
  // going to move around a little bit at a time so we don't hit
  // the target.
  var moveAngle = Math.acos( vecDot( mToT, mToP ) );
  if ( moveAngle > Math.PI * 0.25 )
    moveAngle = Math.PI * 0.25;

  // We're not, decide which way to go around.
  if ( cross( mToT, mToP ) > 0 ) {
    // Try to go around to the right.
    moveTo( p, vecAdd( m.pos, mToP.rotate( moveAngle ).scale( 4 ) ) );
  } else {
    // Try to go around to the left.
    moveTo( p, vecAdd( m.pos, mToP.rotate( -moveAngle ).scale( 4 ) ) );
  }
  
  return false;
}

// Scanner to parse input from the game engine.
var input = new Scanner( System[ 'in' ] );

// current score for each player.
var score = [ 0, 0 ];

// Read the list of vertex locations.
var n = input.nextInt();
// List of points in the map.
var vertexList = new Array();
for ( var i = 0; i < n; i++ ) {
  vertexList[ i ] = new Vertex3D();
  vertexList[ i ].x = input.nextInt();
  vertexList[ i ].y = input.nextInt();
  vertexList[ i ].z = input.nextInt();
 }

// Read the list of region outlines.
n = input.nextInt();
// List of regions in the map
var regionList = new Array();
for ( var i = 0; i < n; i++ ) {
  var m = input.nextInt();;
  regionList[ i ] = new Array();
  for ( var j = 0; j < m; j++ )
    regionList[ i ][ j ] = input.nextInt();
}

// List of current region colors
var regionColors = new Array( regionList.length );

// List of pusher objects.
var pList = new Array( 2 * PCOUNT );
for ( var i = 0; i < pList.length; i++ )
  pList[ i ] = new Pusher();

// List of marker objects
var mList = new Array( MCOUNT );
for ( var i = 0; i < mList.length; i++ )
  mList[ i ] = new Marker();

var turnNum = input.nextInt();
while ( turnNum >= 0 ) {
  score[ RED ] = input.nextInt();
  score[ BLUE ] = input.nextInt();

  // Read all the region colors.
  n = input.nextInt();
  for ( var i = 0; i < regionList.length; i++ )
    regionColors[ i ] = input.nextInt();

  // Read all the pusher locations.
  n = input.nextInt();
  for ( var i = 0; i < pList.length; i++ ) {
    pList[ i ].pos.x = input.nextDouble();
    pList[ i ].pos.y = input.nextDouble();
    pList[ i ].vel.x = input.nextDouble();
    pList[ i ].vel.y = input.nextDouble();
  }

  // Read all the marker locations.
  n = input.nextInt();
  for ( var i = 0; i < n; i++ ) {
    mList[ i ].pos.x = input.nextDouble();
    mList[ i ].pos.y = input.nextDouble();
    mList[ i ].vel.x = input.nextDouble();
    mList[ i ].vel.y = input.nextDouble();
    mList[ i ].color = input.nextInt();
  }
    
  // Compute a bit vector for the region colors incident on each
  // vertex.
  var vertexColors = new Array( vertexList.length );
  for ( var i = 0; i < vertexColors.length; i++ )
    vertexColors[ i ] = 0;
  for ( var i = 0; i < regionList.length; i++ )
    for ( var j = 0; j < regionList[ i ].length; j++ )
      vertexColors[ regionList[ i ][ j ] ] |= ( 1 << regionColors[ i ] );
  
  // Candidate vertices for putting a marker on, vertices that have
  // some red but are not all red.
  var candidates = new Array();
  for ( var i = 0; i < vertexList.length; i++ )
    if ( ( vertexColors[ i ] & 0x1 ) == 1 &&
         vertexColors[ i ] != 1  )
      candidates.push( i );

  // Choose a next action for each pusher, each pusher is responsible
  // for the marker with the same index.
  for ( var pdex = 0; pdex < PCOUNT; pdex++ ) {
    var p = pList[ pdex ];
      
    // See how long this pusher has been doing its job.
    if ( p.busy ) {
      // Go to idle if we work to long on the same job.
      p.jobTime++;
      if ( p.jobTime >= 60 )
        p.busy = false;
    }

    // If we lose our marker, then just sit idle.
    if ( mList[ pdex ].color != RED ) {
      p.busy = false;
    }
        
    // Otherwise, try to find a new place to push our marker.
    if ( mList[ pdex ].color == RED &&
         !p.busy ) {
      if ( candidates.length > 0 ) {
        var choice = Math.floor( Math.random() * candidates.length );
        p.targetVertex = candidates[ choice ];
        candidates.splice( choice, 1 );
        p.busy = true;
      }
    }

    // Choose a move direction in support of our current goal.
    if ( p.busy ) {
      // Get behind our marker and push it toward its destination.
      var v = vertexList[ p.targetVertex ];
      var dest = new Vector2D( v.x, v.y );
      if ( moveAround( p, mList[ pdex ], dest ) ) {
        var mToD = vecSub( dest, mList[ pdex ].pos ).norm();
        moveTo( p, vecSub( mList[ pdex ].pos, mToD ) );
      }
    } else
      System.out.print( "0.0 0.0" );

    // Print a space or a newline depending on whether we're at
    // the last pusher.
    if ( pdex + 1 < PCOUNT )
      System.out.print( " " );
    else
      System.out.println();
  }

  turnNum = input.nextInt();
}

