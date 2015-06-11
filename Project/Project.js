//3456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_
// (JT: why the numbers? counts columns, helps me keep 80-char-wide listings)
//
// Chapter 5: ColoredTriangle.js (c) 2012 matsuda  AND
// Chapter 4: RotatingTriangle_withButtons.js (c) 2012 matsuda AND
// Chapter 2: ColoredPoints.js (c) 2012 matsuda
//
// merged and modified to became:
//
// ControlMulti.js for EECS 351-1, 
//									Northwestern Univ. Jack Tumblin

//		--converted from 2D to 4D (x,y,z,w) vertices
//		--demonstrate how to keep & use MULTIPLE colored shapes 
//			in just one Vertex Buffer Object(VBO).
//		--demonstrate several different user I/O methods: 
//				--Webpage pushbuttons 
//				--Webpage edit-box text, and 'innerHTML' for text display
//				--Mouse click & drag within our WebGL-hosting 'canvas'
//				--Keyboard input: alphanumeric + 'special' keys (arrows, etc)
//
// Vertex shader program----------------------------------
var VSHADER_SOURCE = 
  'uniform mat4 u_ModelMatrix;\n' +
  'attribute vec4 a_Position;\n' +
  'attribute vec4 a_Color;\n' +
  'varying vec4 v_Color;\n' +
  'void main() {\n' +
  '  gl_Position = u_ModelMatrix * a_Position;\n' +
  '  gl_PointSize = 10.0;\n' +
  '  v_Color = a_Color;\n' +
  '}\n';

// Fragment shader program----------------------------------
var FSHADER_SOURCE = 
  '#ifdef GL_ES\n' +
  'precision mediump float;\n' +
  '#endif GL_ES\n' +
  'varying vec4 v_Color;\n' +
  'void main() {\n' +
  '  gl_FragColor = v_Color;\n' +
  '}\n';

// Global Variables for the spinning tetrahedron:
var ANGLE_STEP = 45.0;  // default rotation angle rate (deg/sec)
var floatsPerVertex = 7;
// Global vars for mouse click-and-drag for rotation.
var isDrag=false;		// mouse-drag: true when user holds down mouse button
var xMclik=-0.4;			// last mouse button-down position (in CVV coords)
var yMclik=0.3;   
var xMdragTot=0.0;	// total (accumulated) mouse-drag amounts (in CVV coords).
var yMdragTot=0.0;  
var lrRotation=0.50;
var xpp=0.0;
var ypp=0.0;
var follow=true;
function main() {
//==============================================================================
  // Retrieve <canvas> element
  var canvas = document.getElementById('webgl');

  // Get the rendering context for WebGL
  var gl = getWebGLContext(canvas);
  if (!gl) {
    console.log('Failed to get the rendering context for WebGL');
    return;
  }

  // Initialize shaders
  if (!initShaders(gl, VSHADER_SOURCE, FSHADER_SOURCE)) {
    console.log('Failed to intialize shaders.');
    return;
  }

  // Initialize a Vertex Buffer in the graphics system to hold our vertices
  var n = initVertexBuffer(gl);
  if (n < 0) {
    console.log('Failed to set the vertex information');
    return;
  }

	// Register the Mouse & Keyboard Event-handlers-------------------------------
	// If users move, click or drag the mouse, or they press any keys on the 
	// the operating system will sense them immediately as 'events'.  
	// If you would like your program to respond to any of these events, you must // tell JavaScript exactly how to do it: you must write your own 'event 
	// handler' functions, and then 'register' them; tell JavaScript WHICH 
	// events should cause it to call WHICH of your event-handler functions.
	//
	// First, register all mouse events found within our HTML-5 canvas:
  canvas.onmousedown	=	function(ev){myMouseDown( ev, gl, canvas) }; 
  
  					// when user's mouse button goes down call mouseDown() function
  canvas.onmousemove = 	function(ev){myMouseMove( ev, gl, canvas) };
  
											// call mouseMove() function					
  canvas.onmouseup = 		function(ev){myMouseUp(   ev, gl, canvas)};
  					// NOTE! 'onclick' event is SAME as on 'mouseup' event
  					// in Chrome Brower on MS Windows 7, and possibly other 
  					// operating systems; use 'mouseup' instead.
  					
  // Next, register all keyboard events found within our HTML webpage window:
	window.addEventListener("keydown", myKeyDown, false);
	window.addEventListener("keyup", myKeyUp, false);
	window.addEventListener("keypress", myKeyPress, false);
  // The 'keyDown' and 'keyUp' events respond to ALL keys on the keyboard,
  // 			including shift,alt,ctrl,arrow, pgUp, pgDn,f1,f2...f12 etc. 
  //			I find these most useful for arrow keys; insert/delete; home/end, etc.
  // The 'keyPress' events respond only to alpha-numeric keys, and sense any 
  //  		modifiers such as shift, alt, or ctrl.  I find these most useful for
  //			single-number and single-letter inputs that include SHIFT,CTRL,ALT.

	// END Mouse & Keyboard Event-Handlers-----------------------------------
	
  // Specify the color for clearing <canvas>
  gl.clearColor(0.0, 0.0, 0.0, 1.0);

	// NEW!! Enable 3D depth-test when drawing: don't over-draw at any pixel 
	// unless the new Z value is closer to the eye than the old one..
	gl.depthFunc(gl.LESS);
	gl.enable(gl.DEPTH_TEST); 	  
	
  // Get handle to graphics system's storage location of u_ModelMatrix
  var u_ModelMatrix = gl.getUniformLocation(gl.program, 'u_ModelMatrix');
  if (!u_ModelMatrix) { 
    console.log('Failed to get the storage location of u_ModelMatrix');
    return;
  }
  // Create a local version of our model matrix in JavaScript 
  var modelMatrix = new Matrix4();
  
  // Create, init current rotation angle value in JavaScript
  var currentAngle = 0.0;

  // ANIMATION: create 'tick' variable whose value is this function:
  //----------------- 
  var tick = function() {
    currentAngle = animate(currentAngle);  // Update the rotation angle
    draw(gl, n, currentAngle, modelMatrix, u_ModelMatrix);   // Draw shapes
//    console.log('currentAngle=',currentAngle); // put text in console.

//	Show some always-changing text in the webpage :  
//		--find the HTML element called 'CurAngleDisplay' in our HTML page,
//			 	(a <div> element placed just after our WebGL 'canvas' element)
// 				and replace it's internal HTML commands (if any) with some
//				on-screen text that reports our current angle value:
//		--HINT: don't confuse 'getElementByID() and 'getElementById()
	
		// Also display our current mouse-dragging state:
		//--------------------------------
    requestAnimationFrame(tick, canvas);   
    									// Request that the browser re-draw the webpage
    									// (causes webpage to endlessly re-draw itself)
  };
  tick();							// start (and continue) animation: draw current image
	
}
function makeCube(){
cubVerts=new Float32Array([
	 0.00, 0.00, 0.10, 1.00, 0.7, 0.7, 0.7,		// rectangular prism 
     0.20, 0.00, 0.10, 1.00, 0.3, 0.5, 0.3,
     0.20, 0.50, 0.10, 1.00, 0.7, 0.7, 0.7,
     0.00, 0.50, 0.10, 1.00, 0.3, 0.5, 0.3,
     0.00, 0.50,-0.05, 1.00, 0.7, 0.7, 0.7,
     0.20, 0.50,-0.05, 1.00, 0.3, 0.5, 0.3,
     0.20, 0.00,-0.05, 1.00, 0.7, 0.7, 0.7,    //7
     0.20, 0.00, 0.10, 1.00, 0.3, 0.5, 0.3,
     0.20, 0.50, 0.10, 1.00, 0.7, 0.7, 0.7,
     0.20, 0.50,-0.05, 1.00, 0.3, 0.5, 0.3,
     0.20, 0.00,-0.05, 1.00, 0.7, 0.7, 0.7,   //11
     0.00, 0.00,-0.05, 1.00, 0.3, 0.5, 0.3,
     0.00, 0.50,-0.05, 1.00, 0.7, 0.7, 0.7,
     0.00, 0.50, 0.10, 1.00, 0.3, 0.5, 0.3,
     0.00, 0.00, 0.10, 1.00, 0.7, 0.7, 0.7,
     0.00, 0.00,-0.05, 1.00, 0.3, 0.5, 0.3,
     ])
}
function makeSphere(r,g,b) {
//==============================================================================
// Make a sphere from one OpenGL TRIANGLE_STRIP primitive.   Make ring-like 
// equal-lattitude 'slices' of the sphere (bounded by planes of constant z), 
// and connect them as a 'stepped spiral' design (see makeCylinder) to build the
// sphere from one triangle strip.
  var slices = 12;		// # of slices of the sphere along the z axis. >=3 req'd
											// (choose odd # or prime# to avoid accidental symmetry)
  var sliceVerts	= 12;	// # of vertices around the top edge of the slice
											// (same number of vertices on bottom of slice, too)
  var topColr = new Float32Array([r, g, b]);	// North Pole: light gray
  var equColr = new Float32Array([r, g, b]);	// Equator:    bright green
  var botColr = new Float32Array([r, g, b]);	// South Pole: brightest gray.
  var sliceAngle = Math.PI/slices;	// lattitude angle spanned by one slice.

	// Create a (global) array to hold this sphere's vertices:
  sphVerts = new Float32Array(  ((slices * 2* sliceVerts) -2) * floatsPerVertex);
										// # of vertices * # of elements needed to store them. 
										// each slice requires 2*sliceVerts vertices except 1st and
										// last ones, which require only 2*sliceVerts-1.
										
	// Create dome-shaped top slice of sphere at z=+1
	// s counts slices; v counts vertices; 
	// j counts array elements (vertices * elements per vertex)
	var cos0 = 0.0;					// sines,cosines of slice's top, bottom edge.
	var sin0 = 0.0;
	var cos1 = 0.0;
	var sin1 = 0.0;	
	var j = 0;							// initialize our array index
	var isLast = 0;
	var isFirst = 1;
	for(s=0; s<slices; s++) {	// for each slice of the sphere,
		// find sines & cosines for top and bottom of this slice
		if(s==0) {
			isFirst = 1;	// skip 1st vertex of 1st slice.
			cos0 = 1.0; 	// initialize: start at north pole.
			sin0 = 0.0;
		}
		else {					// otherwise, new top edge == old bottom edge
			isFirst = 0;	
			cos0 = cos1;
			sin0 = sin1;
		}								// & compute sine,cosine for new bottom edge.
		cos1 = Math.cos((s+1)*sliceAngle);
		sin1 = Math.sin((s+1)*sliceAngle);
		// (Note we don't initialize j; grows with each new attrib,vertex, and slice)
		// go around the entire slice, generating TRIANGLE_STRIP verts
		if(s==slices-1) isLast=1;	// skip last vertex of last slice.
		for(v=isFirst; v< 2*sliceVerts-isLast; v++, j+=floatsPerVertex) {	
			if(v%2==0)
			{				// put even# vertices at the the slice's top edge
							// (why PI and not 2*PI? because 0 <= v < 2*sliceVerts
							// and thus we can simplify cos(2*PI(v/2*sliceVerts))  
				sphVerts[j  ] = sin0 * Math.cos(Math.PI*(v)/sliceVerts); 	
				sphVerts[j+1] = sin0 * Math.sin(Math.PI*(v)/sliceVerts);	
				sphVerts[j+2] = cos0;		
				sphVerts[j+3] = 1.0;			
			}
			else { 	// put odd# vertices around the slice's lower edge;
							// x,y,z,w == cos(theta),sin(theta), 1.0, 1.0
							// 					theta = 2*PI*((v-1)/2)/capVerts = PI*(v-1)/capVerts
				sphVerts[j  ] = sin1 * Math.cos(Math.PI*(v-1)/sliceVerts);		// x
				sphVerts[j+1] = sin1 * Math.sin(Math.PI*(v-1)/sliceVerts);		// y
				sphVerts[j+2] = cos1;																				// z
				sphVerts[j+3] = 1.0;																				// w.		
			}
			if(s==0) {	// finally, set some interesting colors for vertices:
				sphVerts[j+4]=topColr[0]; 
				sphVerts[j+5]=topColr[1]; 
				sphVerts[j+6]=topColr[2];	
				}
			else if(s==slices-1) {
				sphVerts[j+4]=botColr[0]; 
				sphVerts[j+5]=botColr[1]; 
				sphVerts[j+6]=botColr[2];	
			}
			else {
					sphVerts[j+4]=r;// equColr[0]; 
					sphVerts[j+5]=g;// equColr[1]; 
					sphVerts[j+6]=b;// equColr[2];					
			}
		}
	}
	return sphVerts
}
function makeTorus() {
//==============================================================================
// 		Create a torus centered at the origin that circles the z axis.  
// Terminology: imagine a torus as a flexible, cylinder-shaped bar or rod bent 
// into a circle around the z-axis. The bent bar's centerline forms a circle
// entirely in the z=0 plane, centered at the origin, with radius 'rbend'.  The 
// bent-bar circle begins at (rbend,0,0), increases in +y direction to circle  
// around the z-axis in counter-clockwise (CCW) direction, consistent with our
// right-handed coordinate system.
// 		This bent bar forms a torus because the bar itself has a circular cross-
// section with radius 'rbar' and angle 'phi'. We measure phi in CCW direction 
// around the bar's centerline, circling right-handed along the direction 
// forward from the bar's start at theta=0 towards its end at theta=2PI.
// 		THUS theta=0, phi=0 selects the torus surface point (rbend+rbar,0,0);
// a slight increase in phi moves that point in -z direction and a slight
// increase in theta moves that point in the +y direction.  
// To construct the torus, begin with the circle at the start of the bar:
//					xc = rbend + rbar*cos(phi); 
//					yc = 0; 
//					zc = -rbar*sin(phi);			(note negative sin(); right-handed phi)
// and then rotate this circle around the z-axis by angle theta:
//					x = xc*cos(theta) - yc*sin(theta) 	
//					y = xc*sin(theta) + yc*cos(theta)
//					z = zc
// Simplify: yc==0, so
//					x = (rbend + rbar*cos(phi))*cos(theta)
//					y = (rbend + rbar*cos(phi))*sin(theta) 
//					z = -rbar*sin(phi)
// To construct a torus from a single triangle-strip, make a 'stepped spiral' along the length of the bent bar; successive rings of constant-theta, using the same design used for cylinder walls in 'makeCyl()' and for 'slices' in makeSphere().  Unlike the cylinder and sphere, we have no 'special case' for the first and last of these bar-encircling rings.
//
var rbend = 1.5;										// Radius of circle formed by torus' bent bar
var rbar = 0.15;											// radius of the bar we bent to form torus
var barSlices = 23;									// # of bar-segments in the torus: >=3 req'd;
																		// more segments for more-circular torus
var barSides = 13;										// # of sides of the bar (and thus the 
																		// number of vertices in its cross-section)
																		// >=3 req'd;
																		// more sides for more-circular cross-section
// for nice-looking torus with approx square facets, 
//			--choose odd or prime#  for barSides, and
//			--choose pdd or prime# for barSlices of approx. barSides *(rbend/rbar)
// EXAMPLE: rbend = 1, rbar = 0.5, barSlices =23, barSides = 11.

	// Create a (global) array to hold this torus's vertices:
 torVerts = new Float32Array(floatsPerVertex*(2*barSides*barSlices +2));
//	Each slice requires 2*barSides vertices, but 1st slice will skip its first 
// triangle and last slice will skip its last triangle. To 'close' the torus,
// repeat the first 2 vertices at the end of the triangle-strip.  Assume 7

var phi=0, theta=0;										// begin torus at angles 0,0
var thetaStep = 2*Math.PI/barSlices;	// theta angle between each bar segment
var phiHalfStep = Math.PI/barSides;		// half-phi angle between each side of bar
																			// (WHY HALF? 2 vertices per step in phi)
	// s counts slices of the bar; v counts vertices within one slice; j counts
	// array elements (Float32) (vertices*#attribs/vertex) put in torVerts array.
	for(s=0,j=0; s<barSlices; s++) {		// for each 'slice' or 'ring' of the torus:
		for(v=0; v< 2*barSides; v++, j+=7) {		// for each vertex in this slice:
			if(v%2==0)	{	// even #'d vertices at bottom of slice,
				torVerts[j  ] = (rbend + rbar*Math.cos((v)*phiHalfStep)) * 
																						 Math.cos((s)*thetaStep);
							  //	x = (rbend + rbar*cos(phi)) * cos(theta)
				torVerts[j+1] = (rbend + rbar*Math.cos((v)*phiHalfStep)) *
																						 Math.sin((s)*thetaStep);
								//  y = (rbend + rbar*cos(phi)) * sin(theta) 
				torVerts[j+2] = -rbar*Math.sin((v)*phiHalfStep);
								//  z = -rbar  *   sin(phi)
				torVerts[j+3] = 1.0;		// w
			}
			else {				// odd #'d vertices at top of slice (s+1);
										// at same phi used at bottom of slice (v-1)
				torVerts[j  ] = (rbend + rbar*Math.cos((v-1)*phiHalfStep)) * 
																						 Math.cos((s+1)*thetaStep);
							  //	x = (rbend + rbar*cos(phi)) * cos(theta)
				torVerts[j+1] = (rbend + rbar*Math.cos((v-1)*phiHalfStep)) *
																						 Math.sin((s+1)*thetaStep);
								//  y = (rbend + rbar*cos(phi)) * sin(theta) 
				torVerts[j+2] = -rbar*Math.sin((v-1)*phiHalfStep);
								//  z = -rbar  *   sin(phi)
				torVerts[j+3] = 1.0;		// w
			}
			torVerts[j+4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j+5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j+6] = Math.random();		// random color 0.0 <= B < 1.0
		}
	}
	// Repeat the 1st 2 vertices of the triangle strip to complete the torus:
			torVerts[j  ] = rbend + rbar;	// copy vertex zero;
						  //	x = (rbend + rbar*cos(phi==0)) * cos(theta==0)
			torVerts[j+1] = 0.0;
							//  y = (rbend + rbar*cos(phi==0)) * sin(theta==0) 
			torVerts[j+2] = 0.0;
							//  z = -rbar  *   sin(phi==0)
			torVerts[j+3] = 1.0;		// w
			torVerts[j+4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j+5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j+6] = Math.random();		// random color 0.0 <= B < 1.0
			j+=7; // go to next vertex:
			torVerts[j  ] = (rbend + rbar) * Math.cos(thetaStep);
						  //	x = (rbend + rbar*cos(phi==0)) * cos(theta==thetaStep)
			torVerts[j+1] = (rbend + rbar) * Math.sin(thetaStep);
							//  y = (rbend + rbar*cos(phi==0)) * sin(theta==thetaStep) 
			torVerts[j+2] = 0.0;
							//  z = -rbar  *   sin(phi==0)
			torVerts[j+3] = 1.0;		// w
			torVerts[j+4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j+5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j+6] = Math.random();		// random color 0.0 <= B < 1.0
}
function makeTetrahedron(){
	var c30 = Math.sqrt(0.75);					// == cos(30deg) == sqrt(3) / 2
	var sq2	= Math.sqrt(2.0);
tetVerts = new Float32Array([
  // Vertex coordinates(x,y,z,w) and color (R,G,B) for a color tetrahedron:
	//		Apex on +z axis; equilateral triangle base at z=0
/*	Nodes:
		 0.0,	 0.0, sq2, 1.0,			1.0, 	1.0,	1.0,	// Node 0 (apex, +z axis;  white)
     c30, -0.5, 0.0, 1.0, 		0.0,  0.0,  1.0, 	// Node 1 (base: lower rt; red)
     0.0,  1.0, 0.0, 1.0,  		1.0,  0.0,  0.0,	// Node 2 (base: +y axis;  grn)
    -c30, -0.5, 0.0, 1.0, 		0.0,  1.0,  0.0, 	// Node 3 (base:lower lft; blue)

*/
			// Face 0: (left side)  
     0.0,	 0.0, sq2, 1.0,			1.0, 	1.0,	1.0,	// Node 0
     c30, -0.5, 0.0, 1.0, 		0.0,  0.0,  1.0, 	// Node 1
     0.0,  1.0, 0.0, 1.0,  		1.0,  0.0,  0.0,	// Node 2
			// Face 1: (right side)
		 0.0,	 0.0, sq2, 1.0,			1.0, 	1.0,	1.0,	// Node 0
     0.0,  1.0, 0.0, 1.0,  		1.0,  0.0,  0.0,	// Node 2
    -c30, -0.5, 0.0, 1.0, 		0.0,  1.0,  0.0, 	// Node 3
    	// Face 2: (lower side)
		 0.0,	 0.0, sq2, 1.0,			1.0, 	1.0,	1.0,	// Node 0 
    -c30, -0.5, 0.0, 1.0, 		0.0,  1.0,  0.0, 	// Node 3
     c30, -0.5, 0.0, 1.0, 		0.0,  0.0,  1.0, 	// Node 1 
     	// Face 3: (base side)  
    -c30, -0.5,  0.0, 1.0, 		0.0,  1.0,  0.0, 	// Node 3
     0.0,  1.0,  0.0, 1.0,  	1.0,  0.0,  0.0,	// Node 2
     c30, -0.5,  0.0, 1.0, 		0.0,  0.0,  1.0, 	// Node 1
  ]);
}
function makeGroundGrid() {
//==============================================================================
// Create a list of vertices that create a large grid of lines in the x,y plane
// centered at x=y=z=0.  Draw this shape using the GL_LINES primitive.

	var xcount = 100;			// # of lines to draw in x,y to make the grid.
	var ycount = 100;		
	var xymax	= 50.0;			// grid size; extends to cover +/-xymax in x and y.
 	var xColr = new Float32Array([1.0, 1.0, 0.3]);	// bright yellow
 	var yColr = new Float32Array([0.5, 1.0, 0.5]);	// bright green.
 	
	// Create an (global) array to hold this ground-plane's vertices:
	gndVerts = new Float32Array(floatsPerVertex*2*(xcount+ycount));
						// draw a grid made of xcount+ycount lines; 2 vertices per line.
						
	var xgap = xymax/(xcount-1);		// HALF-spacing between lines in x,y;
	var ygap = xymax/(ycount-1);		// (why half? because v==(0line number/2))
	
	// First, step thru x values as we make vertical lines of constant-x:
	for(v=0, j=0; v<2*xcount; v++, j+= floatsPerVertex) {
		if(v%2==0) {	// put even-numbered vertices at (xnow, -xymax, 0)
			gndVerts[j  ] = -xymax + (v  )*xgap;	// x
			gndVerts[j+1] = -xymax;								// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		else {				// put odd-numbered vertices at (xnow, +xymax, 0).
			gndVerts[j  ] = -xymax + (v-1)*xgap;	// x
			gndVerts[j+1] = xymax;								// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		gndVerts[j+4] = xColr[0];			// red
		gndVerts[j+5] = xColr[1];			// grn
		gndVerts[j+6] = xColr[2];			// blu
	}
	// Second, step thru y values as wqe make horizontal lines of constant-y:
	// (don't re-initialize j--we're adding more vertices to the array)
	for(v=0; v<2*ycount; v++, j+= floatsPerVertex) {
		if(v%2==0) {		// put even-numbered vertices at (-xymax, ynow, 0)
			gndVerts[j  ] = -xymax;								// x
			gndVerts[j+1] = -xymax + (v  )*ygap;	// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		else {					// put odd-numbered vertices at (+xymax, ynow, 0).
			gndVerts[j  ] = xymax;								// x
			gndVerts[j+1] = -xymax + (v-1)*ygap;	// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		gndVerts[j+4] = yColr[0];			// red
		gndVerts[j+5] = yColr[1];			// grn
		gndVerts[j+6] = yColr[2];			// blu
	}
}
function makeSpaceship(){
ssVerts = new Float32Array ([
     0.00, 0.00, .25 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.50, 0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.25, 0.00, 0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.00, -0.10 ,1.00, Math.random(), Math.random(), Math.random(),
	 
	 0.00, -0.25, 0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 -0.25, 0.00, 0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.00, 0.25 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.50, 0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 
	 -0.25, 0.00, 0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.00, -0.10 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.25, 0.00,  0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, -0.25,  0.00 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.00,  0.25 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.25, 0.00,  0.00 ,1.00, Math.random(), Math.random(), Math.random(),
  ]);


}
function makeUFO(){
ufoVerts = new Float32Array ([
	 0.00, 0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 -1.00, 0.0, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 -0.81,0.0,-0.59,1.00, Math.random(), 0.4*Math.random(), Math.random(),//triangle 1 
	 0.00, -0.5, 0 ,1.00, 0.3, 0.5, 0.3,//triangle 2 
	 -0.31, 0.0, -0.95,1.00, Math.random(), Math.random(), Math.random(),//triangle 3
	 -0.81,0.0,-0.59,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),//triangle 4 
	 -0.31, 0.0, -0.95,1.00, Math.random(), Math.random(), Math.random(),//triangle 5
	 -0.81,0.0,-0.59,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, -0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),//triangle 6
	 0.31, 0.0, -0.95, 1.00, Math.random(), Math.random(), Math.random(),//triangle 7
     -0.81,0.0,-0.59,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),//triangle 8
	 0.31, 0.0, -0.95, 1.00, Math.random(), Math.random(), Math.random(),//triangle 9 
	 0.81, 0.0, -0.59, 1.00, Math.random(), Math.random(), Math.random(),//triangle 6
	 0.00, -0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),//triangle 6
	 0.81, 0.0, -0.59,1.00, Math.random(), Math.random(), Math.random(),//triangle 6
	 1.00, 0.0, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.81, 0, 0.59,1.0,Math.random(), Math.random(), Math.random(),
	 1.00, 0.0, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, -0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.31, 0.0, 0.95 ,1.00, Math.random(), Math.random(), Math.random(),
	 1.00, 0.0, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 -0.31, 0.0, 0.95 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.31, 0.0, 0.95 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, -0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 -0.81, 0.0, 0.59 ,1.00, Math.random(), Math.random(), Math.random(),
	 -0.31, 0.0, 0.95 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, 0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 -1.00, 0.0, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 -0.81, 0.0, 0.59 ,1.00, Math.random(), Math.random(), Math.random(),
	 0.00, -0.5, 0 ,1.00, Math.random(), Math.random(), Math.random(),
	 ])
}
function makeSolarSystem(){
sunVerts= makeSphere(1,1,0);
earthVerts= makeSphere(0,1,0);
moonVerts= makeSphere(1,1,1);
blueVerts= makeSphere(0,0,1);
}
function initVertexBuffer(gl) {
//==============================================================================
	makeTetrahedron();
	makeSpaceship();
	makeCube();
	makeTorus();
	makeSolarSystem();
	makeGroundGrid();	
	makeUFO();
	
	
	
	var mySiz = (tetVerts.length + ssVerts.length + cubVerts.length + torVerts.length +sphVerts.length + gndVerts.length + ufoVerts.length + sphVerts.length + sphVerts.length + sphVerts.length);
	var nn = mySiz / floatsPerVertex;
	console.log('nn is', nn, 'mySiz is', mySiz, 'floatsPerVertex is', floatsPerVertex);
	var colorShapes = new Float32Array(mySiz);
	tetStart=0 
	for(i=0,j=0; j< tetVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = tetVerts[j];
		}
	ssStart=i
	for(j=0; j< ssVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = ssVerts[j];
		}
	cubStart=i
	for(j=0; j< cubVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = cubVerts[j];
		}
	sunStart = i;						// next, we'll store the sphere;
	for(j=0; j< sphVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = sunVerts[j];
		}
		torStart = i;						// next, we'll store the torus;
	for(j=0; j< torVerts.length; i++, j++) {
		colorShapes[i] = torVerts[j];
		}
		gndStart = i;						// next we'll store the ground-plane;
	for(j=0; j< gndVerts.length; i++, j++) {
		colorShapes[i] = gndVerts[j];
		}
		ufoStart=i;
	for(j=0; j< ufoVerts.length; i++, j++) {
		colorShapes[i] = ufoVerts[j];
		}
		earthStart=i;
	for(j=0; j< sphVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = earthVerts[j];
		}
		moonStart=i;
	for(j=0; j< sphVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = moonVerts[j];
		}
		blueStart=i;
	for(j=0; j< sphVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = blueVerts[j];
		}
		
  // Create a buffer object
  var shapeBufferHandle = gl.createBuffer();  
  if (!shapeBufferHandle) {
    console.log('Failed to create the shape buffer object');
    return false;
  }

  // Bind the the buffer object to target:
  gl.bindBuffer(gl.ARRAY_BUFFER, shapeBufferHandle);
  // Transfer data from Javascript array colorShapes to Graphics system VBO
  // (Use sparingly--may be slow if you transfer large shapes stored in files)
  gl.bufferData(gl.ARRAY_BUFFER, colorShapes, gl.STATIC_DRAW);

  var FSIZE = colorShapes.BYTES_PER_ELEMENT; // how many bytes per stored value?
    
  //Get graphics system's handle for our Vertex Shader's position-input variable: 
  var a_Position = gl.getAttribLocation(gl.program, 'a_Position');
  if (a_Position < 0) {
    console.log('Failed to get the storage location of a_Position');
    return -1;
  }
  // Use handle to specify how to retrieve position data from our VBO:
  gl.vertexAttribPointer(
  		a_Position, 	// choose Vertex Shader attribute to fill with data
  		4, 						// how many values? 1,2,3 or 4.  (we're using x,y,z,w)
  		gl.FLOAT, 		// data type for each value: usually gl.FLOAT
  		false, 				// did we supply fixed-point data AND it needs normalizing?
  		FSIZE * 7, 		// Stride -- how many bytes used to store each vertex?
  									// (x,y,z,w, r,g,b) * bytes/value
  		0);						// Offset -- now many bytes from START of buffer to the
  									// value we will actually use?
  gl.enableVertexAttribArray(a_Position);  
  									// Enable assignment of vertex buffer object's position data

  // Get graphics system's handle for our Vertex Shader's color-input variable;
  var a_Color = gl.getAttribLocation(gl.program, 'a_Color');
  if(a_Color < 0) {
    console.log('Failed to get the storage location of a_Color');
    return -1;
  }
  
  
  // Use handle to specify how to retrieve color data from our VBO:
  gl.vertexAttribPointer(
  	a_Color, 				// choose Vertex Shader attribute to fill with data
  	3, 							// how many values? 1,2,3 or 4. (we're using R,G,B)
  	gl.FLOAT, 			// data type for each value: usually gl.FLOAT
  	false, 					// did we supply fixed-point data AND it needs normalizing?
  	FSIZE * 7, 			// Stride -- how many bytes used to store each vertex?
  									// (x,y,z,w, r,g,b) * bytes/value
  	FSIZE * 4);			// Offset -- how many bytes from START of buffer to the
  									// value we will actually use?  Need to skip over x,y,z,w
  									
  gl.enableVertexAttribArray(a_Color);  
  									// Enable assignment of vertex buffer object's position data

	//--------------------------------DONE!
  // Unbind the buffer object 
  gl.bindBuffer(gl.ARRAY_BUFFER, null);

  return nn;
}
function draw(gl, n, currentAngle, modelMatrix, u_ModelMatrix) {
//==============================================================================
  // Clear <canvas>  colors AND the depth buffer
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  
// Great question from student:
// "?How can I get the screen-clearing color (or any of the many other state
//		variables of OpenGL)?  'glGet()' doesn't seem to work..."
// ANSWER: from WebGL specification page: 
//							https://www.khronos.org/registry/webgl/specs/1.0/
//	search for 'glGet()' (ctrl-f) yields:
//  OpenGL's 'glGet()' becomes WebGL's 'getParameter()'

	clrColr = new Float32Array(4);
	clrColr = gl.getParameter(gl.COLOR_CLEAR_VALUE);
	console.log("clear value:", clrColr);
	

// NEXT, create different drawing axes, and...
	modelMatrix.setTranslate(0, 0.5, 0.0);  // 'set' means DISCARD old matrix,
	modelMatrix.scale(1,1,-1);							// convert to left-handed coord sys																			// to match WebGL display canvas.
	modelMatrix.scale(0.25, 0.25, 0.25);				// Make it smaller.
	var dist = Math.sqrt(xMdragTot*xMdragTot + yMdragTot*yMdragTot);
	modelMatrix.rotate(dist*120.0, -yMdragTot+0.0001, xMdragTot+0.0001, 0.0);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements); 
	//gl.drawArrays(gl.TRIANGLE_STRIP, tetStart/floatsPerVertex, tetVerts.length/floatsPerVertex);

  // Draw Space Ship that follows mouse
	modelMatrix.setTranslate(xpp  ,ypp  ,0.0);  
	modelMatrix.scale(1,1,-1);		
	modelMatrix.scale(0.2, 0.2, 0.2);
	 modelMatrix.rotate(currentAngle, 0, 1, 0);  // Spin on XY diagonal axis
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.drawArrays(gl.TRIANGLE_STRIP, ssStart/floatsPerVertex, ssVerts.length/floatsPerVertex);

//DRAW SUN
    modelMatrix.setTranslate( 0, 0, 0.0); // 'set' means DISCARD old matrix,
    modelMatrix.scale(1,1,-1);							// convert to left-handed coord sys
    modelMatrix.scale(0.2, 0.2, 0.2);
    modelMatrix.rotate(currentAngle, 0, 1, 0);  // Spin on XY diagonal axis
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
    gl.drawArrays(gl.TRIANGLE_STRIP,sunStart/floatsPerVertex,sphVerts.length/floatsPerVertex);
	    
		
	pushMatrix(modelMatrix);


	//Draw Planet
	modelMatrix.translate( 4,0, 0, 0.0); // 'set' means DISCARD old matrix,
    modelMatrix.scale(0.5, 0.5, .5);
    modelMatrix.rotate(currentAngle*2, 0, 0, 1);  // Spin on XY diagonal axis
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
    gl.drawArrays(gl.TRIANGLE_STRIP,earthStart/floatsPerVertex,sphVerts.length/floatsPerVertex);
   
   // DRAW RINGS
    modelMatrix.translate(0.0, 0, 0.0);				
    modelMatrix.scale(1,1,1);							
    modelMatrix.rotate(currentAngle, 0, 1, 1);  
    gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
    gl.drawArrays(gl.TRIANGLES, torStart/floatsPerVertex,torVerts.length/floatsPerVertex);	

	// Draw Moon
	modelMatrix.rotate(currentAngle*4, 1,0, 1);  // Spin on XY diagonal axis
	modelMatrix.translate( 0,2, 0, 0.0); // 'set' means DISCARD old matrix,
    modelMatrix.scale(.25,.25,.25);							
    modelMatrix.translate( 0,-2, 0, 0.0);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
    gl.drawArrays(gl.TRIANGLE_STRIP,moonStart/floatsPerVertex,sphVerts.length/floatsPerVertex);
   	
	
	//back to sun axes 
	modelMatrix = popMatrix();
	pushMatrix(modelMatrix);
	//Draw Planet 2 
	modelMatrix.rotate(currentAngle, 1, 1, 1);  // Spin on XY diagonal axis
	modelMatrix.translate( 0,2, 0, 0.0); // 'set' means DISCARD old matrix,
    modelMatrix.scale(0.35, 0.35, .35);
    modelMatrix.rotate(currentAngle*5, 0, 0, 1);  // Spin on XY diagonal axis
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
    gl.drawArrays(gl.TRIANGLE_STRIP,blueStart/floatsPerVertex,sphVerts.length/floatsPerVertex);
	
	
	modelMatrix.translate( 0,2, 0., 0.0)
	modelMatrix.scale(0.5, 0.5, .5);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
    gl.drawArrays(gl.TRIANGLE_STRIP, ufoStart/floatsPerVertex, ufoVerts.length/floatsPerVertex);	
	
	modelMatrix = popMatrix();


//Draw UFO 
    modelMatrix.setTranslate( 0, 0.7, 0.0);	
	modelMatrix.scale(0.25, 0.25, 0.25);				// shrink by 10X:
	var dist = Math.sqrt(xMdragTot*xMdragTot + yMdragTot*yMdragTot);
	modelMatrix.rotate(dist*120.0, -yMdragTot+0.0001, xMdragTot+0.0001, 0.0);
// Pass our current matrix to the vertex shaders:
  gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
  // Draw just the ground-plane's vertices
  gl.drawArrays(gl.TRIANGLE_STRIP, ufoStart/floatsPerVertex, ufoVerts.length/floatsPerVertex);	
							  
}

// Last time that this function was called:  (used for animation timing)
var g_last = Date.now();

function animate(angle) {
//==============================================================================
  // Calculate the elapsed time
  var now = Date.now();
  var elapsed = now - g_last;
  g_last = now;
  
  // Update the current rotation angle (adjusted by the elapsed time)
  //  limit the angle to move smoothly between +20 and -85 degrees:
//  if(angle >  120.0 && ANGLE_STEP > 0) ANGLE_STEP = -ANGLE_STEP;
//  if(angle < -120.0 && ANGLE_STEP < 0) ANGLE_STEP = -ANGLE_STEP;
  
  var newAngle = angle + (ANGLE_STEP * elapsed) / 1000.0;
  if(newAngle > 180.0) newAngle = newAngle - 360.0;
  if(newAngle <-180.0) newAngle = newAngle + 360.0;
  return newAngle;
}

//==================HTML Button Callbacks======================

function angleSubmit() {
// Called when user presses 'Submit' button on our webpage
//		HOW? Look in HTML file (e.g. ControlMulti.html) to find
//	the HTML 'input' element with id='usrAngle'.  Within that
//	element you'll find a 'button' element that calls this fcn.

// Read HTML edit-box contents:
	var UsrTxt=document.getElementById('usrAngle').value;	
// Display what we read from the edit-box: use it to fill up
// the HTML 'div' element with id='Result':
  document.getElementById('Result').innerHTML ='You Typed: '+UsrTxt;
};

function clearDrag() {
// Called when user presses 'Clear' button in our webpage
	xMdragTot = 0.0;
	yMdragTot = 0.0;
}

function spinUp() {
// Called when user presses the 'Spin >>' button on our webpage.
// ?HOW? Look in the HTML file (e.g. ControlMulti.html) to find
// the HTML 'button' element with onclick='spinUp()'.
  ANGLE_STEP += 25; 
}

function spinDown() {
// Called when user presses the 'Spin <<' button
 ANGLE_STEP -= 25; 
}

function runStop() {
// Called when user presses the 'Run/Stop' button
  if(ANGLE_STEP*ANGLE_STEP > 1) {
    myTmp = ANGLE_STEP;
    ANGLE_STEP = 0;
  }
  else {
  	ANGLE_STEP = myTmp;
  }
}

//===================Mouse and Keyboard event-handling Callbacks

function myMouseDown(ev, gl, canvas) {
//==============================================================================
// Called when user PRESSES down any mouse button;
// 									(Which button?    console.log('ev.button='+ev.button);   )
// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
  var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
  var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
  var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
//  console.log('myMouseDown(pixel coords): xp,yp=\t',xp,',\t',yp);
	// Convert to Canonical View Volume (CVV) coordinates too:
    x = (xp - canvas.width/2)  / 		// move origin to center of canvas and
  						 (canvas.width/2);			// normalize canvas to -1 <= x < +1,
	y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
							 (canvas.height/2);
//	console.log('myMouseDown(CVV coords  ):  x, y=\t',x,',\t',y);
	
	if (!follow) {follow = true;}
	
	else if (follow) {follow = false;}
    
  
	isDrag = true;											// set our mouse-dragging flag
	xMclik = x;												// record where mouse-dragging began
	yMclik = y;
};



function myMouseMove(ev, gl, canvas) {
//==============================================================================
// Called when user MOVES the mouse with a button already pressed down.
// 									(Which button?   console.log('ev.button='+ev.button);    )
// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

				// IGNORE all mouse-moves except 'dragging'

	// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
  var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
  var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
	var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
//  console.log('myMouseMove(pixel coords): xp,yp=\t',xp,',\t',yp);
  
	// Convert to Canonical View Volume (CVV) coordinates too:
  var x = (xp - canvas.width/2)  / 		// move origin to center of canvas and
  						 (canvas.width/2);			// normalize canvas to -1 <= x < +1,
	var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
							 (canvas.height/2);
//	console.log('myMouseMove(CVV coords  ):  x, y=\t',x,',\t',y);

	if(follow){
		xpp=x;
		ypp=y;}
	if(isDrag==false) return;	
	// find how far we dragged the mouse:
	xMdragTot += (x - xMclik);					// Accumulate change-in-mouse-position,&
	yMdragTot += (y - yMclik);
	xMclik = x;													// Make next drag-measurement from here.
	yMclik = y;
	
	};

function myMouseUp(ev, gl, canvas) {
//==============================================================================
// Called when user RELEASES mouse button pressed previously.
// 									(Which button?   console.log('ev.button='+ev.button);    )
// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
  var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
  var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
	var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
//  console.log('myMouseUp  (pixel coords): xp,yp=\t',xp,',\t',yp);
  
	// Convert to Canonical View Volume (CVV) coordinates too:
  var x = (xp - canvas.width/2)  / 		// move origin to center of canvas and
  						 (canvas.width/2);			// normalize canvas to -1 <= x < +1,
	var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
							 (canvas.height/2);
	console.log('myMouseUp  (CVV coords  ):  x, y=\t',x,',\t',y);
	
	isDrag = false;											// CLEAR our mouse-dragging flag, and
	// accumulate any final bit of mouse-dragging we did:
	xMdragTot += (x - xMclik);
	yMdragTot += (y - yMclik);
	console.log('myMouseUp: xMdragTot,yMdragTot =',xMdragTot,',\t',yMdragTot);
};

function myKeyDown(ev) {
//===============================================================================
// Called when user presses down ANY key on the keyboard, and captures the 
// keyboard's scancode or keycode(varies for different countries and alphabets).
//  CAUTION: You may wish to avoid 'keydown' and 'keyup' events: if you DON'T 
// need to sense non-ASCII keys (arrow keys, function keys, pgUp, pgDn, Ins, 
// Del, etc), then just use the 'keypress' event instead.
//	 The 'keypress' event captures the combined effects of alphanumeric keys and // the SHIFT, ALT, and CTRL modifiers.  It translates pressed keys into ordinary
// ASCII codes; you'll get the ASCII code for uppercase 'S' if you hold shift 
// and press the 's' key.
// For a light, easy explanation of keyboard events in JavaScript,
// see:    http://www.kirupa.com/html5/keyboard_events_in_javascript.htm
// For a thorough explanation of the messy way JavaScript handles keyboard events
// see:    http://javascript.info/tutorial/keyboard-events
//

	switch(ev.keyCode) {			// keycodes !=ASCII, but are very consistent for 
	//	nearly all non-alphanumeric keys for nearly all keyboards in all countries.
		case 37:		// left-arrow key
			spinDown();
			break;
		case 38:		// up-arrow key
			
			break;
		case 39:		// right-arrow key
					spinUp();
  		break;
		case 40:		// down-arrow key
			
  		break;
		case 32:
			runStop()
			break;
		case 113: //F2
		document.getElementById('Help').innerHTML= 
			'<br>Controls: Left Arrow - Decrease Speed/Reverse</br><br>          Right Arrow - Increase Speed</br><br>          Space Bar - Start/Stop </br><br>          Mouse Click - Toggle Spaceship Follow Mouse </br><br>          Click and Drag - Rotate UFO </br>';
		break
		default:
			console.log('myKeyDown()--keycode=', ev.keyCode, ', charCode=', ev.charCode);
  		document.getElementById('Result').innerHTML =
  			'myKeyDown()--keyCode='+ev.keyCode;
			break;
	}
}

function myKeyUp(ev) {
//===============================================================================
// Called when user releases ANY key on the keyboard; captures scancodes well

	console.log('myKeyUp()--keyCode='+ev.keyCode+' released.');
}


function myKeyPress(ev) {
//===============================================================================
// Best for capturing alphanumeric keys and key-combinations such as 
// CTRL-C, alt-F, SHIFT-4, etc.

	console.log('myKeyPress():keyCode='+ev.keyCode  +', charCode=' +ev.charCode+
												', shift='    +ev.shiftKey + ', ctrl='    +ev.ctrlKey +
												', altKey='   +ev.altKey   +
												', metaKey(Command key or Windows key)='+ev.metaKey);
												
												
}