/* DRAWING BOARD HANDLERS */
// Top level function for tracking mouse movement and clicks
function main() {
    document.onmousemove = function(e){
        var pos = getXY(canvas, e);
    		cursorX = pos.x;
    		cursorY = pos.y;
    };
    canvas.addEventListener("mousewheel", MouseWheelHandler, false);
    canvas.addEventListener("click", checkNode, false);
}

function panBoard(pan, e) {
  if (nodeClicked !== 1) {
    // Update existing drawings
    // redraw any edits
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    dtx.clearRect(0, 0, canvas.width, canvas.height);
    translateBy(circlesArr, 0, pan);
    refreshLines(ctx);
    drawPaths(edges);
    drawCircles(circlesArr);
  }
}

function MouseWheelHandler(e) {
// https://stackoverflow.com/questions/25909341/disable-mouse-wheel-to-scroll-up-or-down
    // cross-browser wheel delta
    var e = window.event || e; // old IE support
    var delta = Math.max(-1, Math.min(1, (e.wheelDelta || -e.detail)));

    if(delta==1)         // if mouse scrolls up
    {
        e.preventDefault();
        e.stopPropagation();
        panBoard(20, e);
        return false;
    }
    if(delta==-1)        // if mouse scrolls down, we disable scrolling.
    {
        e.preventDefault();
        e.stopPropagation();
        panBoard(-20, e);
        return false;
    }
    return false;
}

function checkNode(e) {
// main handler for mouse clicks
	var pos = getXY(canvas, e);
	clickPos = [pos.x, pos.y];
  //start the loop
  for (i = 0; i < circlesArr.length; i++) {
      // check if mouse clicks on any nodes
      if (isHit(circlesArr[i], pos.x, pos.y)) {
          if (nodeClicked == 0) {
          // if no nodes are already selected, select the clicked node
              highlightNode(circlesArr[i]);
              nodeClicked = 1;
              lastClicked = circlesArr[i];
          }
          else if (nodeClicked == 1) {
          // if a node is already selected, draw a line to the clicked node
              drawLineBtwnCircles(lastClicked, circlesArr[i]);
              unlightNode(lastClicked);
              nodeClicked = 0;
              lastClicked = null;
          }
          // stop searching the list if a node was clicked
          break;
      }
  }
  // redraw any edits
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  refreshLines(ctx);
}

// Helper function to drawLineBtwnCircles to know which two circles to draw the center-to-center line
// Determines which circles the current line spans
function whichCircles(currLine) {
	for (i = 0; i < circlesArr.length; i++) {
		if (isHit(circlesArr[i], currLine[0][0], currLine[0][1])) {
			var hit1 = circlesArr[i];
		}
		if (isHit(circlesArr[i], currLine[1][0], currLine[1][1])) {
			var hit2 = circlesArr[i];
		}
	}
	return [hit1, hit2];
}

// Draws a center-to-center line between two circles
function drawLineBtwnCircles(circle1, circle2) {
	pt1 = [circle1.x, circle1.y];
	pt2 = [circle2.x, circle2.y];
	// Remove if the line already exists
	if (lineExists(circle1, circle2)) {
		linesRemove([pt1, pt2]);
	}
	// Add if it doesn't
	else {
		lines.push([pt1, pt2]);
		connections.push([circle1.index, circle2.index]);
	}
	console.log(connections);
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	refreshLines(ctx);
}

// Helper function to lines array to avoid adding duplicates
function lineExists(circle1, circle2) {
	var pt1 = [circle1.x, circle1.y];
	var pt2 = [circle2.x, circle2.y];

	for (i = 0; i < lines.length; i++) {
		if (linesEqual(lines[i], [pt1, pt2])) {
			return true;
		}
	}
	return false;
}

// Helper function to lines array for removing an element
function linesRemove(line) {
	for (i = 0; i < lines.length; i++) {
		if (linesEqual(line, lines[i])) {
			lines.splice(i, 1);
			connections.splice(i, 1);
		}
	}
}

// Helper function to lines array for checking equality
function linesEqual(line1, line2) {
	var line1x1 = line1[0][0];
	var line1y1 = line1[0][1];
	var line1x2 = line1[1][0];
	var line1y2 = line1[1][1];
	var line2x1 = line2[0][0];
	var line2y1 = line2[0][1];
	var line2x2 = line2[1][0];
	var line2y2 = line2[1][1];
	if ((line1x1 == line2x1 && line1y1 == line2y1 && line1x2 == line2x2 && line1y2 == line2y2) || (line1x1 == line2x2 && line1y1 == line2y2 && line1x2 == line2x1 && line1y2 == line2y1)) {
		return true;
	}
return false;
}

// Circle object
var Circle = function(x, y, radius, index) {
	this.x = x;
	this.y = y;
	this.radius = radius;
	this.index = index;

	this.render = function(ctx, color) {
		ctx.beginPath();
		ctx.arc(this.x, this.y, this.radius, 0, 2 * Math.PI, false);
		ctx.fillStyle = color;
		ctx.fill();
		ctx.lineWidth = 3;
		ctx.strokeStyle = '#003300';
		ctx.stroke();
	}
}

// Translates an array of circles
function translateBy(c_arr, x, y) {
	for (i = 0; i < c_arr.length; i++) {
		circle = c_arr[i];
		circle.x += x;
		circle.y += y;
	}
  for (i = 0; i < lines.length; i++) {
    lines[i][0][1] += y;
    lines[i][1][1] += y;
  }
}

// Renders an array of circle objects
function drawCircles(arr) {
	for (let i = 0; i < arr.length; i++) {
		this_circle = arr[i];
		this_circle.render(dtx, '#c7e67f');
	}
}

// Draws a list of edges between nodes
function drawPaths(arr) {
	for (let i = 0; i < arr.length; i++) {
			idx1 = arr[i][0];
			idx2 = arr[i][1];
			pt1 = [circlesArr[idx1].x, circlesArr[idx1].y];
			pt2 = [circlesArr[idx2].x, circlesArr[idx2].y];
			dtx.globalAlpha = 0.5;
			dtx.beginPath();
			dtx.moveTo(...pt1);
			dtx.lineTo(...pt2);
			dtx.lineWidth = 6;
			dtx.strokeStyle = '#800000';
			dtx.stroke();
	}
	dtx.globalAlpha = 1.0;
}

// See https://codeboxsystems.com/tutorials/en/how-to-drag-and-drop-objects-javascript-canvas/
// Helper function to check if any object intersects with circle object
function isHit(shape, x, y) {
  if (shape.constructor.name === 'Circle') {
    var dx = shape.x - x;
    var dy = shape.y - y;
    if (dx * dx + dy * dy < shape.radius * shape.radius) {
      return true
    }
  }
  return false;
}

// Redraws all saved lines
function refreshLines(ctx) {
	for (let i = 0; i < lines.length; i++) {
    // Draws the line
		ctx.beginPath();
    ctx.moveTo(...lines[i][0]);
    ctx.lineTo(...lines[i][1]);
		ctx.lineWidth = 3;
    ctx.strokeStyle = '#000000';
    ctx.stroke();
    // Draws the joint position
    ctx.beginPath();
		ctx.arc(...lines[i][0], 3, 0, 2 * Math.PI, false);
		ctx.fillStyle = '#f5bf42';
		ctx.fill();
		ctx.lineWidth = 2;
		ctx.strokeStyle = '#000000';
		ctx.stroke();

    ctx.beginPath();
		ctx.arc(...lines[i][1], 3, 0, 2 * Math.PI, false);
		ctx.fillStyle = '#f5bf42';
		ctx.fill();
		ctx.lineWidth = 2;
		ctx.strokeStyle = '#000000';
		ctx.stroke();
	}
}

// From: https://stackoverflow.com/questions/29501447/why-does-css-centering-mess-up-canvas-mouse-coordinates
function getXY(canvas, event) {
    var rect = canvas.getBoundingClientRect();  // absolute position of canvas
    return {
        x: event.clientX - rect.left,
        y: event.clientY - rect.top
    }
}

/* USER INTERFACE */
// Button function for removing the last drawn line
function undoLines() {
	lines.pop();
	connections.pop();
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	refreshLines(ctx);
}

// Button function for clearing all the saved lines
function clearLines() {
	lines = [];
	connections = [];
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	refreshLines(ctx);
}

function highlightNode(circle) {
    circle.render(dtx, 'yellow');
}

function unlightNode(circle) {
    circle.render(dtx, '#c7e67f');
}

function setError(m) {
	document.getElementById("error-text").innerHTML = m;
}

/* PRE-PROCESSING */
// Helper function for converting the input coordinates to the displayed nodes on canvas "diagram"
function convertCoordsToCircles(coords, radius){
	const circlesArr = [];
	for (let i = 0; i < coords.length; i++) {
		var circleX = coords[i][0],
			circleY = coords[i][1];
		var newCircle = new Circle(circleX, circleY, radius, i);
		circlesArr.push(newCircle);
	}
	return circlesArr;
}

// Scales the coordinates to the canvas dimensions
function evalCoordsBox(coords) {
	var coordsX = [];
	var coordsY = [];
	for (let i = 0; i < coords.length; i++) {
		coordsX.push(coords[i][0]);
		coordsY.push(coords[i][1]);
	}
	let padding = 1;
	var xlim = [Math.min(...coordsX) - padding, Math.max(...coordsX) + padding];
	var ylim = [Math.min(...coordsY) - padding, Math.max(...coordsY) + padding];
	return [xlim, ylim];
}

// Scales drawn circles to canvas dimensions
function fitCoordsToCanvas(coords, board) {
  	const newCoords = [];

    // Box the coords
    // This finds the x and y maximum and minimum points
    var box = evalCoordsBox(coords);
    var minx = box[0][0],
        maxx = box[0][1],
        miny = box[1][0],
        maxy = box[1][1];
    // Convert to range for the WxH dimensions
    var rangex = maxx - minx,
        rangey = maxy - miny;

    // Hardcoded for scaling 0-600bps to 800px canvas
    let pctScale = 25;

    // Find the final pixel dimensions of the box
    // For some reason X needs to be spaced a little bit more otherwise
    // some gaps are too tight
    rangex = Math.round(rangex * pctScale);
    rangey = Math.round(rangey * pctScale);

    // Determine the translated coordinates to center the box
    let mdptx = rangex / 2;
    let mdpty = rangey / 2;
    var tx = board.width/2 - mdptx;
    var ty = board.height/2 - mdpty;

    // Calculate the scaled coordinates
  	const pctCoords = coordsToPct(coords);
    var boxheight = rangey;
    var boxwidth = rangex;
    for (let i = 0; i < coords.length; i++) {
  		var nx = pctCoords[i][0] * boxwidth,
  			ny = pctCoords[i][1] * boxheight;
  		newCoords.push([nx, ny]);
  	}
  	return [newCoords, tx, ty];
}

// Scales the coordinates to the canvas dimensions
function coordsToPct(coords) {
  	const pctCoords = [];
  	var box = evalCoordsBox(coords);
    var minx = box[0][0],
        maxx = box[0][1],
        miny = box[1][0],
        maxy = box[1][1];
	for (let i = 0; i < coords.length; i++) {
		var px = linPct(minx, maxx, coords[i][0]),
		    py = 1.0 - linPct(miny, maxy, coords[i][1]);
		pctCoords.push([px, py])
	}
	return pctCoords;
}

// Helper math function for converting min/max to percentage dimensions
function linPct(min, max, pt) {
	return (pt - min)/(max - min);
}

// Helper function for scaling radius of circle objects
function scaleCircleRadiusToCanvas(board, coords) {
	area = board.width * board.height;
  var box = evalCoordsBox(coords);
    // if x range > y range
    var minx = box[0][0],
        maxx = box[0][1],
        miny = box[1][0],
        maxy = box[1][1];
    var rangex = maxx - minx,
        rangey = maxy - miny;
	return area / 20000;
}

// Sets the line array according to a preset connections array
function applyConnections(vars) {
    connections = vars;
    for (let i = 0; i < vars.length; i++) {
        line = vars[i];
        idx1 = line[0];
        idx2 = line[1];
        c1 = circlesArr[idx1];
        c2 = circlesArr[idx2];
        lines.push([[c1.x, c1.y],[c2.x, c2.y]]);
    }
    refreshLines(ctx);
}

/* GLOBALS */

// Gets the canvas elements
var diagram = document.getElementById("diagram");
var dtx = diagram.getContext("2d");
var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");

// Initializes mouse position
var cursorX;
var cursorY;
var firstClick = [0,0];
var intervalLoop = null;

// Initialize data variables
var lines = [];  // Array for lines. Implemented as stack.
var connections = [];  // Array for output

// Initialize UI variables
var currLine = [];  // XY coordinates for two points of line
var nodeClicked = 0;
var lastClicked = null;
var ypan = 0;

main(); // Begins to monitor mouse movement
