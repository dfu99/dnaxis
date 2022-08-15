import * as THREE from "./three.min.js";
import * as PREVIEW from "./preview.js";
// Top level function for tracking mouse movement and clicks
function main() {
    document.onmousemove = function(e){
        bdraw.clearRect(0, 0, board_draw.width, board_draw.height);
        var pos = getXY(board_draw, e);

    		cursorX = pos.x;
    		cursorY = pos.y;
        if (inBounds(board_draw, e) && nodeClicked != 1) {
          // Prints position of mouse
          printSnap(cursorX, cursorY, XOVERMODE);
        }
        // if coordinates empty, place the first node anywhere
        if (inBounds(board_draw, e) && coords.length == 0) {
          newCircle.setPos(cursorX, cursorY);
          newCircle.render(bdraw, "white");

          // Prints position of new node
          printSnap(...newCircle.getPos(), XOVERMODE);
        }
        // if nodes exist, only place nodes relative to existing nodes
        if (nodeClicked == 1) {
            // Draw a new node that follow the mouse direction
            // but stays within a certain distance of the last node
            var [px, py] = lastClicked.getPos();

            let vector = [cursorX - px, cursorY - py];
            let nx = vector[0] / math.norm(vector),
                ny = vector[1] / math.norm(vector);
            let xx = px + nx * INTERHELICAL * 10 * uiscale,
                yy = py + ny * INTERHELICAL * 10 * uiscale;
            newCircle.setPos(xx, yy);
            newCircle.render(bdraw, "white");

            // Prints position of new node
            printSnap(...newCircle.getPos(), XOVERMODE);
        }
    };

    board_draw.addEventListener("mousewheel", MouseWheelHandler, false);
    board_draw.addEventListener("click", checkNode, false);
    board_draw.addEventListener("contextmenu", RMBHandler, false);
}

function panBoard(pan, e) {
  if (nodeClicked !== 1) {
    ypan += pan;
    // Update axes
    baxes.clearRect(0, 0, board_axes.width, board_axes.height);
    drawAxes();

    // Update the coordinates display
    var pos = getXY(board_draw, e);
    cursorX = pos.x;
    cursorY = pos.y;
    if (inBounds(board_draw, e) && nodeClicked != 1) {
      // Prints position of mouse
      printSnap(cursorX, cursorY, XOVERMODE);
    }

    // Update any existing Nodes
    updateDraw();
  }
}

function RMBHandler(e) {
  var isRightMB;
  e = e || window.event;
  e.preventDefault();
  var pos = getXY(board_draw, e);
	var clickPos = [pos.x, pos.y];
  for (let i=0; i<circlesArr.length; i++) {
    if (isHit(circlesArr[i], ...clickPos)) {
      circlesArr[i].flip(bdisp);
      updateCoords();
    }
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
        panBoard(10, e);
        return false;
    }
    if(delta==-1)        // if mouse scrolls down, we disable scrolling.
    {
        e.preventDefault();
        e.stopPropagation();
        panBoard(-10, e);
        return false;
    }
    return false;
}

// Updates the saved coordinates and the preview when the text box inputs
// is changed
$('#mesh_txt_input').on('input propertychange paste', function() {
    var text = document.getElementById("mesh_txt_input").value;
    const lines = text.split(/\r?\n/);
    coords = [];
    for (let i=0; i<lines.length; i++) {
      var values = lines[i].split(/,/);
      // Enforce input format
        if (!(values.length == 3))
          { console.log("Exited for length"); continue; }
        if (!(Number.isInteger(parseFloat(values[0]))))
          { console.log("bp is invalid "); continue; }
        if (!(Number.isFinite(parseFloat(values[1]))))
          { console.log("height is invalid"); continue; }
        if (!(intIsBool(parseInt(values[2]))))
          { console.log("dirBit is invalid"); continue; }

          let bp = 0.332 * parseFloat(values[0]) / (2 * Math.PI);
          let h = parseFloat(values[1]);
          let db = parseFloat(values[2]);
          coords.push([bp, h, db]);
    }
    updateDraw();
    PREVIEW.preview(coords);
});

function intIsBool(v) {
    if (v === 1 || v === 0) return true;
    else return false;
}

// Update the text area with the coordinates
function updateText() {
    var coordsString = "";
    for (let i=0; i<coords.length; i++) {
        let x = coords[i][0],
            y = coords[i][1],
            dirBit = coords[i][2];
            var bp = radius2bp(x * 10);
            let cf = XOVERMODE;
            var round = bp % cf;
            bp = round > cf/2 ? bp - (bp % cf) + cf : bp - (bp % cf);
            var height = Math.round(y * 10)/10;
        coordsString = coordsString + String(bp) + "," + String(height) + "," + String(dirBit) + "\n";
    }
    document.getElementById("mesh_txt_input").value = coordsString;
}

// Update the node diagram with coordinates
function updateDraw() {
  idx = 0;
  circlesArr = [];
  bdisp.clearRect(0, 0, board_display.width, board_display.height);
  for (let i=0; i<coords.length; i++) {
      var x = coords[i][0] * 10,
          y = coords[i][1] * 10,
          db = coords[i][2];
      [x, y] = cart2canvas(x, y, board_display);
      var refreshCircle = new Circle(x, y, circleRadius, idx, db);
      circlesArr.push(refreshCircle);
      refreshCircle.render(bdisp);
  }
}

function updateCoords() {
  coords = [];
  for (let i=0; i<circlesArr.length; i++) {
    var x = circlesArr[i].x,
    y = circlesArr[i].y,
    db = circlesArr[i].dirB;
    coords.push([x/10, y/10, db]);
  }
  PREVIEW.preview(coords);
  updateText();
}

function checkNode(e) {
  	var pos = getXY(board_draw, e);
  	var clickPos = [pos.x, pos.y];
    if (nodeClicked == 0) {
        // Place new node if canvas is blank
        if (circlesArr.length == 0) {
            var [xx, yy] = snapToFactor(...clickPos, XOVERMODE);
            var saveCircle = new Circle(xx, yy, circleRadius, idx, dirBit);
            idx = idx + 1;
            saveCircle.render(bdisp);
            circlesArr.push(saveCircle);
            console.log(circlesArr);
            var x, y;
            [x, y] = canvas2cart(xx, yy, board_display);
            x = x/10;
            y = y/10;

            bdraw.clearRect(0, 0, board_draw.width, board_draw.height);
            coords.push([x, y, dirBit]);
            PREVIEW.preview(coords);

        }
        else {
            // Selects a node if not placing first node
            for (let i=0; i<circlesArr.length; i++) {

                if (isHit(circlesArr[i], ...clickPos)) {
                    nodeClicked = 1;
                    lastClicked = circlesArr[i];
                    circlesArr[i].highlight(bdisp);
                }
            }
        }
    }
    else if (nodeClicked == 1) {
      // Place the next node from a selected node
        // Get the snapped position
        var [px, py] = lastClicked.getPos();
        var vector = [pos.x - px, pos.y - py];
        var nx = vector[0] / math.norm(vector),
        ny = vector[1] / math.norm(vector);
        var xx = px + nx * INTERHELICAL * 10 * uiscale,
            yy = py + ny * INTERHELICAL * 10 * uiscale;
        [xx, yy] = snapToFactor(xx, yy, XOVERMODE);
        if (collides(xx, yy)) {
          // error if nodes overlap
            printError("ERROR: Helices cannot overlap. Try a different placement.");
        }
        else {
          // good placement
          // flip dirBit from selected node
          dirBit = lastClicked.dirB ^ 1;
          var saveCircle = new Circle(xx, yy, circleRadius, idx, dirBit);
          nodeClicked = 0;
          idx = idx + 1;
          lastClicked.unlight(bdisp);
          lastClicked = null;
          saveCircle.render(bdisp);
          circlesArr.push(saveCircle);
          console.log(circlesArr);
          var x, y;
          [x, y] = canvas2cart(xx, yy, board_display);
          x = x/10;
          y = y/10;
          bdraw.clearRect(0, 0, board_draw.width, board_draw.height);
          coords.push([x, y, dirBit]);
          PREVIEW.preview(coords);
          // After placement, check if we are over the scaf len limit
          checkLength();
          updateText();
        }
    }
}

function collides(x, y) {
  for (let i=0; i<circlesArr.length; i++) {
    var [xx, yy] = circlesArr[i].getPos();
    var vector = [x - xx, y - yy];
    var distance = math.norm(vector);
    if (distance < 25 * uiscale) {
      return true;
    }
  }
  return false;
}

function snapToFactor(x, y, m) {
    var xx, yy;
    [xx, yy] = canvas2cart(x, y, board_display);
    // Round the x-coordinate
    var bp = radius2bp(xx);
    var round = bp % m;
    var bpr = round > m/2 ? bp - (bp % m) + m : bp - (bp % m);
    xx = bp2radius(bpr);

    // Round the y-coordinate
    var yyr = Math.round(yy * 10) / 10;

    // Convert back to canvas coordinates
    [xx, yy] = cart2canvas(xx, yyr, board_display);
    return [xx, yy];
}

async function printError(s) {
    bmsg.fillStyle="red";
    bmsg.font = "16px Helvetica";
    bmsg.fillText(s, 5, 50);
    await new Promise(r => setTimeout(r, 5000));
    bmsg.clearRect(0, 0, board_xy.width, board_xy.height);
}

function printSnap(x, y, m) {
  // Prints the coordinates on the top right corner of the canvas
  // Converts (x, y) in (bp, height)
  // Rounds bp to m

  // Enforce integer inputs
  x = Math.round(x);
  y = Math.round(y);
  // Constantly refresh the board
  bxy.clearRect(0, 0, board_xy.width, board_xy.height);
  binfo.clearRect(0, 0, board_xy.width, board_xy.height);


  var xx, yy;
  [xx, yy] = canvas2cart(x, y, board_display);
  // Round the x-coordinate
  var bp = radius2bp(xx);
  var round = bp % m;
  var bpr = round > m/2 ? bp - (bp % m) + m : bp - (bp % m);

  // Round the y-coordinate
  var yyr = Math.round(yy * 10) / 100;

  bxy.font = "15px Helvetica";
  bxy.fillStyle = "black";
  bxy.fillText("Base pairs:"+String(bpr), 5, 20);
  bxy.fillText("Height:"+String(yyr), 125, 20);
  bxy.fillText("Est. Scaffold Length: "+measureScafLength(), 225, 20);
}

// See https://codeboxsystems.com/tutorials/en/how-to-drag-and-drop-objects-javascript-canvas/
// Helper function to check if any object intersects with circle object
function isHit(shape, x, y) {
  if (shape.constructor.name === 'Circle') {
    var [sx, sy] = shape.getPos();
    var dx = sx - x;
    var dy = sy - y;
    if (dx * dx + dy * dy < shape.radius * shape.radius) {
      return true
    }
  }
  return false;
}

// Circle object
var Circle = function(cursorX, cursorY, radius, index, dirB, board) {

    this.setPos = function(cursorX, cursorY) {
      [this.x, this.y] = canvas2cart(cursorX, cursorY, board_display);
    }

    this.getPos = function() {
      return cart2canvas(this.x, this.y, board_display);
    }

  // inputs are cursor positions on canvas
  // save attributes in cartesian
  this.setPos(cursorX, cursorY);
	this.radius = radius;
	this.index = index;
  this.dirB = dirB;
  this.selected = false;
  this.setColor = function() {
    if (this.selected == true) {
      this.color = '#F0E442';
    }
    else if (this.dirB == 1) {
      this.color = FIVECOLOR;
    }
    else if (this.dirB == 0) {
      this.color = THREECOLOR;
    }
  }
  this.setColor();

	this.render = function(ctx, color=this.color) {
		ctx.beginPath();
    var [xx, yy] = this.getPos();
		ctx.arc(xx, yy, this.radius, 0, 2 * Math.PI, false);
		ctx.fillStyle = color;
		ctx.fill();

    if (this.selected == false) {
		    ctx.strokeStyle = '#000000';
        ctx.lineWidth = 3;
    }
    else {
      ctx.strokeStyle = '#ff0000';
      ctx.lineWidth = 1;
    }
		ctx.stroke();
	}

  this.flip = function(ctx) {
    this.dirB ^= 1;
    this.setColor();
    this.render(ctx);
  }

  this.highlight = function(ctx) {
    this.selected = true;
    this.setColor();
    this.render(ctx);
  }
  this.unlight = function(ctx) {
    this.selected = false;
    this.setColor();
    this.render(ctx);
  }
}

// Button function for removing the last drawn line
function undoNodes() {
	circlesArr.pop();
	coords.pop();
  nodeClicked = 0;
	bdisp.clearRect(0, 0, board_display.width, board_display.height);
	updateDraw();
  updateText();
  PREVIEW.preview(coords);
}

// Button function for clearing all the saved lines
function clearNodes() {
	circlesArr = [];
	coords = [];
  nodeClicked = 0;
	bdisp.clearRect(0, 0, board_display.width, board_display.height);
  updateDraw();
  updateText();
  PREVIEW.preview(coords);
}

// From: https://stackoverflow.com/questions/29501447/why-does-css-centering-mess-up-canvas-mouse-coordinates
function getXY(canvas, event) {
    var rect = canvas.getBoundingClientRect();  // absolute position of canvas
    return {
        x: event.clientX - rect.left,
        y: event.clientY - rect.top
    }
}

function inBounds(canvas, event) {
  var rect = canvas.getBoundingClientRect();
  if (event.clientX > rect.left && event.clientX < rect.right && event.clientY < rect.bottom && event.clientY > rect.top) {
    return true;
  }

  else {
    return false;
  }
}

function setError(m) {
	document.getElementById("error-text").innerHTML = m;
}

function bp2radius(bp) {
  return 0.332 * 10 * bp / ( 2 * Math.PI);
}

function radius2bp(r) {
  return Math.round(2 * Math.PI * r / 0.332 / 10);
}

// converts from cartesian coordinates to their canvas position
function cart2canvas(x, y, board) {
  x = Math.round(x * uiscale);
  y = -Math.round((y - ypan) * uiscale);
  return [x, y];
}

// converts a canvas position to cartesian coordinates
function canvas2cart(x, y, board) {
  x = Math.round(x / uiscale);
  y = -Math.round(y / uiscale) + ypan;
  return [x, y];
}

// See https://www.freecodecamp.org/news/render-3d-objects-in-browser-drawing-a-box-with-threejs/
// AJAX call to Flask to convert textarea data into preview diagram
$(document).ready(function() {
    $('#swap').on('click', function() {
        if (INPUTMODE == 'draw') {
            document.getElementById("mesh_txt_input").style.zIndex = "10";
            INPUTMODE = 'text';
            updateText();
            document.getElementById("input_label").innerHTML = "Input (Text)";
        }
        else {
            document.getElementById("mesh_txt_input").style.zIndex = "0";
            INPUTMODE = 'draw';
            document.getElementById("input_label").innerHTML = "Input (Nodes)";
        }
    })

    $('#undo').on('click', function() {
      undoNodes();
    })
    $('#clear').on('click', function() {
      clearNodes();
    })

    $('input[type=radio][name=opt_xovercount]').change(function(){
      clearNodes();
      XOVERMODE = parseInt(document.querySelector('input[name=opt_xovercount]:checked').value);
    })

    $('input[type=number][name=opt_interdist]').change(function(){
      clearNodes();
      INTERHELICAL = parseFloat(document.getElementById('opt_interdist').value);
    })
})

// Checks the total length of the structure against available scaffold length
function checkLength() {
  var totalLength = measureScafLength();
  if (totalLength > SCAFLENLIMIT) {
    printError("Warning: Scaffold length limit reached.")
  }
}

function measureScafLength() {
  var totalLength = 0;
  for (let i=0; i<circlesArr.length; i++) {
    let bp = radius2bp(circlesArr[i].x);
    let cf = XOVERMODE;
    var round = bp % cf;
    bp = round > cf/2 ? bp - (bp % cf) + cf : bp - (bp % cf);
    totalLength += bp;
  }
  return totalLength;
}

function drawAxes() {
  let padding = 30;
  let margin = 40;
  let w = board_axes.width;
  let h = board_axes.height;

  // y axis
  baxes.beginPath();
  baxes.moveTo(padding, margin);
  baxes.lineTo(padding, h);
  baxes.lineWidth = 3;
  baxes.strokeStyle = '#000000';
  baxes.stroke();

  // x axis
  baxes.beginPath();
  baxes.moveTo(0, h - padding);
  baxes.lineTo(w, h - padding);
  baxes.lineWidth = 3;
  baxes.strokeStyle = '#000000';
  baxes.stroke();

  // x axis tick marks
  let maxbp = radius2bp(w / uiscale);
  for (let i=72; i<maxbp; i=i+XOVERMODE * 20) {
    let x = bp2radius(i) * uiscale;
    baxes.beginPath();
    baxes.moveTo(x, h - padding);
    baxes.lineTo(x, h - padding - 5);
    baxes.lineWidth = 3;
    baxes.strokeStyle = '#000000';
    baxes.stroke();
    baxes.font = "15px Helvetica";
    baxes.fillStyle = "black";
    baxes.fillText(i, x - 10, h - padding + 20);
  }

  // y axis tick marks
  let miny = -((h - margin) / 10 / uiscale) + ypan/10;
  let maxy = -(margin / 10 / uiscale) + ypan/10;
  for (let i=miny; i<maxy; i=i+2) {
    let y = -((i * 10 - ypan) * uiscale);
    baxes.beginPath();
    baxes.moveTo(padding, y);
    baxes.lineTo(padding+5, y);
    baxes.lineWidth = 3;
    baxes.strokeStyle = '#000000';
    baxes.stroke();
    baxes.font = "15px Helvetica";
    baxes.fillStyle = "black";
    baxes.fillText(i, padding - 25, y+5);
  }
}

// Gets the canvas elements
// displays existing nodes
var board_display = document.getElementById("layer-show");
var bdisp = board_display.getContext("2d");
// displays cursor node (constantly refreshes)
var board_draw = document.getElementById("layer-draw");
var bdraw = board_draw.getContext("2d");
// displays errors
var board_info = document.getElementById("layer-info");
var binfo = board_info.getContext("2d");
// displays coordinates
var board_xy = document.getElementById("layer-xy");
var bxy = board_xy.getContext("2d");
// displays the axes
var board_axes = document.getElementById("layer-axes");
var baxes = board_axes.getContext("2d");
// displays messages
var board_msg = document.getElementById("layer-msg");
var bmsg = board_msg.getContext("2d");

// Initializes mouse position
var cursorX;
var cursorY;
var firstClick = [0,0];
var intervalLoop = null;
var nodeClicked = 0;
var lastClicked = null;


var defaultBp = 128;
var XOVERMODE = 4;
var INPUTMODE = 'draw';
var INTERHELICAL = 2.6;
var SCAFLENLIMIT = 8064 + 5386 - 500;
const THREECOLOR = "#56B4E9";
const FIVECOLOR = "#009E73";

var circlesArr = [];
var coords = [];
var idx = 0;
var dirBit = 1;
var ypan = 390;

// Initiate a circle
/*
var startCircle = new Circle(...cart2canvas(bp2radius(defaultBp), 0, board_display), 10, idx);
idx = idx + 1;
startCircle.render(bdisp, '#NODECOLOR');
circlesArr.push(startCircle);
var coords = [[0.332*defaultBp/(2*Math.PI), 0, 1]];
*/

/*
// For panning around with middle mouse button
// Controls the offset values for canvas2cart and cart2canvas
var moveOffsetX = 0;
var moveOffsetY = 0;
*/

// display settings
var uiscale = 2;
const circleRadius = 10 * uiscale;

// If continuing from STL nodes, pre-load nodes into page
if (fromstl) {
  coords = existing;
  updateDraw();
}

// placeholder circle for placement of circles in canvas
var newCircle = new Circle(0, 0, circleRadius, 1, dirBit)

// if non-empty
if (coords.length != 0) {
  PREVIEW.preview(coords);
}

drawAxes();
main(); // Begins to monitor mouse movement
