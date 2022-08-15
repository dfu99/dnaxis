// Draw the nodes when the coordinates are returned from Flask
export function preview(coords) {
    if (coords.length == 0) {
      camera.position.set(0, -1 * (midh+90), midh+10);
      camera.lookAt(0, 0, midh);
      renderer.render(scene, camera);
      return;
    }
    scene.remove.apply(scene, scene.children);
    var minh=coords[0][1],
        maxh=coords[0][1],
        minr=coords[0][0],
        maxr=coords[0][0];

    // Find the height center
    for (let i = 0; i < coords.length; i++) {
        var r = coords[i][0],
            h = coords[i][1],
            db = coords[i][2];
        const THREECOLOR = "#56B4E9";
        const FIVECOLOR = "#009E73";
        let c;
        if (db === 0) {c = THREECOLOR;}
        else {c = FIVECOLOR;}
        addHelix(r, h, c);

        minh = Math.min(minh, h);
        maxh = Math.max(maxh, h);
        minr = Math.min(minr, r);
        maxr = Math.max(maxr, r);
    }

    var midh = (minh + maxh) / 2;
    // console.log("Camera set to: "+midh);
    var midr = (minr + maxr) / 2;

    camera.position.set(0, -1 * (midh+90), midh+10);
    camera.lookAt(0, 0, midh);
    renderer.render(scene, camera);
}

function addHelix(r, h, c) {
    // Add circle, torus, and wireframe to represent the helix
    const cgeometry = new THREE.CircleGeometry(1.3, 32);
    const cmaterial = new THREE.MeshBasicMaterial( {color: c} );
    const circle = new THREE.Mesh( cgeometry, cmaterial );
    circle.position.set(r, 0, h);
    circle.rotateX(1.571);
    scene.add (circle);

    const tgeometry = new THREE.TorusGeometry(r, 1.2, 16, 16, Math.PI);
    const tmaterial = new THREE.MeshBasicMaterial( { color: 0xe7e7e7, side: THREE.FrontSide} );
    const torus = new THREE.Mesh( tgeometry, tmaterial );

    const wmaterial = new THREE.LineBasicMaterial( { color: 0x000000 } );
    const wireframe = new THREE.WireframeGeometry( tgeometry );
    const line = new THREE.LineSegments ( wireframe, wmaterial );

    line.material.depthTest = true;
    line.material.opacity = 0.5;
    line.material.transparent = true;

    line.position.set(0, 0, h);
    scene.add( line );

    torus.position.set(0, 0, h);
    scene.add( torus );

    renderer.render(scene, camera);
}

var layer = document.getElementById("layer2");
var myHeight = layer.height;
var myWidth = layer.width;
// Scene
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xffffff);

// Set up lights
const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
// scene.add(ambientLight);

const directionalLight = new THREE.DirectionalLight(0xffffff, 0.2);
directionalLight.position.set(10, 20, 0); // x, y, z
//scene.add(directionalLight);

// Camera
const width = 75;
const height = width * (myHeight / myWidth);
const camera = new THREE.OrthographicCamera(
  width / -2, // left
  width / 2, // right
  height / 2, // top
  height / -2, // bottom
  0.1, // near
  10000 // far
);

// Perspective camera
const aspect = window.innerWidth / window.innerHeight;
const cameraP = new THREE.PerspectiveCamera(
  75, // field of view in degrees
  aspect, // aspect ratio
  1, // near plane
  100 // far plane
);

camera.position.set(0, -100, 20);
camera.lookAt(0, 0, 10);

cameraP.position.set(20, -20, 0);
cameraP.lookAt(0, 0, 0);

// Renderer
const renderer = new THREE.WebGLRenderer({canvas: layer2});
renderer.setSize(myWidth, myHeight);
renderer.render(scene, camera);
