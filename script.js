// Import necessary components from Three.js
// (If using modules, otherwise include via script tags)
// import * as THREE from 'three';
// import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
// import * as dat from 'dat.gui';

// Set up the Three.js scene, camera, and renderer
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xeeeeee); // Light grey background

const sliceDelta = 0.15; // Thickness of the slicing plane

// Create a perspective camera
const camera = new THREE.PerspectiveCamera(
  45, // Field of view
  window.innerWidth / window.innerHeight, // Aspect ratio
  0.1, // Near clipping plane
  1000 // Far clipping plane
);
camera.position.set(0, 0, 20); // Position the camera
camera.lookAt(0, 0, 0); // Make the camera look at the origin

// Create the WebGL renderer
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight); // Set the render size
document.body.appendChild(renderer.domElement); // Add the renderer to the HTML document

// Add OrbitControls for interactive rotation, panning, and zooming
const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableRotate = true;
controls.rotateSpeed = 1.0;
controls.enablePan = true;
controls.enableZoom = true;

// Set controls to rotate with right mouse button
controls.mouseButtons = {
  LEFT: THREE.MOUSE.PAN,
  MIDDLE: THREE.MOUSE.DOLLY,
  RIGHT: THREE.MOUSE.ROTATE,
};

// Parameters for the simulation
const params = {
  extent: 5, // Half-size of the grid
  scale: 3.0, // Scaling factor for the grid
  n: 3, // Principal quantum number
  l: 1, // Azimuthal quantum number (orbital angular momentum quantum number)
  m: 1, // Magnetic quantum number
  amplitudeScale: 2.0, // Amplitude scaling factor for amplitude mapping
  autoRotate: true, // Automatically animate the wave motion
  colorScheme: 'Amplitude', // Color scheme: 'Amplitude', 'Phase'
  omega: 1.0, // Angular frequency for time dependence
  timeScale: 2.0, // Time scaling factor
  t: 0, // Time variable
  sliceAxis: 'None', // Axis along which to slice: 'None', 'X', 'Y', 'Z'
  slicePosition: 0.0, // Position of the slicing plane
  pointSize: 1.0, // Base size of points in the point cloud
  vectorScale: 2, // Scaling factor for vector lengths
  harmonicType: 'Scalar',
  skipPoints: 4, // Sampling rate for vector field visualization
  threshold: 0.01, // Minimum vector length to visualize
  displacementScale: 0.1, // Scaling factor for displacement

};
params.showInfo = function() {
  document.getElementById('info-popup').style.display = 'block';
};

// Global variables
let pointCloud; // The point cloud representing the scalar field (spherical harmonic)
let arrowGroup; // Group to hold vector arrows
let arrowHelpers = []; // Array to store arrow references
let arrowIndices = []; // Array to store indices corresponding to arrows
let highlightedPointIndex = null;
let highlightedPointColor = new THREE.Color(1, 0, 0); // Bright red color

// Create GUI
const gui = new dat.GUI({ autoPlace: false });
document.getElementById('gui-container').appendChild(gui.domElement);

// Declare controllers for l and m
let lController;
let mController;
let scaleController; // Controller for scale parameter

// Function to update m options based on l
function updateMOptions() {
  const mOptions = [];
  for (let i = -params.l; i <= params.l; i++) {
    mOptions.push(i);
  }

  // Remove and recreate mController
  if (mController) {
    gui.remove(mController);
  }
  mController = gui
    .add(params, 'm', mOptions)
    .name('m: Magnetic Quantum Number');
  mController.onChange((value) => {
    params.m = parseInt(value, 10);
    recreateObjects();
  });
}

document.getElementById('info-close').onclick = function() {
  document.getElementById('info-popup').style.display = 'none';
};

window.onclick = function(event) {
  const infoPopup = document.getElementById('info-popup');
  if (event.target === infoPopup) {
    infoPopup.style.display = 'none';
  }
};

// Function to update l and m options based on n
function updateQuantumNumbers() {
  // Update l options based on n
  const lOptions = [];
  for (let i = 0; i <= params.n - 1; i++) {
    lOptions.push(i);
  }

  // Remove and recreate lController
  if (lController) {
    gui.remove(lController);
  }
  lController = gui
    .add(params, 'l', lOptions)
    .name('l: Azimuthal Quantum Number');
  lController.onChange((value) => {
    params.l = parseInt(value, 10);
    // Reset m to 0 if it's out of the new l range
    if (params.m > params.l || params.m < -params.l) {
      params.m = 0;
    }
    updateMOptions();
    recreateObjects();
  });

  // Update m options based on l
  updateMOptions();
}

 gui.add(params, 'showInfo').name('Info');
params.highlightRandomPoint = function() {
  const pointCount = pointCloud.geometry.attributes.position.count;
  highlightedPointIndex = Math.floor(Math.random() * pointCount);
  console .log("highlighted point index "+highlightedPointIndex+
      ", x: "+pointCloud.geometry.attributes.position.array[highlightedPointIndex*3]+", y: "+pointCloud.geometry.attributes.position.array[highlightedPointIndex*3+1]+", z: "+pointCloud.geometry.attributes.position.array[highlightedPointIndex*3+2]);
};

gui.add(params, 'highlightRandomPoint').name('Highlight Random Point');

// Add other GUI controls
gui.add(params, 'amplitudeScale', 0.1, 10).name('Amplitude Scale');
gui.add(params, 'displacementScale', 0.0, 1).name('Displacement Scale');

gui.add(params, 'omega', 0.1, 10).name('Angular Frequency');
gui.add(params, 'timeScale', 1, 10).name('Time Scale');
gui.add(params, 'autoRotate').name('Animate Wave');
gui.add(params, 'colorScheme', ['Amplitude', 'Phase']).name('Color Scheme');
gui
  .add(params, 'sliceAxis', ['None', 'X', 'Y', 'Z'])
  .name('Slice Axis')
  .onChange(recreateObjects);
gui.add(params, 'slicePosition', -10, 10).name('Slice Position');
gui
  .add(params, 'pointSize', 0.1, 5)
  .name('Point Size')
  .onChange(recreateObjects);
gui
  .add(params, 'extent', 4, 15, 1)
  .name('Grid Extent')
  .onChange(recreateObjects);
scaleController = gui
  .add(params, 'scale', 0.2, 5.0)
  .name('Scaling Factor')
  .onChange(() => {
    precomputeWavefunctionData();
    recreateObjects();
  });

// Extend GUI controls for Vector Spherical Harmonics (VSH)
gui
  .add(params, 'vectorScale', 0.1, 5.0)
  .name('Vector Scale')
  .onChange(recreateObjects);
// Add harmonic type control
gui.add(params, 'harmonicType', ['Scalar', 'Vector'])
    .name('Harmonic Type')
    .onChange(() => {
//      adjustGUIControls();
      recreateObjects();
    });

gui
  .add(params, 'skipPoints', 1, 10, 1)
  .name('Skip Points')
  .onChange(recreateObjects);
gui
  .add(params, 'threshold', 0.001, 0.1)
  .name('Threshold')
  .onChange(recreateObjects);

// Add n controller and attach event listener
gui
  .add(params, 'n', 1, 6, 1)
  .name('n: Principal Quantum Number')
  .onChange(() => {
    // Reset l and m to valid values
    if (params.l > params.n - 1) {
      params.l = params.n - 1;
    }
    updateQuantumNumbers();
    recreateObjects();
  });

  
 
// Initialize l and m controllers
updateQuantumNumbers();

// Function to recreate objects when parameters change
function recreateObjects() {
  // Remove existing point cloud
  if (pointCloud) {
    scene.remove(pointCloud);
    pointCloud.geometry.dispose();
    pointCloud.material.dispose();
    pointCloud = null;
  }

  // Remove existing vector field arrows
  if (arrowGroup) {
    scene.remove(arrowGroup);
    arrowGroup.traverse((child) => {
      if (child.geometry) child.geometry.dispose();
      if (child.material) child.material.dispose();
    });
    arrowGroup = null;
  }

  // Reset time
  params.t = 0;

  // Clear caches
  factorialCache = {};
  associatedLegendrePCache = {};
  generalizedLaguerreCache = {};
  sphericalBesselCache = {};
  derivativePthetaCache = {};
  derivativeRhoZnCache = {};

  // Create new objects
  createObjects();
}

// Function to create the point cloud and initialize fields
function createObjects() {
  createPointCloud(); // Create the point cloud representing the scalar field
  precomputeWavefunctionData(); // Precompute per-point data
  if (params.harmonicType === 'Vector' ) {
    computeAndVisualizeVectorField(); // Compute and visualize the vector field
  }
}

// Function to create the point cloud representing the scalar spherical harmonic field
function createPointCloud() {
  const gridSize = params.extent; // Half-size of the grid
  const gridStep = 0.25; // Distance between grid points
  const positions = []; // Array to store positions of points
  const colors = []; // Array to store colors of points
  const sizes = []; // Array to store sizes of points
  const alphas = []; // Array to store alpha (transparency) of points
  const pointData = []; // Array to store per-point data

  // Generate a 3D grid of points
  for (let x = -gridSize; x <= gridSize; x += gridStep) {
    for (let y = -gridSize; y <= gridSize; y += gridStep) {
      for (let z = -gridSize; z <= gridSize; z += gridStep) {
        positions.push(x, y, z);
        colors.push(1, 1, 1); // Initial color (will be updated)
        sizes.push(params.pointSize); // Initial size (will be updated)
        alphas.push(1.0); // Initial alpha (will be updated)

        // Compute and store per-point data
        const r = Math.sqrt(x * x + y * y + z * z);
        const cosTheta = y / (r || 1e-8);
        const theta = Math.acos(cosTheta);
        const phi = Math.atan2(z, x);

        pointData.push({
          x0: x,
          y0: y,
          z0: z,
          r: r,
          theta: theta,
          phi: phi,
          cosTheta: cosTheta,
        });
      }
    }
  }

  const pointCount = positions.length / 3;

  const positionsArray = new Float32Array(positions);
  const colorsArray = new Float32Array(colors);
  const sizesArray = new Float32Array(sizes);
  const alphasArray = new Float32Array(alphas);

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute(
    'position',
    new THREE.BufferAttribute(positionsArray, 3)
  );
  geometry.setAttribute('color', new THREE.BufferAttribute(colorsArray, 3));
  geometry.setAttribute('size', new THREE.BufferAttribute(sizesArray, 1));
  geometry.setAttribute('alpha', new THREE.BufferAttribute(alphasArray, 1));

  // Custom shader material for the point cloud
  const vertexShader = `
    attribute float size;
    attribute vec3 color;
    attribute float alpha;
    varying vec3 vColor;
    varying float vAlpha;
    void main() {
      vColor = color;
      vAlpha = alpha;
      vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
      gl_PointSize = size * (300.0 / -mvPosition.z); // Adjust point size based on depth
      gl_Position = projectionMatrix * mvPosition;
    }
  `;

  const fragmentShader = `
  varying vec3 vColor;
  varying float vAlpha;

  void main() {
    vec2 coord = gl_PointCoord - vec2(0.5);
    float distance = length(coord);

    if (distance > 0.5) {
      discard;
    }

    gl_FragColor = vec4(vColor, vAlpha);
  }
`;

  const material = new THREE.ShaderMaterial({
    vertexShader: vertexShader,
    fragmentShader: fragmentShader,
    transparent: true,
    depthTest: true,
    depthWrite: false,
    blending: THREE.NormalBlending,
  });

  pointCloud = new THREE.Points(geometry, material);


  scene.add(pointCloud);


  // Copy original positions and store per-point data
  pointCloud.geometry.userData = {
    originalPositions: positionsArray.slice(),
    pointData: pointData, // Store per-point data
  };
  // Select a random point to highlight

  highlightedPointIndex = Math.floor(Math.random() * pointCount);
}

// Factorial function with memoization to improve performance
let factorialCache = {};

function factorial(n) {
  if (n < 0) return NaN; // Factorial is not defined for negative numbers
  if (n === 0 || n === 1) return 1; // Base cases
  if (factorialCache[n]) return factorialCache[n]; // Return cached value if available
  let res = 1;
  for (let i = 2; i <= n; i++) {
    res *= i;
  }
  factorialCache[n] = res; // Cache the computed factorial
  return res;
}

// Associated Legendre Polynomial P_l^m(x)
// Computes the value of the associated Legendre polynomial of degree l and order m at x
let associatedLegendrePCache = {};

function associatedLegendreP(l, m, x) {
  // Quantize x to 4 decimal places
  const xKey = x.toFixed(4);
  const key = `${l},${m},${xKey}`;

  if (associatedLegendrePCache[key]) {
    return associatedLegendrePCache[key];
  }

  // Handle the case when m is negative
  m = Math.abs(m);

  let pmm = 1.0; // Initial value P_m^m(x)
  if (m > 0) {
    const somx2 = Math.sqrt((1.0 - x) * (1.0 + x)); // sqrt(1 - x^2)
    let fact = 1.0;
    for (let i = 1; i <= m; i++) {
      pmm *= -fact * somx2; // Recurrence relation
      fact += 2.0;
    }
  }
  if (l === m) {
    associatedLegendrePCache[key] = pmm;
    return pmm; // Return if l equals m
  }

  let pmmp1 = x * (2.0 * m + 1.0) * pmm; // Compute P_{m+1}^m(x)
  if (l === m + 1) {
    associatedLegendrePCache[key] = pmmp1;
    return pmmp1; // Return if l equals m + 1
  }

  let pll = 0.0;
  // Use recurrence relation to compute P_l^m(x) for l > m + 1
  for (let ll = m + 2; ll <= l; ll++) {
    pll =
      ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }

  associatedLegendrePCache[key] = pll; // Cache the result
  return pll;
}

// Derivative of Associated Legendre Polynomial with respect to theta (analytical expression)
let derivativePthetaCache = {};

function derivative_associatedLegendreP_theta(l, m, theta) {
  const cosTheta = Math.cos(theta);
  const sinTheta = Math.sin(theta);

  // Avoid division by zero
  if (Math.abs(sinTheta) < 1e-10) {
    return 0;
  }

  const key = `${l},${m},${theta.toFixed(4)}`;
  if (derivativePthetaCache[key]) {
    return derivativePthetaCache[key];
  }

  // Analytical expression for the derivative
  const P_lm = associatedLegendreP(l, m, cosTheta);
  const P_lm_plus = associatedLegendreP(l + 1, m, cosTheta);

  const dP_dtheta = (l * cosTheta * P_lm - (l + m) * P_lm_plus) / sinTheta;

  derivativePthetaCache[key] = dP_dtheta;
  return dP_dtheta;
}

// Spherical Harmonic Y_l^m(theta, phi)
// Computes the spherical harmonic function Y_l^m at angles theta and phi
function computeSphericalHarmonic(l, m, P_lm, theta, phi, t) {
  const normalization = Math.sqrt(
    ((2 * l + 1) / (4 * Math.PI)) *
      (factorial(l - Math.abs(m)) / factorial(l + Math.abs(m)))
  ); // Normalization constant

  let Y_lm = normalization * P_lm; // Compute the scalar part

  const phase = m * phi - params.omega * t; // Introduce time dependence through phase
  if (m > 0) {
    Y_lm *= Math.cos(phase); // For positive m
  } else if (m < 0) {
    Y_lm *= Math.sin(-phase); // For negative m
  } else {
    Y_lm *= Math.cos(-params.omega * t); // For m = 0
  }
  return Y_lm;
}

// Generalized Laguerre Polynomial L_n^alpha(x)
let generalizedLaguerreCache = {};

function generalizedLaguerreL(n, alpha, x) {
  const xKey = x.toFixed(4);
  const key = `${n},${alpha},${xKey}`;

  if (generalizedLaguerreCache[key]) {
    return generalizedLaguerreCache[key];
  }

  let L = new Array(n + 1);
  L[0] = 1;
  if (n > 0) {
    L[1] = 1 + alpha - x;
  }
  for (let k = 2; k <= n; k++) {
    L[k] =
      ((2 * k - 1 + alpha - x) * L[k - 1] -
        (k - 1 + alpha) * L[k - 2]) /
      k;
  }

  generalizedLaguerreCache[key] = L[n];
  return L[n];
}

// Radial wave function R_{nl}(r)
// Computes the radial part of the hydrogen atom wave function for quantum numbers n and l at radius r
function radialWaveFunction(n, l, r) {
  const Z = 1; // Atomic number for hydrogen
  const a0 = 1; // Bohr radius in atomic units
  const rho = (2 * Z * r) / (n * a0); // Scaled radial coordinate

  const nMinusLMinus1 = n - l - 1;
  const nPlusL = n + l;

  // Check for invalid quantum numbers
  if (nMinusLMinus1 < 0 || nPlusL < 0) {
    return 0;
  }

  // Compute the normalization constant N_{nl}
  const nMinusLMinus1Factorial = factorial(nMinusLMinus1);
  const nPlusLFactorial = factorial(nPlusL);

  const normalization = Math.sqrt(
    ((2 * Z) / (n * a0)) ** 3 *
      (nMinusLMinus1Factorial / (2 * n * nPlusLFactorial))
  );

  // Compute the associated Laguerre polynomial L_{n - l - 1}^{2l + 1}(rho)
  const laguerre = generalizedLaguerreL(nMinusLMinus1, 2 * l + 1, rho);

  // Compute the radial wave function
  const radialPart = normalization * Math.exp(-rho / 2) * rho ** l * laguerre;

  return radialPart;
}

// Spherical Bessel function of the first kind j_n(x)
let sphericalBesselCache = {};

function sphericalBessel_jn(n, x) {
  const xKey = x.toFixed(4);
  const key = `${n},${xKey}`;

  if (sphericalBesselCache[key]) {
    return sphericalBesselCache[key];
  }

  if (n === 0) {
    const result = Math.sin(x) / x; // j_0(x)
    sphericalBesselCache[key] = result;
    return result;
  } else if (n === 1) {
    const result = (Math.sin(x) / (x * x)) - (Math.cos(x) / x); // j_1(x)
    sphericalBesselCache[key] = result;
    return result;
  } else {
    // Use recurrence relation for j_n(x)
    let jn_minus2 = Math.sin(x) / x; // j_0(x)
    let jn_minus1 = (Math.sin(x) / (x * x)) - (Math.cos(x) / x); // j_1(x)
    let jn;
    for (let l = 2; l <= n; l++) {
      jn = ((2 * l - 1) / x) * jn_minus1 - jn_minus2;
      jn_minus2 = jn_minus1;
      jn_minus1 = jn;
    }
    sphericalBesselCache[key] = jn;
    return jn;
  }
}

// Derivative of rho * j_n(rho) with respect to rho
let derivativeRhoZnCache = {};

function derivative_rho_zn(n, rho) {
  const key = `${n},${rho.toFixed(4)}`;
  if (derivativeRhoZnCache[key]) {
    return derivativeRhoZnCache[key];
  }

  const jn = sphericalBessel_jn(n, rho);
  const jn_minus1 = sphericalBessel_jn(n - 1, rho);
  const result = rho * jn_minus1 - n * jn;

  derivativeRhoZnCache[key] = result;
  return result;
}

// Function to precompute wavefunction data
function precomputeWavefunctionData() {
  const pointData = pointCloud.geometry.userData.pointData;
  const n = params.n;
  const l = params.l;
  const m = params.m;
  const scale = params.scale;

  for (let i = 0; i < pointData.length; i++) {
    const data = pointData[i];
    const cosTheta = data.cosTheta;
    const theta = data.theta;
    const r_scaled = data.r * scale;

    // Compute and store associated Legendre polynomial P_l^m(cosTheta)
    data.P_lm = associatedLegendreP(l, m, cosTheta);

    // Compute and store derivative of P_l^m with respect to theta
    data.dP_lm_dtheta = derivative_associatedLegendreP_theta(l, m, theta);

    // Compute and store trigonometric values
    data.sinTheta = Math.sin(theta);
    data.cosTheta = cosTheta;
    data.sinPhi = Math.sin(data.phi);
    data.cosPhi = Math.cos(data.phi);

    // Compute and store radial wave function R_nl(r_scaled)
    data.R_nl = radialWaveFunction(n, l, r_scaled);

    // Compute and store spherical Bessel functions
    const k = 1; // Wave number (can be adjusted)
    const rho = k * data.r;
    data.rho = rho;
    if (rho > 1e-5) {
      data.jn = sphericalBessel_jn(n, rho);
      data.derivative_rhoz_n = derivative_rho_zn(n, rho);
    } else {
      data.jn = 0;
      data.derivative_rhoz_n = 0;
    }
  }
}

// Function to compute the vector at a point (Vector Spherical Harmonics)
// Represents the vector field associated with the spherical harmonic of degree n and order m
function computeVectorAtPoint(n, m, data, t) {
  const r = data.r;
  const theta = data.theta;
  const phi = data.phi;

  const epsilon = 1e-5; // Small value to avoid division by zero
  const maxVectorMagnitude = 1.0; // Maximum allowed vector magnitude

  const rho = data.rho;
  if (rho < epsilon) {
    return new THREE.Vector3(0, 0, 0); // Return zero vector if rho is too small
  }

  // Retrieve precomputed values
  const jn = data.jn;
  const derivative_rhoz_n = data.derivative_rhoz_n;

  // Introduce time dependence in angular functions through phase
  const phase = m * phi - params.omega * t;

  // Compute angular functions with time dependence
  const cos_m_phi = Math.cos(phase);
  const sin_m_phi = Math.sin(phase);

  let sin_theta = data.sinTheta;
  if (Math.abs(sin_theta) < epsilon) {
    sin_theta = epsilon; // Avoid division by zero
  }

  // Compute components in spherical coordinates
  const P_n_m = data.P_lm; // Use precomputed P_n_m
  const dP_n_m_dtheta = data.dP_lm_dtheta; // Use precomputed derivative

  // Radial component N_r
  const N_r = ((n * (n + 1) * jn) / rho) * P_n_m * cos_m_phi;

  // Theta component N_theta
  const N_theta = (derivative_rhoz_n / rho) * dP_n_m_dtheta * cos_m_phi;

  // Phi component N_phi
  const N_phi =
    (-m * (P_n_m / sin_theta) * derivative_rhoz_n) / rho * sin_m_phi;

  // Convert to Cartesian coordinates using precomputed trigonometric values
  const sin_theta_cos_phi = sin_theta * data.cosPhi;
  const sin_theta_sin_phi = sin_theta * data.sinPhi;
  const cos_theta = data.cosTheta;

  const e_r = new THREE.Vector3(
    sin_theta_cos_phi,
    sin_theta_sin_phi,
    cos_theta
  ); // Radial unit vector
  const e_theta = new THREE.Vector3(
    cos_theta * data.cosPhi,
    cos_theta * data.sinPhi,
    -sin_theta
  ); // Theta unit vector
  const e_phi = new THREE.Vector3(-data.sinPhi, data.cosPhi, 0); // Phi unit vector

  const N_vector = new THREE.Vector3(0, 0, 0);
  N_vector.addScaledVector(e_r, N_r);
  N_vector.addScaledVector(e_theta, N_theta);
  N_vector.addScaledVector(e_phi, N_phi);

  // Cap the vector magnitude to prevent excessively large vectors
  N_vector.clampLength(0, maxVectorMagnitude);

  return N_vector;
}
function updatePointCloud() {

    const positionAttribute = pointCloud.geometry.attributes.position;
    const colorAttribute = pointCloud.geometry.attributes.color;
    const sizeAttribute = pointCloud.geometry.attributes.size;
    const alphaAttribute = pointCloud.geometry.attributes.alpha;
    const pointCount = positionAttribute.count;

    const n = params.n;
    const l = params.l;
    const m = params.m;
    const pointData = pointCloud.geometry.userData.pointData;

    for (let i = 0; i < pointCount; i++) {
      const data = pointData[i];
      const x0 = data.x0;
      const y0 = data.y0;
      const z0 = data.z0;

      const P_lm = data.P_lm;
      const R_nl = data.R_nl;

      // Compute the phase
      const phase = m * data.phi - params.omega * params.t;

      let amplitude;
      let displacementVector = new THREE.Vector3(0, 0, 0);
      let color = new THREE.Color();
      let size;
      if (params.harmonicType === 'Scalar') {
        // Scalar Harmonic: Compute Psi_nlm
        let Y_lm = computeSphericalHarmonic(
            l,
            m,
            P_lm,
            data.theta,
            data.phi,
            params.t
        );

        const Psi_nlm = R_nl * Y_lm;
        amplitude = Math.abs(Psi_nlm) * 50; // Scale amplitude for visualization

        // Radial displacement
        const displacementMagnitude = amplitude * params.displacementScale;
        const r = data.r;
        let unitRadialVector = new THREE.Vector3(x0, y0, z0);
        if (r > 1e-5) {
          unitRadialVector.divideScalar(r);
        } else {
          unitRadialVector.set(0, 0, 0);
        }
        displacementVector = unitRadialVector.multiplyScalar(displacementMagnitude);

        // Color based on amplitude or phase
        if (params.colorScheme === 'Amplitude') {
          const hue = 0.7 - amplitude * 0.7; // Hue from blue (0.7) to red (0.0)
          const lightness = 0.4 + amplitude * 0.6; // Lightness from 0.4 to 1.0
          color.setHSL(hue, 1, lightness);
        } else if (params.colorScheme === 'Phase') {
          const phaseValue =
              (Math.atan2(
                      Math.sin(m * data.phi - params.omega * params.t),
                      Math.cos(m * data.phi - params.omega * params.t)
                  ) +
                  Math.PI) /
              (2 * Math.PI); // Normalize to [0,1]
          color.setHSL(phaseValue, 1, 0.5); // Hue based on phase
        }
        size = params.pointSize * (1 + amplitude * params.amplitudeScale);
      } else if (params.harmonicType === 'Vector') {
        // Vector Harmonic: Compute vector displacement
        const vectorDisplacement = computeVectorAtPoint(n, m, data, params.t)
            .multiplyScalar(params.displacementScale);
        displacementVector.copy(vectorDisplacement);

        // Compute amplitude based on vector magnitude
        amplitude = vectorDisplacement.length() *5; // Scale for visualization

        // Color based on vector properties (e.g., direction)
        const phi = data.phi;
        const normalizedPhi = (phi + Math.PI) / (2 * Math.PI); // Normalize phi to [0,1]
        const hue = normalizedPhi; // Hue varies from 0 to 1 based on phi
        color.setHSL(hue, 1, 0.5);
        size = 0.5*params.pointSize * (1 + amplitude * params.amplitudeScale);
      }

      // Limit the size of the displacement vector
      displacementVector.clampLength(0, 1.0);

      // Update position
      const newX = x0 + displacementVector.x;
      const newY = y0 + displacementVector.y;
      const newZ = z0 + displacementVector.z;
      positionAttribute.setXYZ(i, newX, newY, newZ);

      // Adjust point size based on amplitude

      sizeAttribute.setX(i, size);

      // Apply slicing to determine visibility
      let visible = true;
      if (params.sliceAxis !== 'None') {
        const sliceValue = params.slicePosition;
        switch (params.sliceAxis) {
          case 'X':
            visible = Math.abs(x0 - sliceValue) < sliceDelta;
            break;
          case 'Y':
            visible = Math.abs(y0 - sliceValue) < sliceDelta;
            break;
          case 'Z':
            visible = Math.abs(z0 - sliceValue) < sliceDelta;
            break;
        }
      }

      // Set transparency based on visibility and amplitude
      let alpha = visible ? Math.min(1.0, amplitude - 0.1) : 0;
      alphaAttribute.setX(i, alpha);

      // Highlighted point handling remains the same
      if (i === highlightedPointIndex) {
        // Highlighted point: set color to bright red and increase size
        colorAttribute.setXYZ(i, highlightedPointColor.r, highlightedPointColor.g, highlightedPointColor.b);
        size *= 1.5; // Increase size by 50% for visibility
        size = Math.max(size, 0.01);
        sizeAttribute.setX(i, size);
        alphaAttribute.setX(i, 0.6); // Set transparency to 0.6
      } else {
        // Set color for regular points
        colorAttribute.setXYZ(i, color.r, color.g, color.b);
      }
    }

    // Mark attributes as needing updates
    positionAttribute.needsUpdate = true;
    colorAttribute.needsUpdate = true;
    sizeAttribute.needsUpdate = true;
    alphaAttribute.needsUpdate = true;
  }



// Function to compute and visualize the vector field
function computeAndVisualizeVectorField() {
  const positionsArray = pointCloud.geometry.attributes.position.array;
  const numPoints = positionsArray.length / 3;

  const pointData = pointCloud.geometry.userData.pointData;
  const n = params.n;
  const m = params.m;

  arrowGroup = new THREE.Group();
  arrowHelpers = [];
  arrowIndices = [];

  let skip = params.skipPoints;
  // Show all vectors if slicing is applied
  if (params.sliceAxis !== 'None') {
    skip = 1;
  }

  const threshold = params.threshold;
  const scaleFactor = params.vectorScale;
  const maxVectorMagnitude = 1.0; // Should match the one in computeVectorAtPoint

  const spacing = params.skipPoints * 0.25; // Use 'skipPoints' to control the grid spacing

  for (let i = 0; i < numPoints; i++) {
    const data = pointData[i];
    const x = data.x0;
    const y = data.y0;
    const z = data.z0;

    // Check if the point is evenly spaced based on grid coordinates
    if (
      Math.abs(x % spacing) > 0.01 ||
      Math.abs(y % spacing) > 0.01 ||
      Math.abs(z % spacing) > 0.01
    ) {
      continue; // Skip this point if it does not align with the desired grid spacing
    }

    // Compute the vector at this point
    const vector = computeVectorAtPoint(n, m, data, params.t);

    const magnitude = vector.length();
    const length = magnitude * scaleFactor;

    // Skip vectors below the threshold
    if (length > threshold) {
      const dir = vector.clone().normalize();
      const origin = new THREE.Vector3(x, y, z);

      // Compute color based on vector direction
      const phi = data.phi;
      const normalizedPhi = (phi + Math.PI) / (2 * Math.PI); // Normalize phi to [0,1]
      const hue = normalizedPhi; // Hue varies from 0 to 1 based on phi
      const color = new THREE.Color();
      color.setHSL(hue, 1, 0.5);

      // Create arrow helper with the computed direction, origin, length, and color
      const arrowHelper = new THREE.ArrowHelper(
        dir,
        origin,
        length,
        color.getHex()
      );

      // Apply slicing
      let visible = true;

      if (params.sliceAxis !== 'None') {
        const sliceValue = params.slicePosition;
        switch (params.sliceAxis) {
          case 'X':
            visible = Math.abs(x - sliceValue) < sliceDelta;
            break;
          case 'Y':
            visible = Math.abs(y - sliceValue) < sliceDelta;
            break;
          case 'Z':
            visible = Math.abs(z - sliceValue) < sliceDelta;
            break;
        }
      }

      arrowHelper.visible = visible;

      arrowGroup.add(arrowHelper);
      arrowHelpers.push(arrowHelper);
      arrowIndices.push(i); // Store the index corresponding to this arrow
    }
  }

  scene.add(arrowGroup);
}

// Function to update the arrow helpers in each animation frame
function updateArrowHelpers() {
  const scaleFactor = params.vectorScale;
  const maxVectorMagnitude = 1.0; // Should match the one in computeVectorAtPoint
  const pointData = pointCloud.geometry.userData.pointData;
  const n = params.n;
  const m = params.m;

  for (let i = 0; i < arrowHelpers.length; i++) {
    const arrowHelper = arrowHelpers[i];
    const index = arrowIndices[i];

    const data = pointData[index];
    const x = data.x0;
    const y = data.y0;
    const z = data.z0;

    // Compute the vector at this point with current time
    const vector = computeVectorAtPoint(n, m, data, params.t);
    const magnitude = vector.length();
    const length = magnitude * scaleFactor;

    // Update arrow direction and length
    if (length > params.threshold) {
      const dir = vector.clone().normalize();

      arrowHelper.setDirection(dir);
      arrowHelper.setLength(length);

      // Compute color based on vector direction (phi angle)
      const phi = data.phi;
      const normalizedPhi = (phi + Math.PI) / (2 * Math.PI); // Normalize phi to [0,1]
      const hue = normalizedPhi; // Hue varies from 0 to 1 based on phi
      const color = new THREE.Color();
      color.setHSL(hue, 1, 0.5);

      // Update arrow color
      if (arrowHelper.setColor) {
        arrowHelper.setColor(color);
      } else {
        // Access the materials directly if setColor is not available
        arrowHelper.line.material.color.copy(color);
        arrowHelper.cone.material.color.copy(color);
      }

      // Update visibility based on slicing
      let visible = true;
      if (params.sliceAxis !== 'None') {
        const sliceValue = params.slicePosition;
        switch (params.sliceAxis) {
          case 'X':
            visible = Math.abs(x - sliceValue) < sliceDelta;
            break;
          case 'Y':
            visible = Math.abs(y - sliceValue) < sliceDelta;
            break;
          case 'Z':
            visible = Math.abs(z - sliceValue) < sliceDelta;
            break;
        }
      }
      arrowHelper.visible = visible;
    } else {
      arrowHelper.visible = false;
    }
  }
}

// Animation loop
function animate() {
  requestAnimationFrame(animate);

  // Update time variable
  if (params.autoRotate) {
    const dt = 0.01 * params.timeScale;
    params.t += dt;
  }


    updatePointCloud(); // Update the scalar field and positions


  // Update arrow helpers
  if (params.harmonicType === 'Vector'  && arrowHelpers.length > 0) {
    updateArrowHelpers(); // Update the vector field
  }

  // Render the scene
  controls.update();
  renderer.render(scene, camera);
}


// Start the simulation
createObjects();
animate();

// Handle window resize events
window.addEventListener('resize', onWindowResize, false);

function onWindowResize() {
  // Update camera aspect ratio and renderer size
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
}
