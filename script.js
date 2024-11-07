// Declare 'params' as a global variable
let params = {};

// Begin waveFunctions.js content

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

// Derivative of Associated Legendre Polynomial with respect to theta
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
    const P_n_m = data.P_lm; // Use precomputed P_n_m
    const dP_n_m_dtheta = data.dP_lm_dtheta; // Use precomputed derivative

    // Introduce time dependence through phase
    const phase = m * phi - params.omega * t;

    // Compute angular functions with time dependence
    const cos_m_phi = Math.cos(phase);
    const sin_m_phi = Math.sin(phase);

    let sin_theta = data.sinTheta;
    if (Math.abs(sin_theta) < epsilon) {
        sin_theta = epsilon; // Avoid division by zero
    }

    let N_r = 0, N_theta = 0, N_phi = 0;

    if (params.harmonicType === 'Magnetic Vector Harmonics') {
        // Magnetic Vector Spherical Harmonics

        // Radial component is zero
        N_r = 0;

        // Theta component
        N_theta = (data.jn / sin_theta) * m * P_n_m * sin_m_phi;

        // Phi component
        N_phi = -data.jn * dP_n_m_dtheta * cos_m_phi;

    } else if (params.harmonicType === 'Electric Vector Harmonics') {
        // Electric Vector Spherical Harmonics

        // Compute derivative of spherical Bessel function
        const jn = data.jn;
        const djn_drho = data.djn_drho;

        // Radial component
        N_r = n * (n + 1) * (jn / rho) * P_n_m * cos_m_phi;

        // Theta component
        N_theta = djn_drho * dP_n_m_dtheta * cos_m_phi;

        // Phi component
        N_phi = djn_drho * (m * P_n_m / sin_theta) * sin_m_phi;

    } else {
        // For other harmonic types, return zero vector
        return new THREE.Vector3(0, 0, 0);
    }

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


function clearCaches() {
    // Clear caches
    factorialCache = {};
    associatedLegendrePCache = {};
    generalizedLaguerreCache = {};
    sphericalBesselCache = {};
    derivativePthetaCache = {};
    derivativeRhoZnCache = {};
}

// End of waveFunctions.js content
// ========================================================================================================
// Set up the Three.js scene, camera, and renderer
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xffffff); // Light grey background

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

// Global variables
let pointCloud; // The point cloud representing the scalar field (spherical harmonic)
let arrowGroup; // Group to hold vector arrows
let arrowHelpers = []; // Array to store arrow references
let arrowIndices = []; // Array to store indices corresponding to arrows

let highlightedPointColor = new THREE.Color(9.9,0.2, 0); // Bright red color
let highlightedPointIndices = []; // Array of indices
let highlightCubes = []; // Array of cube meshes
let markerTrails = []; // Array to store trail lines for each highlighted point
function createHighlightCubes() {
    // Remove existing highlight cubes and trails
    if (highlightCubes.length > 0) {
        for (let cube of highlightCubes) {
            scene.remove(cube);
            cube.geometry.dispose();
            cube.material.dispose();
        }
        highlightCubes = [];
    }

    if (markerTrails.length > 0) {
        for (let trail of markerTrails) {
            scene.remove(trail);
            trail.geometry.dispose();
            trail.material.dispose();
        }
        markerTrails = [];
    }

    // Create new cubes and trails
    for (let i = 0; i < highlightedPointIndices.length; i++) {
        // Create cube
        const geometry = new THREE.BoxGeometry(0.2, 0.2, 0.2); // Adjust size as needed
        const material = new THREE.MeshBasicMaterial({
            color: highlightedPointColor,
            opacity: 0.7,
            transparent: true,
        });
        const cube = new THREE.Mesh(geometry, material);
        scene.add(cube);
        highlightCubes.push(cube);

        // Create trail for the cube
        const trailMaterial = new THREE.LineBasicMaterial({ color: highlightedPointColor });
        const trailGeometry = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3()]);
        const trail = new THREE.Line(trailGeometry, trailMaterial);
        trail.userData.positions = []; // To store trail positions
        scene.add(trail);
        markerTrails.push(trail);
    }
}

function highlightRandomPoints() {
    if (!pointCloud) return;

    const pointCount = pointCloud.geometry.attributes.position.count;
    highlightedPointIndices = [];

    // Generate unique random indices
    const indicesSet = new Set();
    while (indicesSet.size < params.numHighlightedPoints && indicesSet.size < pointCount) {
        const index = Math.floor(Math.random() * pointCount);
        indicesSet.add(index);
    }

    highlightedPointIndices = Array.from(indicesSet);

    // Create the highlight cubes
    createHighlightCubes();
}
// Function to recreate objects when parameters change
function recreateObjects() {
    // Remove existing point cloud
    if (pointCloud) {
        scene.remove(pointCloud);
        pointCloud.geometry.dispose();
        pointCloud.material.dispose();
        pointCloud = null;
    }
    // Remove existing highlight cubes and trails
    if (highlightCubes.length > 0) {
        for (let cube of highlightCubes) {
            scene.remove(cube);
            cube.geometry.dispose();
            cube.material.dispose();
        }
        highlightCubes = [];
    }

    if (markerTrails.length > 0) {
        for (let trail of markerTrails) {
            scene.remove(trail);
            trail.geometry.dispose();
            trail.material.dispose();
        }
        markerTrails = [];
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
    clearCaches();

    // Create new objects
    createObjects();
}

// Function to create the point cloud and initialize fields
function createObjects() {
    createPointCloud(); // Create the point cloud representing the scalar field
    highlightRandomPoints();
    precomputeWavefunctionData(); // Precompute per-point data
    if (
        params.harmonicType === 'Magnetic Vector Harmonics' ||
        params.harmonicType === 'Electric Vector Harmonics'
    ) {
        computeAndVisualizeVectorField(); // Compute and visualize the vector field
    }
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
          const sliceDelta = params.sliceWidth;
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
function computeKernelQuaternion(originalPosition, radialUnitVector) {
    let angle = computeChargeStrength(originalPosition);

    // Compute an axis orthogonal to the radial unit vector
    let rotationAxis = new THREE.Vector3();
    const xAxis = new THREE.Vector3(1, 0, 0);
    rotationAxis.crossVectors(radialUnitVector, xAxis);

    // If radial vector is parallel to xAxis, use yAxis instead
    if (rotationAxis.lengthSq() < 1e-8) {
        const yAxis = new THREE.Vector3(0, 1, 0);
        rotationAxis.crossVectors(radialUnitVector, yAxis);
    }
    rotationAxis.normalize();

    // Create the kernel rotation quaternion k(r)
    const kernelQuaternion = new THREE.Quaternion();
    kernelQuaternion.setFromAxisAngle(rotationAxis, angle);

    return kernelQuaternion;
}

function exponentialSmoothing(r, radius, width) {
  return Math.exp(-((r - radius) ** 2) / (2 * width ** 2));
}

function computeTimeQuaternion(radialUnitVector, phase) {
    // Create the time-based rotation quaternion q(t)
    const timeQuaternion = new THREE.Quaternion();
    timeQuaternion.setFromAxisAngle(radialUnitVector, phase);

    return timeQuaternion;
}
function transformCharge(originalPosition) {
    // Compute radial distance from the origin
    const r = originalPosition.length();

    // Define wave parameters
    const k =  4*Math.PI / params.waveNumber;  // Wave number, adjust to control wavelength
    const omega = params.timeScale ;           // Angular frequency, adjust speed as needed

    // Determine pulsation direction and calculate scaling factor
    let scalingFactor=1;
    // function if distance, it should diminish the further away it is
    let stretch = 0.5 / (1 + r);

    if (params.chargeMode === 'Plus Charge') {
        // Outward pulsation
        scalingFactor = 1.0 + stretch * Math.sin(k * r - omega * params.phase);
    } else {
        // Inward pulsation
        scalingFactor = 1.0 + stretch * Math.sin(k * r + omega * params.phase);
    }
  //  scalingFactor= scalingFactor*params.displacementScale;
    // Ensure scaling factor remains positive
    scalingFactor = Math.max(0.01, scalingFactor);
    // also make sure it doesn't go over 2
    if (scalingFactor >1.5) {
        scalingFactor = 1.5;
    }

    // Apply the scaling to the original position to create a pulsating effect
    const newPosition = originalPosition.clone().multiplyScalar(scalingFactor);

    return {newPosition, scalingFactor};
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
          originalPosition: new THREE.Vector3(x, y, z),
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


}
// Modify precomputeWavefunctionData to compute the derivative of the spherical Bessel function
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
            data.jn_plus1 = sphericalBessel_jn(n + 1, rho); // For derivative
            data.djn_drho = (n / rho) * data.jn - data.jn_plus1; // Derivative of jn with respect to rho
            data.derivative_rhoz_n = derivative_rho_zn(n, rho);
        } else {
            data.jn = 0;
            data.djn_drho = 0;
            data.derivative_rhoz_n = 0;
        }
    }
}

// Function to apply the spin 1/2 transformation to a position using quaternions
function spinHalfTransform(originalPosition, timeQuaternion, inverseTimeQuaternion) {
  const radius = originalPosition.length(); // Distance from the origin
  let kernelRotationAngle = 0; // Rotation angle for the kernel
  const angle = (params.maxKernelAngle / 180.0) * Math.PI; // Convert max angle to radians
  const power = params.decayPower; // Decay power for the kernel function

  // Original kernel function: Decaying rotation angle with distance
  const adjustedRadius = Math.max(radius, 0);
  // Fraction decreases with distance using a power function
  const fraction = 1 / Math.pow(adjustedRadius + 1, power);
  kernelRotationAngle = angle * fraction;

  // Kernel rotation quaternion k(r) rotating around X-axis
  const kernelQuaternion = new THREE.Quaternion();
  kernelQuaternion.setFromAxisAngle(
      new THREE.Vector3(1, 0, 0), // Rotation axis (X-axis)
      kernelRotationAngle          // Rotation angle based on kernel
  );

  // Combined rotation: q(t) * k(r) * q⁻¹(t)
  const combinedQuaternion = new THREE.Quaternion();
  combinedQuaternion.multiplyQuaternions(timeQuaternion, kernelQuaternion);
  combinedQuaternion.multiply(inverseTimeQuaternion);

  // Apply the combined rotation to the original position
  const newPosition = originalPosition.clone().applyQuaternion(combinedQuaternion);

  return { newPosition, combinedQuaternion };
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

  // Time-based rotation quaternion q(t) around Y-axis for spin 1/2
  const timeQuaternion = new THREE.Quaternion();
  timeQuaternion.setFromAxisAngle(
      new THREE.Vector3(0, 1, 0), // Rotation axis (Y-axis)
      params.phase                // Rotation angle (phase)
  );

  // Inverse of q(t), which is q⁻¹(t) or conjugate in quaternion terms
  const inverseTimeQuaternion = timeQuaternion.clone().conjugate();

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

    // Original position vector
    let originalPosition = data.originalPosition.clone();
	const isCharge = params.harmonicType === 'Charge Only' ;
    // Apply spin 1/2 transformation if enabled
    if (params.spinorMode !== 'Off' || params.harmonicType === 'Spinor Only') {
      const transformResult = spinHalfTransform(originalPosition, timeQuaternion, inverseTimeQuaternion);
      originalPosition = transformResult.newPosition;
    }
	 if (isCharge) {
	  const transformResult = transformCharge(originalPosition);
	  originalPosition = transformResult.newPosition;
      let scalingFactor = transformResult.scalingFactor;
     const minScale = 0.75; // Adjust this range as necessary
     const maxScale = 1.25;
     const normalizedScale = (scalingFactor - minScale) / (maxScale - minScale);

     // Clamp normalizedScale to [0, 1] to avoid overflows
     const clampedScale = Math.max(0, Math.min(1, normalizedScale));
     const hue = clampedScale * 0.7; // Map normalized scale to hue (0 is red, 0.7 is greenish)

         color.setHSL(hue, 1.0, 0.5); // Adjust saturation and lightness as needed
	}
	 // Apply spin 1/2 transformation if enabled
        if (params.applySpinHalf) {
            const transformResult = spinHalfTransform(originalPosition, timeQuaternion, inverseTimeQuaternion);
            originalPosition = transformResult.newPosition;
        }

        // Apply charge transformation if enabled
        if (params.applyChargeTransformation) {
            const transformResult = transformCharge(originalPosition);
            originalPosition = transformResult.newPosition;
            let scalingFactor = transformResult.scalingFactor;

            // Compute color based on scalingFactor
            const minScale = 0.75; // Adjust this range as necessary
            const maxScale = 1.25;
            const normalizedScale = (scalingFactor - minScale) / (maxScale - minScale);

            // Clamp normalizedScale to [0, 1] to avoid overflows
            const clampedScale = Math.max(0, Math.min(1, normalizedScale));
            const hue = clampedScale * 0.7; // Map normalized scale to hue (0 is red, 0.7 is greenish)

            color.setHSL(hue, 1.0, 0.5); // Adjust saturation and lightness as needed
        }
    if (params.harmonicType === 'Scalar') {
      // Scalar Harmonic: Compute Psi_nlm
      let Y_lm = computeSphericalHarmonic(
          l,  m,     P_lm,      data.theta,      data.phi, params.t
      );

      const Psi_nlm = R_nl * Y_lm;
      amplitude = Math.abs(Psi_nlm) * 50; // Scale amplitude for visualization

      // Radial displacement
      const displacementMagnitude = amplitude * params.displacementScale*2;
      const r = originalPosition.length();
      let unitRadialVector = originalPosition.clone().normalize();
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
    } else if (
            params.harmonicType === 'Magnetic Vector Harmonics' ||
            params.harmonicType === 'Electric Vector Harmonics'
        ) {
            // Vector Harmonics: Compute vector displacement
            const vectorDisplacement = computeVectorAtPoint(n, m, data, params.t).multiplyScalar(
                params.displacementScale * 2
            );
            displacementVector.copy(vectorDisplacement);

            // Compute amplitude based on vector magnitude
            amplitude = vectorDisplacement.length() * 5; // Scale for visualization

            // Color based on vector properties (e.g., direction)
            const phi = data.phi;
            const normalizedPhi = (phi + Math.PI) / (2 * Math.PI); // Normalize phi to [0,1]
            const hue = normalizedPhi; // Hue varies from 0 to 1 based on phi
            color.setHSL(hue, 1, 0.5);
            size = 0.5 * params.pointSize * (1 + amplitude * params.amplitudeScale);
       
    } else if (params.harmonicType === 'Spinor Only') {
      // Spinor Only: Color based on position after spin 1/2 transformation
      amplitude = originalPosition.length() * 0.3;
      let normAmp =(amplitude + params.extent) / (2 * params.extent);
      color.setHSL(normAmp, 1, 0.5);
      size = params.pointSize*0.5;
      amplitude =normAmp;
    
    } else if (params.harmonicType === 'Charge Only') {
      // Charge Only: Color based on position after transformation
      amplitude = originalPosition.length() * 0.5;
      let normAmp =(amplitude + params.extent) / (2 * params.extent);
      color.setHSL(normAmp, 1, 0.5);
      size = 0.5 * params.pointSize * (1 + amplitude * params.amplitudeScale);
     // amplitude =normAmp;
    }

    // Limit the size of the displacement vector
    displacementVector.clampLength(0, 1.0);

    // Update position
    const newX = originalPosition.x + displacementVector.x;
    const newY = originalPosition.y + displacementVector.y;
    const newZ = originalPosition.z + displacementVector.z;
    positionAttribute.setXYZ(i, newX, newY, newZ);

    // Adjust point size based on amplitude
    sizeAttribute.setX(i, size);

    // Apply slicing to determine visibility
    let visible = true;
    if (params.sliceAxis !== 'None') {
      const sliceValue = params.slicePosition;
      const sliceDelta = params.sliceWidth;
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
     let alpha;
    if (params.harmonicType === 'Spinor Only' || params.harmonicType === 'Charge Only' ) {
		// make inside less transparent
		if (visible) {
            alpha = 0.3;
        }
        else {
            alpha = 0;
        }
      // In Spinor Only mode, set alpha to a lower value for transparency

    } else {
      // Set transparency based on visibility and amplitude for other modes
      alpha = visible ? Math.min(1.0, amplitude - 0.1) : 0;
    }
    alphaAttribute.setX(i, alpha);

  // Set color for regular points
  colorAttribute.setXYZ(i, color.r, color.g, color.b);


  }
    // After updating positions, update the cubes' positions and trails
    if (highlightedPointIndices.length > 0) {
        const positions = positionAttribute.array;
        const trailLength = params.trailLength;

        for (let i = 0; i < highlightedPointIndices.length; i++) {
            const index = highlightedPointIndices[i] * 3;
            const x = positions[index];
            const y = positions[index + 1];
            const z = positions[index + 2];
            const newPosition = new THREE.Vector3(x, y, z);

            // Update cube position
            highlightCubes[i].position.copy(newPosition);

            // Update trail
            const trail = markerTrails[i];
            const positionsArray = trail.userData.positions;

            positionsArray.push(newPosition.clone());
            if (positionsArray.length > trailLength) {
                positionsArray.shift(); // Remove oldest position to maintain trail length
            }

            // Update trail geometry
            trail.geometry.setFromPoints(positionsArray);
        }
    }
  // Mark attributes as needing updates
  positionAttribute.needsUpdate = true;
  colorAttribute.needsUpdate = true;
  sizeAttribute.needsUpdate = true;
  alphaAttribute.needsUpdate = true;
}

// Function to update the arrow helpers in each animation frame
function updateArrowHelpers(timeQuaternion, inverseTimeQuaternion) {
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

      // Apply spin 1/2 transformation if enabled
      let origin = new THREE.Vector3(x, y, z);
      if (
          params.spinHalf &&
          timeQuaternion !== undefined &&
          inverseTimeQuaternion !== undefined
      ) {
        origin = spinHalfTransform(origin, timeQuaternion, inverseTimeQuaternion);
        dir.applyQuaternion(
            timeQuaternion.clone().multiply(inverseTimeQuaternion)
        );
      }

      arrowHelper.position.copy(origin);
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
        const sliceDelta = params.sliceWidth;
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

    // Update phase angle for spin 1/2
     params.phase = params.t;
  }

  updatePointCloud(); // Update the scalar field and positions

  // Update arrow helpers
  if (
        (params.harmonicType === 'Magnetic Vector Harmonics' ||
            params.harmonicType === 'Electric Vector Harmonics') &&
        arrowHelpers.length > 0
    ) {
        updateArrowHelpers(); // Update the vector field
    }

  // Render the scene
  controls.update();
  renderer.render(scene, camera);
}

// Handle window resize events
window.addEventListener('resize', onWindowResize, false);

function onWindowResize() {
  // Update camera aspect ratio and renderer size
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
}

// End of script.js content

// Begin guicontrols.js content

let gui;
let quantumFolder, visualizationFolder, spinFolder, chargeFolder, sphericalFolder, vectorFolder;
let mController, lController, scaleController;

// Function to update m options based on l
function updateMOptions() {
    if (!params) return; // Ensure params is initialized

    const mOptions = [];
    for (let i = -params.l; i <= params.l; i++) {
        mOptions.push(i);
    }

    if (mController) {
        quantumFolder.remove(mController);
    }
    mController = quantumFolder
        .add(params, 'm', mOptions)
        .name('m: Magnetic Quantum Number');
    mController.onChange((value) => {
        params.m = parseInt(value, 10);
        recreateObjects();
    });
}

// Function to highlight a random point
function highlightRandomPoint() {
    if (!pointCloud) return; // Ensure pointCloud is available

    const pointCount = pointCloud.geometry.attributes.position.count;
    highlightedPointIndex = Math.floor(Math.random() * pointCount);

}

// Function to update l and m options based on n
function updateQuantumNumbers() {
    if (!params || !quantumFolder) return; // Ensure params and quantumFolder are initialized

    const lOptions = [];
    for (let i = 0; i <= params.n - 1; i++) {
        lOptions.push(i);
    }

    if (lController) {
        quantumFolder.remove(lController);
    }
    lController = quantumFolder
        .add(params, 'l', lOptions)
        .name('l: Azimuthal Quantum Number');
    lController.onChange((value) => {
        params.l = parseInt(value, 10);
        if (params.m > params.l || params.m < -params.l) {
            params.m = 0;
        }
        updateMOptions();
        recreateObjects();
    });

    updateMOptions();

    // if n is 0, make the amplitude scale smaller
    if (params.n === 0) {
        params.amplitudeScale = 0.005;
    }
}

// Adjust the GUI to include the new harmonic type
function initializeGUI() {
    params = {
		  applySpinHalf: false,
    applyChargeTransformation: false,
        extent: 5,
        scale: 3.0,
        n: 3,
        l: 1,
        m: 1,
        amplitudeScale: 2.0,
        autoRotate: true,
        colorScheme: 'Amplitude',
        omega: 1.0,
        timeScale: 2.0,
        t: 0,
        sliceAxis: 'None',
        slicePosition: 0.0,
        sliceWidth: 0.12,
        pointSize: 1.0,
        vectorScale: 2,
        harmonicType: 'Scalar',
        skipPoints: 4,
        threshold: 0.01,
        displacementScale: 0.1,
        numHighlightedPoints: 10,
        trailLength: 50,
        showInfo: function () {
            document.getElementById('info-popup').style.display = 'block';
        },
        highlightRandomPoint,
        enableSpinHalf: false,
        spinorMode: 'Off',
        chargeMode: 'Plus Charge',
        waveNumber: 4,
        maxKernelAngle: 180,
        decayPower: 2,
        phase: 0.0,
        sphereRadius: 2.0,
    };

    gui = new dat.GUI({ autoPlace: false });

    // Group GUI controls into folders to reduce clutter

    visualizationFolder = gui.addFolder('Visualization');

    // Visualization controls
    visualizationFolder.add(params, 'harmonicType', [
        'Scalar',
        'Magnetic Vector Harmonics',
        'Electric Vector Harmonics',
        'Spinor Only',
        'Charge Only',
    ])
        .name('Wave Type')
        .onChange(() => {
            adjustGUIControls();
            recreateObjects();
        });
		
		visualizationFolder.add(params, 'applySpinHalf')
    .name('Apply Spin 1/2 Transformation')
    .onChange(() => {
        adjustGUIControls();
        recreateObjects();
    });

visualizationFolder.add(params, 'applyChargeTransformation')
    .name('Apply Charge Transformation')
    .onChange(() => {
        adjustGUIControls();
        recreateObjects();
    });
	
    quantumFolder = gui.addFolder('Quantum Numbers');
    spinFolder = gui.addFolder('Spin 1/2 Transformation');
    spinFolder.closed = false;

    // Info button
    gui.add(params, 'showInfo').name('Info');
    visualizationFolder.add(params, 'numHighlightedPoints', 0, 100, 1)
        .name('Number of Highlighted Points')
        .onChange(() => {
            highlightRandomPoints();
        });
    visualizationFolder.add(params, 'trailLength', 10, 200, 1)
        .name('Trail Length');
    // Quantum numbers controls
    quantumFolder.add(params, 'n', 1, 6, 1)
        .name('n: Principal Quantum Number')
        .onChange(() => {
            if (params.l > params.n - 1) {
                params.l = params.n - 1;
            }
            // also make sure m is within bounds
            if (params.m > params.l || params.m < -params.l) {
                params.m = 0;
            }
            updateQuantumNumbers();
            recreateObjects();
        });
    // Spin 1/2 Transformation controls
    spinFolder.add(params, 'spinorMode', ['Off', 'Overlay', 'Spinor Only', 'Charge only'])
        .name('Spin 1/2 Mode')
        .onChange(() => {
            adjustGUIControls();
            recreateObjects();
        });

    spinFolder.add(params, 'maxKernelAngle', 0, 180)
        .name('Max Kernel Angle');

    spinFolder.add(params, 'decayPower', 0.5, 3.0)
        .name('Decay Power');
		
		
	chargeFolder = gui.addFolder('Charge Wave Transformation');
	chargeFolder.closed = false;

    sphericalFolder = gui.addFolder('Scalar Harmonics');
    sphericalFolder.closed = false;

    vectorFolder = gui.addFolder('Vector Harmonics');
    vectorFolder.closed = false;

	chargeFolder.add(params, 'chargeMode', ['Off', 'Plus Charge', 'Minus Charge'])
	  .name('Charge Mode')
	  .onChange(() => {
		recreateObjects();
	  });

    chargeFolder.add(params, 'waveNumber', 1, 10)
        .name('Wave number');



    visualizationFolder.add(params, 'amplitudeScale', 0.1, 10).name('Amplitude Scale');
    sphericalFolder.add(params, 'displacementScale', 0.0, 2).name('Displacement Scale');

    sphericalFolder.add(params, 'omega', 0.1, 10).name('Angular Frequency');
    visualizationFolder.add(params, 'timeScale', 1, 10).name('Time Scale');
    visualizationFolder.add(params, 'autoRotate').name('Animate Wave');
    visualizationFolder.add(params, 'colorScheme', ['Amplitude', 'Phase']).name('Color Scheme');
    visualizationFolder
        .add(params, 'sliceAxis', ['None', 'X', 'Y', 'Z'])
        .name('Slice Axis')
        .onChange(recreateObjects);
    visualizationFolder.add(params, 'slicePosition', -10, 10).name('Slice Position');
    visualizationFolder.add(params, 'sliceWidth', 0.12, 5).name('Slice Width');
    visualizationFolder
        .add(params, 'pointSize', 0.1, 3)
        .name('Point Size')
        .onChange(recreateObjects);
    visualizationFolder
        .add(params, 'extent', 3, 7, 1)
        .name('Grid Extent')
        .onChange(recreateObjects);
    scaleController = visualizationFolder
        .add(params, 'scale', 0.2, 5.0)
        .name('Spherical Scale')
        .onChange(() => {
            precomputeWavefunctionData();
            recreateObjects();
        });

    // Extend GUI controls for Vector Spherical Harmonics (VSH)
    vectorFolder
        .add(params, 'vectorScale', 0.1, 5.0)
        .name('Vector Scale')
        .onChange(recreateObjects);

    vectorFolder
        .add(params, 'skipPoints', 1, 10, 1)
        .name('Skip Points')
        .onChange(recreateObjects);
    vectorFolder
        .add(params, 'threshold', 0.001, 0.1)
        .name('Threshold')
        .onChange(recreateObjects);
    document.getElementById('gui-container').appendChild(gui.domElement);
    document.getElementById('info-close').onclick = function () {
        document.getElementById('info-popup').style.display = 'none';
    };
    visualizationFolder.closed = false;
    adjustGUIControls();



    // Initialize l and m controllers
    updateQuantumNumbers();

}
function smoothStep(edge0, edge1, x) {
  const t = Math.min(Math.max((x - edge0) / (edge1 - edge0), 0.0), 1.0);
  return t * t * (3 - 2 * t);
}
function computeChargeStrength(originalPosition) {
    const radius = originalPosition.length(); // Distance from the center
    const R = params.sphereRadius;           // Sphere's radius
    let maxAngle = (params.maxKernelAngle / 180.0) * Math.PI; // Max angle in radians

    // Change direction based on mode
    if (params.chargeMode === 'Minus Charge') {
        maxAngle = -maxAngle;
    }

	const transitionWidth= R*0.5;
	const innerRadius = R - transitionWidth;
    const outerRadius = R + transitionWidth;
  
    if (radius <= innerRadius) {
		// Fully inside the sphere
		maxAngle = maxAngle * (radius / R);
	  } else if (radius >= outerRadius) {
		// Fully outside the sphere
		maxAngle = maxAngle * (R / radius) ** 2;
	  } else {
		// Transition region: interpolate between inside and outside formulas
		const t = smoothStep(innerRadius, outerRadius, radius);
		const insideAngle = maxAngle * (radius / R);
		const outsideAngle = maxAngle * (R / radius) ** 2;
		maxAngle = (1 - t) * insideAngle + t * outsideAngle;
	  }
	  
	  return maxAngle*0.1;
	}
	
	
// Function to adjust GUI controls based on harmonic type and spinor mode
function adjustGUIControls() {
    if (!params || !visualizationFolder || !spinFolder) return;

     const isVector =
        params.harmonicType === 'Magnetic Vector Harmonics' ||
        params.harmonicType === 'Electric Vector Harmonics';
    const isScalar = params.harmonicType === 'Scalar';
    const isSpinor = params.applySpinHalf || params.harmonicType === 'Spinor Only' || params.spinorMode !== 'Off';
	const isCharge = params.applyChargeTransformation ||params.harmonicType === 'Charge Only' ;
    visualizationFolder.__controllers.forEach(controller => {
        if (controller.property === 'vectorScale' ||
            controller.property === 'skipPoints' ||
            controller.property === 'threshold') {
            controller.domElement.parentElement.style.display = isVector ? '' : 'none';
        }
    });



  if (isCharge && params.applyChargeTransformation == false) {
      // use slicing
        params.sliceAxis = 'Z';
        // also refresh the control

        params.slicePosition = 0;
        params.sliceWidth = 0.12;
        params.amplitudeScale=0.1;
  }
  else if (isSpinor && params.applySpinHalf ==false) {
      // use slicing
        params.sliceAxis = 'Y';
        params.slicePosition = 0;
        params.sliceWidth = 0.12;
        params.amplitudeScale=0.1;
        params.pointSize=1.5;
        params.timeScale=5;
  }
  else if (isScalar) {

      // no slicing
        params.sliceAxis = 'None';
  }
  else if (isVector) {
      // slicig
  }
    visualizationFolder.__controllers.forEach(controller => {
        if (controller.property === 'sliceAxis') {
            controller.updateDisplay();
        }
    });
    spinFolder.domElement.style.display = isSpinor ? '' : 'none';
    chargeFolder.domElement.style.display = isCharge ? '' : 'none';
    sphericalFolder.domElement.style.display = (isVector  || isScalar)? '' : 'none';
    vectorFolder.domElement.style.display = (isVector )? '' : 'none';

    // if scalar or vector, open the quantum folder
    if (isScalar || isVector) {
        quantumFolder.open();
    }
    else {
        quantumFolder.close();
    }
}

window.onclick = function (event) {
    const infoPopup = document.getElementById('info-popup');
    if (event.target === infoPopup) {
        infoPopup.style.display = 'none';
    }
};

// End of guicontrols.js content

// Call the initialization functions after the document has loaded
document.addEventListener('DOMContentLoaded', function() {
  initializeGUI();

  // Start the simulation
  createObjects();
  animate();
});
