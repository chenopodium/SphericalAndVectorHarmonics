// Declare 'params' as a global variable
let params = {};
let torusParams= {};
// Begin waveFunctions.js content
// At the top of your file
const tempVector = new THREE.Vector3();
const tempVector1 = new THREE.Vector3();
const tempVector2 = new THREE.Vector3();
const tempVector3 = new THREE.Vector3();
const tempVector4 = new THREE.Vector3();
let displacementVector1 = new THREE.Vector3(0, 0, 0);
let displacementVector2 = new THREE.Vector3(0, 0, 0);
let displacementVector = new THREE.Vector3(0, 0, 0);
const tempQuaternion1 = new THREE.Quaternion();
const tempQuaternion2 = new THREE.Quaternion();
const tempQuaternion3 = new THREE.Quaternion();
const tempColor = new THREE.Color();
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
let sphericalHarmonicCache = {};

const normalizationCache = {};

function getNormalization(l, m) {
    const key = `${l},${m}`;
    if (normalizationCache[key]) {
        return normalizationCache[key];
    }
    const norm = Math.sqrt(
        ((2 * l + 1) / (4 * Math.PI)) *
        (factorial(l - Math.abs(m)) / factorial(l + Math.abs(m)))
    );
    normalizationCache[key] = norm;
    return norm;
}
function computeSphericalHarmonic(l, m, P_lm, theta, phi, t) {
    // Quantize t and phi to reduce cache size
    // const tKey = t.toFixed(4);
    //  const phiKey = phi.toFixed(4);
    //  const key = `${l},${m},${theta.toFixed(4)},${phiKey},${tKey}`;

    //  if (sphericalHarmonicCache[key]) {
    //      return sphericalHarmonicCache[key];
    //  }

    const normalization = getNormalization(l, m);

    let Y_lm = normalization * P_lm; // Compute the scalar part

    const phase = m * phi - params.omega * t; // Introduce time dependence through phase
    if (m > 0) {
        Y_lm *= Math.cos(phase); // For positive m
    } else if (m < 0) {
        Y_lm *= Math.sin(-phase); // For negative m
    } else {
        Y_lm *= Math.cos(-params.omega * t); // For m = 0
    }

    //  sphericalHarmonicCache[key] = Y_lm;
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
let radialWaveFunctionCache = {};

function radialWaveFunction(n, l, r) {
    const rKey = r.toFixed(4);
    const key = `${n},${l},${rKey}`;

    //  if (radialWaveFunctionCache[key]) {
    //       return radialWaveFunctionCache[key];
    //  }

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

  //  radialWaveFunctionCache[key] = radialPart;
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

    let result;

    if (n === 0) {
        result = Math.sin(x) / x; // j_0(x)
    } else if (n === 1) {
        result = (Math.sin(x) / (x * x)) - (Math.cos(x) / x); // j_1(x)
    } else if (n === -1) {
        result = Math.cos(x) / x; // j_{-1}(x)
    } else if (n > 1) {
        // Use upward recurrence for n > 1
        let jn_minus2 = Math.sin(x) / x; // j_0(x)
        let jn_minus1 = (Math.sin(x) / (x * x)) - (Math.cos(x) / x); // j_1(x)
        for (let l = 2; l <= n; l++) {
            result = ((2 * l - 1) / x) * jn_minus1 - jn_minus2;
            jn_minus2 = jn_minus1;
            jn_minus1 = result;
        }
    } else {
        // Use downward recurrence for n < -1
        // Starting from n = -1
        let jn_plus1 = Math.cos(x) / x; // j_{-1}(x)
        let jn = Math.sin(x) / x; // j_0(x)
        for (let l = -1; l >= n; l--) {
            result = ((2 * l + 1) / x) * jn - jn_plus1;
            jn_plus1 = jn;
            jn = result;
        }
    }

    sphericalBesselCache[key] = result;
    return result;
}


// Derivative of rho * j_n(rho) with respect to rho
let derivativeRhoZnCache = {};
function derivative_rho_zn(n, rho) {
    const key = `${n},${rho.toFixed(4)}`;
    if (derivativeRhoZnCache[key]) {
        return derivativeRhoZnCache[key];
    }

    const jn = sphericalBessel_jn(n, rho);
    let result;

    if (n === 0) {
        // Special case for n = 0
        const jn_minus1 = Math.cos(rho) / rho; // j_{-1}(rho) = cos(rho)/rho
        result = rho * jn_minus1;
    } else {
        const jn_minus1 = sphericalBessel_jn(n - 1, rho);
        result = rho * jn_minus1 - n * jn;
    }

    derivativeRhoZnCache[key] = result;
    return result;
}

let vectorAtPointCache = {};

function computeVectorAtPoint(l, m, data, t, waveNumber) {
    // waveNumber: 1 or 2, indicating which wave function data to use

    //const r = data[`r${waveNumber}`];
   // const theta = data[`theta${waveNumber}`];
    const phi = data[`phi${waveNumber}`];

    const epsilon = 1e-5; // Small value to avoid division by zero
    const maxVectorMagnitude = 1.0; // Maximum allowed vector magnitude

    const rho = data[`rho${waveNumber}`];
    if (rho < epsilon) {
        return new THREE.Vector3(0, 0, 0); // Return zero vector if rho is too small
    }

    // Construct a cache key based on parameters
    //  const tKey = t.toFixed(4);
    //    const key = `${n},${m},${r.toFixed(4)},${theta.toFixed(4)},${phi.toFixed(4)},${tKey},${waveNumber}`;

    //   if (vectorAtPointCache[key]) {
    //       return vectorAtPointCache[key];
    //    }

    // Retrieve precomputed values
    const P_lm = data[`P_lm${waveNumber}`]; // Use precomputed P_lm
    const dP_lm_dtheta = data[`dP_lm_dtheta${waveNumber}`]; // Use precomputed derivative

    // Introduce time dependence through phase
    const phase = m * phi - params.omega * t;

    // Compute angular functions with time dependence
    const cos_m_phi = Math.cos(phase);
    const sin_m_phi = Math.sin(phase);

    let sin_theta = data[`sinTheta${waveNumber}`];
    if (Math.abs(sin_theta) < epsilon) {
        sin_theta = epsilon; // Avoid division by zero
    }

    let N_r = 0,
        N_theta = 0,
        N_phi = 0;

    const jn = data[`jn${waveNumber}`];
    const djn_drho = data[`djn_drho${waveNumber}`];

    if (params.harmonicType === 'Magnetic Vector Harmonics') {
        // Magnetic Vector Spherical Harmonics

        // Radial component is zero
        N_r = 0;

        // Theta component
        N_theta = (jn / sin_theta) * m * P_lm * sin_m_phi;

        // Phi component
        N_phi = -jn * dP_lm_dtheta * cos_m_phi;

    } else if (params.harmonicType === 'Electric Vector Harmonics') {
        // Electric Vector Spherical Harmonics

        // Radial component
        N_r = l * (l + 1) * (jn / rho) * P_lm * cos_m_phi;

        // Theta component
        N_theta = djn_drho * dP_lm_dtheta * cos_m_phi;

        // Phi component
        N_phi = djn_drho * (m * P_lm / sin_theta) * sin_m_phi;

    } else {
        // For other harmonic types, return zero vector
        return new THREE.Vector3(0, 0, 0);
    }

    // Convert to Cartesian coordinates using precomputed trigonometric values
    const sinTheta = data[`sinTheta${waveNumber}`];
    const cosTheta = data[`cosTheta${waveNumber}`];
    const sinPhi = data[`sinPhi${waveNumber}`];
    const cosPhi = data[`cosPhi${waveNumber}`];

    const sin_theta_cos_phi = sinTheta * cosPhi;
    const sin_theta_sin_phi = sinTheta * sinPhi;

    const e_r = tempVector1.set(
        sin_theta_cos_phi,
        sin_theta_sin_phi,
        cosTheta
    );
    const e_theta= tempVector2.set(
        cosTheta * cosPhi,
        cosTheta * sinPhi,
        -sinTheta
    ); // Theta unit vector
    const e_phi = tempVector3.set(-sinPhi, cosPhi, 0); // Phi unit vector

    const N_vector = tempVector4.set(0, 0, 0);
    N_vector.addScaledVector(e_r, N_r);
    N_vector.addScaledVector(e_theta, N_theta);
    N_vector.addScaledVector(e_phi, N_phi);

    // Cap the vector magnitude to prevent excessively large vectors
    N_vector.clampLength(0, maxVectorMagnitude);

    // Cache the result
    //   vectorAtPointCache[key] = N_vector.clone();

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
    radialWaveFunctionCache = {};
    sphericalHarmonicCache = {};
    vectorAtPointCache = {};
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
    if ((params.harmonicType === 'Magnetic Vector Harmonics' || params.harmonicType === 'Electric Vector Harmonics') && (params.l === 0 || params.l2 === 0)) {
        updateUserMessage('No vector field can be visualized for l = 0 in magnetic or electric vector harmonics.');
    } else {
        updateUserMessage(''); // Clear message when not applicable
    }
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
function computeAndVisualizeVectorField() {
    const positionsArray = pointCloud.geometry.attributes.position.array;
    const numPoints = positionsArray.length / 3;

    const pointData = pointCloud.geometry.userData.pointData;
    const n = params.n;
    const l = params.l;
    const m = params.m;
    const n2 = params.n2;
    const l2 = params.l2;
    const m2 = params.m2;
    const showSecondWave = params.showSecondWave;

    arrowGroup = new THREE.Group();
    arrowHelpers = [];
    arrowIndices = [];

    const threshold = params.threshold;
    const scaleFactor = params.vectorScale;

    for (let i = 0; i < numPoints; i++) {
        const data = pointData[i];
        const x = data.x0;
        const y = data.y0;
        const z = data.z0;

        // Determine if the point is within the slicing plane
        let inSlicePlane = false;
        if (params.sliceAxis !== 'None') {
            const sliceValue = params.slicePosition;
            const sliceDelta = params.sliceWidth;
            switch (params.sliceAxis) {
                case 'X':
                    inSlicePlane = Math.abs(x - sliceValue) < sliceDelta;
                    break;
                case 'Y':
                    inSlicePlane = Math.abs(y - sliceValue) < sliceDelta;
                    break;
                case 'Z':
                    inSlicePlane = Math.abs(z - sliceValue) < sliceDelta;
                    break;
            }
        }

        // If the point is not in the slicing plane, apply skipping logic
        if (!inSlicePlane && params.skipPoints > 1) {
            if (i % params.skipPoints !== 0) {
                continue; // Skip this point
            }
        } else {
            if (params.skipPoints > 1) {
                if (
                    (100 * x) % params.skipPoints !== 0 ||
                    (100 * y) % params.skipPoints !== 0 ||
                    (100 * z) % params.skipPoints !== 0
                ) {
                    continue; // Skip this point
                }
            }
        }

        // Compute the vector at this point for the first wave
        const vector1 = computeVectorAtPoint(l, m, data, params.t, 1);

        // Initialize the combined vector
        let vector = vector1.clone();

        // If the second wave is enabled, compute the vector for the second wave and add it
        if (showSecondWave) {
            const vector2 = computeVectorAtPoint(l2, m2, data, params.t, 2);
            vector.add(vector2);
        }

        const magnitude = vector.length();
        let length = magnitude * scaleFactor* params.overallScale ;

        // Optionally scale up arrows in the slicing plane
        if (inSlicePlane) {
            length *= 2; // Adjust this factor as needed
        }

        // Decide whether to create an arrow at this point
        if (length > threshold || inSlicePlane) {
            const dir = vector.clone().normalize();
            const origin = new THREE.Vector3(x, y, z);

            // Compute color based on vector direction
            const phi1 = data.phi1;
            let phi2 = 0;
            if (showSecondWave) {
                phi2 = data.phi2;
            }
            const combinedPhi = (phi1 + phi2) / (showSecondWave ? 2 : 1);
            const normalizedPhi = (combinedPhi + Math.PI) / (2 * Math.PI); // Normalize to [0,1]
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

            // Set visibility based on slicing
            let visible = true;
            if (params.sliceAxis !== 'None') {
                visible = inSlicePlane; // Only show arrows in the slicing plane
            }

            arrowHelper.visible = visible;

            arrowGroup.add(arrowHelper);
            arrowHelpers.push(arrowHelper);
            arrowIndices.push(i); // Store the index corresponding to this arrow
        }
    }

    scene.add(arrowGroup);
}


function transformChargeMultiple(originalPosition, spheres) {
    const transformedPositions = [];
    const scalingFactors = [];

    let i = 0;
    spheres.forEach(sphere => {
        i += 1;
        let chargeMode = params.chargeMode;
        if (i>1) {
            chargeMode = params.chargeMode2;
        }
        // Compute displacement from sphere to point
        const displacement = originalPosition.clone().sub(sphere.center);
        const r = displacement.length();

        // Define wave parameters
        const k =  4 * Math.PI / params.waveNumber;  // Wave number
        const omega = params.omega;              // Angular frequency

        // Determine pulsation direction and calculate scaling factor
        let scalingFactor = 1.0;
        let stretch = 0.5 / (1 + r);

        if (chargeMode === 'Plus Charge') {
            // Outward pulsation
            scalingFactor = 1.0 + stretch * Math.sin(k * r - omega * sphere.phase);
        } else {
            // Inward pulsation
            scalingFactor = 1.0 + stretch * Math.sin(k * r + omega * sphere.phase);
        }

        // Ensure scaling factor remains positive and within bounds
        scalingFactor = Math.max(0.01, Math.min(1.5, scalingFactor));

        // Apply the scaling to the displacement
        const newDisplacement = displacement.clone().multiplyScalar(scalingFactor);

        // Compute the new transformed position
        const newPosition = sphere.center.clone().add(newDisplacement);

        transformedPositions.push(newPosition);
        scalingFactors.push(scalingFactor);
    });

    // Combine the transformed positions by averaging
    const combinedPosition = new THREE.Vector3(0, 0, 0);
    transformedPositions.forEach(pos => {
        combinedPosition.add(pos);
    });
    combinedPosition.divideScalar(spheres.length);

    // Average the scaling factors
    const totalScalingFactor = scalingFactors.reduce((sum, val) => sum + val, 0) / spheres.length;

    return { newPosition: combinedPosition, scalingFactor: totalScalingFactor };
}



function createPointCloud() {
    const waveShift = params.showSecondWave ? params.waveSeparation / 2 : 0;

// Compute grid boundaries based on wave separation
    const gridXMin = -params.extent - waveShift;
    const gridXMax = params.extent + waveShift;

  const gridStep = 0.25; // Distance between grid points
  const positions = []; // Array to store positions of points
  const colors = []; // Array to store colors of points
  const sizes = []; // Array to store sizes of points
  const alphas = []; // Array to store alpha (transparency) of points
  const pointData = []; // Array to store per-point data

  // Generate a 3D grid of points
  for (let x = gridXMin; x <= gridXMax; x += gridStep) {
    for (let y = -params.extent; y <= params.extent; y += gridStep) {
      for (let z = -params.extent; z <= params.extent; z += gridStep) {
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

  const positionsArray = new Float32Array(positions);
  const colorsArray = new Float32Array(colors);
  const sizesArray = new Float32Array(sizes);
  const alphasArray = new Float32Array(alphas);

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(positionsArray, 3));
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


function precomputeWavefunctionData() {
  const pointData = pointCloud.geometry.userData.pointData;
  const n = params.n;
  const l = params.l;
  const m = params.m;
  const n2 = params.n2;
  const l2 = params.l2;
  const m2 = params.m2;
  const scale = params.scale;
  const waveSeparation = params.waveSeparation;
  const showSecondWave = params.showSecondWave;


  const k = 1; // Wave number (can be adjusted)

  for (let i = 0; i < pointData.length; i++) {
    const data = pointData[i];

    const halfWaveSeparation = waveSeparation / 2;
    // Original positions for the first wave function
    const x0 = data.x0- halfWaveSeparation;
    const y0 = data.y0;
    const z0 = data.z0;

    // Shifted positions for the second wave function
    const xShifted = x0 + waveSeparation;
    const yShifted = y0;
    const zShifted = z0;

    // Compute r, theta, phi for the first wave function
    const r1 = Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0);
    const cosTheta1 = y0 / (r1 || 1e-8);
    const theta1 = Math.acos(cosTheta1);
    const phi1 = Math.atan2(z0, x0);

    data.r1 = r1;
    data.theta1 = theta1;
    data.phi1 = phi1;
    data.cosTheta1 = cosTheta1;

    // Compute trigonometric values for the first wave function
    data.sinTheta1 = Math.sin(theta1);
    data.sinPhi1 = Math.sin(phi1);
    data.cosPhi1 = Math.cos(phi1);

    // Compute and store associated Legendre polynomials for the first wave function
    data.P_lm1 = associatedLegendreP(l, m, cosTheta1);

    // Compute and store derivative of P_lm with respect to theta for the first wave function
    data.dP_lm_dtheta1 = derivative_associatedLegendreP_theta(l, m, theta1);

    // Compute and store radial wave functions for the first wave function
    const r_scaled1 = r1 * scale;
    data.R_nl1 = radialWaveFunction(n, l, r_scaled1);

    // Compute and store spherical Bessel functions for the first wave function
    const rho1 = k * r_scaled1;
    data.rho1 = rho1;
    if (rho1 > 1e-5) {
        // For the first wave
        data.jn1 = sphericalBessel_jn(l, rho1);
        data.jn_plus1_1 = sphericalBessel_jn(l + 1, rho1);
        data.djn_drho1 = (l / rho1) * data.jn1 - data.jn_plus1_1;
        data.derivative_rhoz_n1 = derivative_rho_zn(l, rho1);
    } else {
      data.jn1 = 0;
      data.djn_drho1 = 0;
      data.derivative_rhoz_n1 = 0;
    }

    if (showSecondWave) {
      // Compute r, theta, phi for the second wave function
      const r2 = Math.sqrt(
        xShifted * xShifted + yShifted * yShifted + zShifted * zShifted
      );
      const cosTheta2 = yShifted / (r2 || 1e-8);
      const theta2 = Math.acos(cosTheta2);
      const phi2 = Math.atan2(zShifted, xShifted);

      data.r2 = r2;
      data.theta2 = theta2;
      data.phi2 = phi2;
      data.cosTheta2 = cosTheta2;

      // Compute trigonometric values for the second wave function
      data.sinTheta2 = Math.sin(theta2);
      data.sinPhi2 = Math.sin(phi2);
      data.cosPhi2 = Math.cos(phi2);

      // Compute and store associated Legendre polynomials for the second wave function
      data.P_lm2 = associatedLegendreP(l2, m2, cosTheta2);

      // Compute and store derivative of P_lm with respect to theta for the second wave function
      data.dP_lm_dtheta2 = derivative_associatedLegendreP_theta(l2, m2, theta2);

      // Compute and store radial wave functions for the second wave function
      const r_scaled2 = r2 * scale;
      data.R_nl2 = radialWaveFunction(n2, l2, r_scaled2);

      // Compute and store spherical Bessel functions for the second wave function
      const rho2 = k * r_scaled2;
      data.rho2 = rho2;
      if (rho2 > 1e-5) {
          data.jn2 = sphericalBessel_jn(l2, rho2);
          data.jn_plus1_2 = sphericalBessel_jn(l2 + 1, rho2);
          data.djn_drho2 = (l2 / rho2) * data.jn2 - data.jn_plus1_2;
          data.derivative_rhoz_n2 = derivative_rho_zn(l2, rho2);
      } else {
        data.jn2 = 0;
        data.djn_drho2 = 0;
        data.derivative_rhoz_n2 = 0;
      }
    }
  }
}

function spinHalfTransformMultiple(originalPosition, spheres) {
    let transformedPosition = originalPosition.clone();

    let i = 0;
    spheres.forEach(sphere => {
        i += 1;
        let direction =1;
        if (i>1) {
            if (params.spinor2Flipped) direction = -1;
        }
        // Compute timeQuaternion q_s(t)
        const phase_s = direction*sphere.phase; // Each sphere has its own phase
        const timeQuaternion = new THREE.Quaternion();
        timeQuaternion.setFromAxisAngle(
            new THREE.Vector3(0, 1, 0), // Y-axis rotation
            phase_s
        );

        // Inverse quaternion
        const inverseTimeQuaternion = timeQuaternion.clone().conjugate();

        // Compute displacement from sphere to point
        const displacement = transformedPosition.clone().sub(sphere.center);
        const radius = displacement.length();

        let kernelRotationAngle = 0;
        const angle = (params.maxKernelAngle / 180.0) * Math.PI;
        const power = params.decayPower;

        // Original kernel function: Decaying rotation angle with distance
        const adjustedRadius = Math.max(radius - params.sphereRadius, 0);
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

        // Apply the combined rotation to the displacement
        const rotatedDisplacement = displacement.clone().applyQuaternion(combinedQuaternion);

        // Compute the new transformed position
        transformedPosition = sphere.center.clone().add(rotatedDisplacement);
    });

    return { newPosition: transformedPosition };
}

function torusSpinHalfTransformMultiple(originalPosition, torusParams) {
    let transformedPosition = originalPosition.clone();

    const torusAxisNormalized = torusParams.torusAxis.clone().normalize();

    // Compute displacement from torusCenter to point
    const displacement = transformedPosition.clone().sub(torusParams.torusCenter);

    // Compute radius from torus center to point
    const radius = displacement.length();

    // Project displacement onto plane perpendicular to torusAxis
    const displacementOnPlane = displacement.clone().projectOnPlane(torusAxisNormalized);

    // Compute rotation axis for kernelQuaternion (normalized vector in the plane)
    const rotationAxis = displacementOnPlane.clone().normalize();

    // Compute timeQuaternion: rotation around torusAxis by phase
    const timeQuaternion = new THREE.Quaternion();
    timeQuaternion.setFromAxisAngle(
        torusAxisNormalized,
        params.phase
    );

    // Compute inverseTimeQuaternion
    const inverseTimeQuaternion = timeQuaternion.clone().conjugate();

    // Compute kernelRotationAngle based on distance from torus center
    let adjustedRadius = Math.max(radius - torusParams.torusRadius, 0);
    // if we are inside the torus, we linearly interpolate the kernel angle from 0 to maxKernelAngle
    let fraction = 0;
    if (adjustedRadius < torusParams.torusTubeRadius) {
        fraction = adjustedRadius / torusParams.torusTubeRadius;

    }
    else fraction = 1 / Math.pow(adjustedRadius + 1, params.decayPower);
    const kernelRotationAngle = (params.maxKernelAngle / 180.0) * Math.PI * fraction;

    // Compute kernelQuaternion: rotation around rotationAxis by kernelRotationAngle
    const kernelQuaternion = new THREE.Quaternion();
    kernelQuaternion.setFromAxisAngle(
        rotationAxis,
        kernelRotationAngle
    );

    // Combined rotation: q(t) * k(r) * q⁻¹(t)
    const combinedQuaternion = new THREE.Quaternion();
    combinedQuaternion.multiplyQuaternions(timeQuaternion, kernelQuaternion);
    combinedQuaternion.multiply(inverseTimeQuaternion);

    // Apply combinedQuaternion to displacement
    const rotatedDisplacement = displacement.clone().applyQuaternion(combinedQuaternion);

    // Compute the new transformed position
    transformedPosition = torusParams.torusCenter.clone().add(rotatedDisplacement);


    if ( originalPosition.x ==3 && originalPosition.y == 0 && originalPosition.z == 0) {
        p("fraction = " + fraction + ", kernelRotationAngle = " + kernelRotationAngle + ", adjustedRadius = " + adjustedRadius +
            ", torusParams.decayPower = " + params.decayPower);
        p("originalPosition = " + originalPosition.x + ", " + originalPosition.y + ", " + originalPosition.z);
        p("transformedPosition = " + transformedPosition.x + ", " + transformedPosition.y + ", " + transformedPosition.z);
    }
    return { newPosition: transformedPosition };
}
function p(s) {
    console.log (s);
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
    const n2 = params.n2;
    const l2 = params.l2;
    const m2 = params.m2;
    const pointData = pointCloud.geometry.userData.pointData;
    const waveSeparation = params.waveSeparation;
    const showSecondWave = params.showSecondWave;

    // Define wave centers
    const center1 = new THREE.Vector3(-waveSeparation/2, 0, 0);
    const center2 = new THREE.Vector3(waveSeparation/2, 0, 0);

    const waveCenters = [];

    const phase1 = params.phase;
    waveCenters.push({ center: center1, phase: phase1 });


    if (showSecondWave) {
        const phase2 = params.phase; // If you have a separate phase for the second wave, adjust here
        waveCenters.push({ center: center2, phase: phase2 });
    }
    // Time-based rotation quaternion q(t) around Y-axis for spin 1/2
    const timeQuaternion = new THREE.Quaternion();
    timeQuaternion.setFromAxisAngle(
        new THREE.Vector3(0, 1, 0), // Rotation axis (Y-axis)
        params.phase // Rotation angle (phase)
    );

    for (let i = 0; i < pointCount; i++) {
        const data = pointData[i];
        const x0 = data.x0;
        const y0 = data.y0;
        const z0 = data.z0;

        // Original positions
        let positionOfPoint = data.originalPosition.clone();

        // Variables for scaling factors (used in Charge Only harmonic type)
        let scalingFactor1 = 1.0;
        let scalingFactor2 = 1.0;

        // Apply transformations to positionOfPoint
        if (params.spinorMode !== 'Off' || params.harmonicType === 'Spinor Only') {
            const transformResult = spinHalfTransformMultiple(positionOfPoint, waveCenters);
	   //      const transformResult = torusSpinHalfTransformMultiple(positionOfPoint, torusParams);
    
            positionOfPoint = transformResult.newPosition;
        }
        if (params.applyChargeTransformation || params.harmonicType === 'Charge Only') {
            const transformResult =
                transformChargeMultiple(positionOfPoint, waveCenters);
            positionOfPoint = transformResult.newPosition;
            scalingFactor1 = transformResult.scalingFactor;
        }

        // Initialize amplitude and displacement vectors
        let amplitude1 = 0,
            amplitude2 = 0;

        let amplitude = 0;
        // set displacement vectors to zero
        displacementVector1.x = 0;
        displacementVector1.y = 0;
        displacementVector1.z = 0;
        displacementVector2.x = 0;
        displacementVector2.y = 0;
        displacementVector2.z = 0;
        displacementVector.x = 0;
        displacementVector.y = 0;
        displacementVector.z = 0;

     //   let color = new THREE.Color();
        let size = params.pointSize;

        // Now, handle each harmonic type
        if (params.harmonicType === 'Scalar') {
            // First wave
            let Y_lm1 = computeSphericalHarmonic(
                l,
                m,
                data.P_lm1,
                data.theta1,
                data.phi1,
                params.t
            );
            const Psi_nlm1 = data.R_nl1 * Y_lm1;
            amplitude1 = Math.abs(Psi_nlm1) * 50; // Scale amplitude

            if (isNaN(amplitude1) || !isFinite(amplitude1)) {
                amplitude1 = 0;
            }

            // Displacement for the first wave
            const displacementMagnitude1 = amplitude1 * params.displacementScale * 2;
            const r1 = positionOfPoint.clone().sub(center1).length(); // Distance from center1
            let unitRadialVector1 = positionOfPoint.clone().sub(center1).normalize();
            if (r1 > 1e-5) {
                unitRadialVector1.divideScalar(r1);
            } else {
                unitRadialVector1.set(0, 0, 0);
            }
            displacementVector1 = unitRadialVector1.multiplyScalar(displacementMagnitude1);

            // Second wave
            if (showSecondWave) {
                let Y_lm2 = computeSphericalHarmonic(
                    l2,
                    m2,
                    data.P_lm2,
                    data.theta2,
                    data.phi2,
                    params.t
                );
                const Psi_nlm2 = data.R_nl2 * Y_lm2;
                amplitude2 = Math.abs(Psi_nlm2) * 50; // Scale amplitude

                if (isNaN(amplitude2) || !isFinite(amplitude2)) {
                    amplitude2 = 0;
                }

                // Displacement for the second wave
                const displacementMagnitude2 = amplitude2 * params.displacementScale * 2;
                const r2 = positionOfPoint.clone().sub(center2).length(); // Distance from center2
                let unitRadialVector2 = positionOfPoint.clone().sub(center2).normalize();
                if (r2 > 1e-5) {
                    unitRadialVector2.divideScalar(r2);
                } else {
                    unitRadialVector2.set(0, 0, 0);
                }
                displacementVector2 = unitRadialVector2.multiplyScalar(displacementMagnitude2);
            }

            // Combine amplitudes and displacements
            amplitude = amplitude1 + amplitude2;
            displacementVector.copy(displacementVector1);
            if (showSecondWave) {
                displacementVector.add(displacementVector2);
            }

            // Color based on amplitude
            if (params.colorScheme === 'Amplitude') {
                const hue = 0.7 - amplitude * 0.7; // Hue from blue (0.7) to red (0.0)
                const lightness = 0.4 + amplitude * 0.6; // Lightness from 0.4 to 1.0
                tempColor.setHSL(hue, 1, lightness);
            } else if (params.colorScheme === 'Phase') {
                // For phase coloring, combine phases from both waves
                const phase1 = m * data.phi1 - params.omega * params.t;
                let phase2 = 0;
                if (showSecondWave) {
                    phase2 = m2 * data.phi2 - params.omega * params.t;
                }
                const combinedPhase = (phase1 + phase2) / (showSecondWave ? 2 : 1);
                const phaseValue =
                    (Math.atan2(Math.sin(combinedPhase), Math.cos(combinedPhase)) + Math.PI) /
                    (2 * Math.PI); // Normalize to [0,1]
                tempColor.setHSL(phaseValue, 1, 0.5); // Hue based on phase
            }

            size =  params.overallScale *params.pointSize * (1 + amplitude * params.amplitudeScale);
            if (isNaN(size) || size <= 0) {
                size = params.pointSize;
            }
        } else if (params.harmonicType === 'Spinor Only') {
            // Compute amplitude and color based on positions after spin 1/2 transformation
            const r1 = positionOfPoint.clone().sub(center1).length();
            amplitude1 = r1 * 0.3;
            let normAmp1 = (amplitude1 + params.extent) / (2 * params.extent);

            if (showSecondWave) {
                const r2 = positionOfPoint.clone().sub(center2).length();
                amplitude2 = r2 * 0.3;
                let normAmp2 = (amplitude2 + params.extent) / (2 * params.extent);
                amplitude = (normAmp1 + normAmp2) / 2; // Average amplitudes
            } else {
                amplitude = normAmp1;
            }

            if (params.colorScheme === 'Amplitude') {
                tempColor.setHSL(amplitude, 1, 0.5);
            } else if (params.colorScheme === 'Phase') {
                const phase1 = m * data.phi1 - params.omega * params.t;
                let phase2 = 0;
                if (showSecondWave) {
                    phase2 = m2 * data.phi2 - params.omega * params.t;
                }
                const combinedPhase = (phase1 + phase2) / (showSecondWave ? 2 : 1);
                const phaseValue =
                    (Math.atan2(Math.sin(combinedPhase), Math.cos(combinedPhase)) + Math.PI) /
                    (2 * Math.PI); // Normalize to [0,1]
                tempColor.setHSL(phaseValue, 1, 0.5); // Hue based on phase
            }
            size = params.pointSize * 0.5;
        } else if (params.harmonicType === 'Charge Only') {
            // Color has been computed during transformations
            const minScale = 0.75; // Adjust this range as necessary
            const maxScale = 1.25;

            let normalizedScale1 = (scalingFactor1 - minScale) / (maxScale - minScale);
            const clampedScale1 = Math.max(0, Math.min(1, normalizedScale1));

            if (showSecondWave) {
                let normalizedScale2 = (scalingFactor2 - minScale) / (maxScale - minScale);
                const clampedScale2 = Math.max(0, Math.min(1, normalizedScale2));

                amplitude = (clampedScale1 + clampedScale2) / 2; // Average amplitudes

                // Set color based on combined scaling factors
                const combinedScalingFactor = (scalingFactor1 + scalingFactor2) / 2;
                const normalizedCombinedScale =
                    (combinedScalingFactor - minScale) / (maxScale - minScale);
                const clampedCombinedScale = Math.max(0, Math.min(1, normalizedCombinedScale));
                const hue = clampedCombinedScale * 0.7;
                tempColor.setHSL(hue, 1.0, 0.5);
            } else {
                amplitude = clampedScale1;
                const hue = clampedScale1 * 0.7;
                tempColor.setHSL(hue, 1.0, 0.5);
            }
            if (params.colorScheme === 'Phase') {

                const phase1 = m * data.phi1 - params.omega * params.t;
                let phase2 = 0;
                if (showSecondWave) {
                    phase2 = m2 * data.phi2 - params.omega * params.t;
                }
                const combinedPhase = (phase1 + phase2) / (showSecondWave ? 2 : 1);
                const phaseValue =
                    (Math.atan2(Math.sin(combinedPhase), Math.cos(combinedPhase)) + Math.PI) /
                    (2 * Math.PI); // Normalize to [0,1]
                tempColor.setHSL(phaseValue, 1, 0.5); // Hue based on phase
            }

            size = 0.5 * params.pointSize * (1 + amplitude * params.amplitudeScale);
        } else if (
            params.harmonicType === 'Magnetic Vector Harmonics' ||
            params.harmonicType === 'Electric Vector Harmonics'
        ) {
            // Compute vector harmonics
            // First wave
            const vector1 = computeVectorAtPoint(l, m, data, params.t, 1);
            const vectorDisplacement1 = vector1.clone().multiplyScalar(
                params.displacementScale * 2*params.overallScale
            );
            amplitude1 = vectorDisplacement1.length() * 5;
            displacementVector1.copy(vectorDisplacement1);

            if (showSecondWave) {
                // Second wave
                const vector2 = computeVectorAtPoint(l2, m2, data, params.t, 2);
                const vectorDisplacement2 = vector2.clone().multiplyScalar(
                    params.displacementScale * 2* params.overallScale
                );
                amplitude2 = vectorDisplacement2.length() * 5;
                displacementVector2.copy(vectorDisplacement2);
            }

            // Combine amplitudes and displacements
            amplitude = amplitude1 + amplitude2;
            displacementVector.copy(displacementVector1);
            if (showSecondWave) {
                displacementVector.add(displacementVector2);
            }

            // Color based on vector properties (e.g., direction)
            const phi1 = data.phi1;
            let phi2 = 0;
            if (showSecondWave) {
                phi2 = data.phi2;
            }
            const combinedPhi = (phi1 + phi2) / (showSecondWave ? 2 : 1);
            const normalizedPhi = (combinedPhi + Math.PI) / (2 * Math.PI); // Normalize to [0,1]
            const hue = normalizedPhi; // Hue varies from 0 to 1 based on phi
            tempColor.setHSL(hue, 1, 0.5);
            if (params.colorScheme === 'Phase') {
                const phase1 = m * data.phi1 - params.omega * params.t;
                let phase2 = 0;
                if (showSecondWave) {
                    phase2 = m2 * data.phi2 - params.omega * params.t;
                }
                const combinedPhase = (phase1 + phase2) / (showSecondWave ? 2 : 1);
                const phaseValue =
                    (Math.atan2(Math.sin(combinedPhase), Math.cos(combinedPhase)) + Math.PI) /
                    (2 * Math.PI); // Normalize to [0,1]
                tempColor.setHSL(phaseValue, 1, 0.5); // Hue based on phase
            }
            size =  params.overallScale * 0.5 * params.pointSize * (1 + amplitude * params.amplitudeScale);
            if (isNaN(size) || size <= 0) {
                size = params.pointSize;
            }
        }

        // Limit the size of the displacement vector
        displacementVector.clampLength(0, 1.0);

        // Compute the final position
        let finalPosition = positionOfPoint.clone().add(displacementVector1);
        if (showSecondWave) {
            const secondPosition = positionOfPoint.clone().add(displacementVector2);
            finalPosition.add(secondPosition).multiplyScalar(0.5); // Average positions
        }

        positionAttribute.setXYZ(i, finalPosition.x, finalPosition.y, finalPosition.z);

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
        if (
            params.harmonicType === 'Spinor Only' ||
            params.harmonicType === 'Charge Only'
        ) {
            alpha = visible ? 0.3 : 0;
        } else {
            alpha = visible ? Math.min(1.0, amplitude - 0.1) : 0;
        }
        if (isNaN(alpha) || alpha < 0 || alpha > 1) {
            alpha = 0;
        }
        alphaAttribute.setX(i, alpha);

        // Set color for points
        colorAttribute.setXYZ(i, tempColor.r, tempColor.g, tempColor.b);
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
function updateArrowHelpers() {
    const scaleFactor = params.vectorScale;
    const pointData = pointCloud.geometry.userData.pointData;
    const n = params.n;
    const l = params.l;
    const m = params.m;
    const n2 = params.n2;
    const l2 = params.l2;
    const m2 = params.m2;
    const showSecondWave = params.showSecondWave;

    for (let i = 0; i < arrowHelpers.length; i++) {
        const arrowHelper = arrowHelpers[i];
        const index = arrowIndices[i];

        const data = pointData[index];
        const x = data.x0;
        const y = data.y0;
        const z = data.z0;

        // Compute the vector at this point with current time
        const vector1 = computeVectorAtPoint(l, m, data, params.t, 1);

        // Initialize the combined vector
        let vector = vector1.clone();

        // If the second wave is enabled, compute the vector for the second wave and add it
        if (showSecondWave) {
            const vector2 = computeVectorAtPoint(l2, m2, data, params.t, 2);
            vector.add(vector2);
        }

        const magnitude = vector.length();
        const length = magnitude * scaleFactor;

        // Update arrow direction and length
        if (length > params.threshold) {
            const dir = tempVector1.copy(vector).normalize();
            const origin = tempVector2.set(x, y, z);

            arrowHelper.position.copy(origin);
            arrowHelper.setDirection(dir);
            arrowHelper.setLength(length);

            // Compute color based on vector direction (phi angle)
            const phi1 = data.phi1;
            let phi2 = 0;
            if (showSecondWave) {
                phi2 = data.phi2;
            }
            const combinedPhi = (phi1 + phi2) / (showSecondWave ? 2 : 1);
            const normalizedPhi = (combinedPhi + Math.PI) / (2 * Math.PI); // Normalize to [0,1]
            const hue = normalizedPhi; // Hue varies from 0 to 1 based on phi
            const color = new THREE.Color();
            color.setHSL(hue, 1, 0.5);

            // Update arrow color
            if (arrowHelper.setColor) {
                arrowHelper.setColor(color);
            } else {
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
    const omega = Math.max(params.omega, 1e-6);
  // Update time variable
    if (params.autoRotate) {
        const dt = 0.01 * params.timeScale;
        params.t += dt;

        // if we are once around reset t to 0
        if (params.t > 2 * Math.PI / omega) {
            params.t = 0;
        }
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
let quantumGui;
let  visualizationFolder, spinFolder, chargeFolder, chargeFolder2, sphericalFolder, vectorFolder;
let experimentalFolder, markerFolder, sliceFolder;
let m1Controller, l1Controller, m2Controller, l2Controller;
let scaleController;

// Function to update m options based on l for Wave 1
function updateMOptions() {
    if (!params) return; // Ensure params is initialized

    const mOptions = [];
    for (let i = -params.l; i <= params.l; i++) {
        mOptions.push(i);
    }

    if (m1Controller) {
        quantumFolder1.remove(m1Controller);
    }
    m1Controller = quantumFolder1
        .add(params, 'm', mOptions)
        .name('m1: Magnetic Quantum Number');
    m1Controller.onChange((value) => {
        params.m = parseInt(value, 10);
        recreateObjects();
    });
}

// Function to update m2 options based on l2 for Wave 2
function updateM2Options() {
    if (!params) return; // Ensure params is initialized

    const m2Options = [];
    for (let i = -params.l2; i <= params.l2; i++) {
        m2Options.push(i);
    }

    if (m2Controller) {
        quantumFolder2.remove(m2Controller);
    }
    m2Controller = quantumFolder2
        .add(params, 'm2', m2Options)
        .name('m2: Magnetic Quantum Number');
    m2Controller.onChange((value) => {
        params.m2 = parseInt(value, 10);
        recreateObjects();
    });
}

// Function to update quantum numbers for Wave 1
function updateQuantumNumbers() {
    if (!params || !quantumFolder1) return; // Ensure params and quantumFolder1 are initialized

    const lOptions = [];
    for (let i = 0; i <= params.n - 1; i++) {
        lOptions.push(i);
    }

    if (l1Controller) {
        quantumFolder1.remove(l1Controller);
    }
    l1Controller = quantumFolder1
        .add(params, 'l', lOptions)
        .name('l1: Azimuthal Quantum Number');
    l1Controller.onChange((value) => {
        params.l = parseInt(value, 10);
        params.l = Math.max(0, Math.min(params.l, params.n - 1));

        if (params.m > params.l || params.m < -params.l) {
            params.m = 0;
        }
        updateMOptions();
        recreateObjects();
    });

    updateMOptions();
}

// Function to update quantum numbers for Wave 2
function updateQuantumNumbers2() {
    if (!params || !quantumFolder2) return; // Ensure params and quantumFolder2 are initialized

    const l2Options = [];
    for (let i = 0; i <= params.n2 - 1; i++) {
        l2Options.push(i);
    }

    if (l2Controller) {
        quantumFolder2.remove(l2Controller);
    }
    l2Controller = quantumFolder2
        .add(params, 'l2', l2Options)
        .name('l2: Azimuthal Quantum Number');
    l2Controller.onChange((value) => {
        params.l2 = parseInt(value, 10);
        if (params.m2 > params.l2 || params.m2 < -params.l2) {
            params.m2 = 0;
        }
        updateM2Options();
        recreateObjects();
    });

    updateM2Options();
}
// Function to highlight a random point
function highlightRandomPoint() {
    if (!pointCloud) return; // Ensure pointCloud is available

    const pointCount = pointCloud.geometry.attributes.position.count;
    highlightedPointIndex = Math.floor(Math.random() * pointCount);

}

function updateUserMessage(message) {
    const messageElement = document.getElementById('user-message'); // Ensure this div exists in your HTML
    if (message) {
        messageElement.style.display = 'block';
        messageElement.innerText = message;
    } else {
        messageElement.style.display = 'none';
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
    n2: 3,
    l2: 1,
    m2: 1,
    showSecondWave: false,
    waveSeparation: 0,
    amplitudeScale: 2.0,
    autoRotate: true,
    colorScheme: 'Amplitude',
    omega: 1.0,
    timeScale: 5.0,
    t: 0,
    sliceAxis: 'None',
    slicePosition: 0.0,
    sliceWidth: 0.12,
    pointSize: 1.0,
    vectorScale: 2,
    harmonicType: 'Scalar',
    skipPoints: 3,
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
      chargeMode2: 'Plus Charge',
    waveNumber: 4,
    maxKernelAngle: 180,
      spinor2Flipped: true,
    decayPower: 2,
    phase: 0.0,
    sphereRadius: 0.1,
      overallScale:1,
  };

  gui = new dat.GUI({ autoPlace: false });
   quantumGui = new dat.GUI({ autoPlace: false });

    // Group GUI controls into folders to reduce clutter

  visualizationFolder = gui.addFolder('Visualization');

  // Visualization controls
  visualizationFolder
    .add(params, 'harmonicType', [
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


    let waveSeparationController = visualizationFolder
        .add(params, 'waveSeparation', 0.0, 2 * params.extent, 0.1)
        .name('Wave Separation')
        .onChange(() => {
            recreateObjects();
        });



    experimentalFolder = quantumGui.addFolder('Experimental');
    experimentalFolder
        .add(params, 'applySpinHalf')
        .name('Add Spin 1/2')
        .onChange(() => {
            adjustGUIControls();
            recreateObjects();
        });

    experimentalFolder
        .add(params, 'applyChargeTransformation')
        .name('Add Charge')
        .onChange(() => {
            adjustGUIControls();
            recreateObjects();
        });

    experimentalFolder.closed = true;

    quantumFolder1 = quantumGui.addFolder('Quantum Numbers Wave 1');
    quantumFolder2 = quantumGui.addFolder('Quantum Numbers Wave 2');

    spinFolder = quantumGui.addFolder('Spin 1/2 Transformation');
    spinFolder.closed = false;


  markerFolder = gui.addFolder('Markers');
    markerFolder
        .add(params, 'numHighlightedPoints', 0, 100, 1)
        .name('Number of Highlighted Points')
        .onChange(() => {
            highlightRandomPoints();
        });

    markerFolder.add(params, 'trailLength', 10, 200, 1).name('Trail Length');
    markerFolder.closed = false;

    // Info button
    gui.add(params, 'showInfo').name('Info');



    visualizationFolder.add(params, 'showSecondWave').name('Second Wave').onChange(() => {
      // change the waveseperation value
        if (params.showSecondWave) {
            params.waveSeparation =  params.extent*1.5;
			params.scale = 3;
        }
        else {
            params.waveSeparation = 0; // Optionally reset when second wave is off
        }
       scaleController.updateDisplay();

      waveSeparationController.updateDisplay();
      adjustGUIControls();
    recreateObjects();
  });

    // Quantum numbers controls for Wave 1
    quantumFolder1
        .add(params, 'n', 1, 6, 1)
        .name('n1: Principal Quantum Number')
        .onChange(() => {
            if (params.l > params.n - 1) {
                params.l = params.n - 1;
            }
            if (params.m > params.l || params.m < -params.l) {
                params.m = 0;
            }
            updateQuantumNumbers();
            recreateObjects();
        });

    // Initialize l1 and m1 controllers
    updateQuantumNumbers();

    // Quantum numbers controls for Wave 2
    quantumFolder2
        .add(params, 'n2', 1, 6, 1)
        .name('n2: Principal Quantum Number')
        .onChange(() => {
            if (params.l2 > params.n2 - 1) {
                params.l2 = params.n2 - 1;
            }
            if (params.m2 > params.l2 || params.m2 < -params.l2) {
                params.m2 = 0;
            }
            updateQuantumNumbers2();
            recreateObjects();
        });

    // Initialize l2 and m2 controllers
    updateQuantumNumbers2();

  // Spin 1/2 Transformation controls
  spinFolder
    .add(params, 'spinorMode', ['Off', 'Overlay', 'Spinor Only', 'Charge only'])
    .name('Spin 1/2 Mode')
    .onChange(() => {
      adjustGUIControls();
      recreateObjects();
    });

  spinFolder.add(params, 'maxKernelAngle', 0, 180).name('Max Kernel Angle');
    spinFolder.add(params, 'spinor2Flipped', 0, 180).name('Flip Spinor 2');
  spinFolder.add(params, 'decayPower', 0.5, 3.0).name('Decay Power');

  chargeFolder = quantumGui.addFolder('Charge Wave 1');
  chargeFolder.closed = true;
    chargeFolder2 = quantumGui.addFolder('Charge Wave 2');
    chargeFolder2.closed = true;

  sphericalFolder = quantumGui.addFolder('Scalar Harmonics');
  sphericalFolder.closed = false;

  vectorFolder = quantumGui.addFolder('Vector Harmonics');
  vectorFolder.closed = false;

  chargeFolder
    .add(params, 'chargeMode', ['Off', 'Plus Charge', 'Minus Charge'])
    .name('Charge Mode 1')
    .onChange(() => {
      recreateObjects();
    });
    chargeFolder2
        .add(params, 'chargeMode2', ['Off', 'Plus Charge', 'Minus Charge'])
        .name('Charge Mode 2')
        .onChange(() => {
            recreateObjects();
        });
  chargeFolder.add(params, 'waveNumber', 1, 10).name('Wave number');

  visualizationFolder
    .add(params, 'amplitudeScale', 0.1, 10)
    .name('Amplitude Scale');
  sphericalFolder
    .add(params, 'displacementScale', 0.0, 2)
    .name('Displacement Scale');

  sphericalFolder.add(params, 'omega', 0.1, 10).name('Angular Frequency');
  visualizationFolder
    .add(params, 'timeScale', 1, 20)
    .name('Time Scale');
  visualizationFolder.add(params, 'autoRotate').name('Animate Wave');
  visualizationFolder
    .add(params, 'colorScheme', ['Amplitude', 'Phase'])
    .name('Color Scheme');



  visualizationFolder
    .add(params, 'pointSize', 0.1, 3)
    .name('Point Size')
    .onChange(recreateObjects);

   let extentController = visualizationFolder
        .add(params, 'extent', 3, 7, 1)
        .name('Grid Extent')
        .onChange(() => {
            // Update the max value of waveSeparationController
            waveSeparationController.max(2 * params.extent);
            recreateObjects();
        });

  scaleController = visualizationFolder
    .add(params, 'scale', 1, 10.0)
    .name('Zoomout')
    .onChange(() => {
      precomputeWavefunctionData();
      recreateObjects();
    });

    sliceFolder = visualizationFolder.addFolder('3D Slicing');
    sliceFolder.closed = false;

    sliceFolder
        .add(params, 'sliceAxis', ['None', 'X', 'Y', 'Z'])
        .name('Slice Axis')
        .onChange(recreateObjects);
    sliceFolder
        .add(params, 'slicePosition', -10, 10)
        .name('Slice Position');
    sliceFolder
        .add(params, 'sliceWidth', 0.12, 5)
        .name('Slice Width');

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
    document.getElementById('quantum-gui-container').appendChild(quantumGui.domElement);
    // move th quantum gui to the right
    quantumGui.domElement.style.position = 'absolute';
    quantumGui.domElement.style.right = '0px';

    // Synchronize the close behavior of both GUI controls
    gui.__closeButton.addEventListener('click', () => {
        if ( quantumGui.domElement.style.display==='')  quantumGui.domElement.style.display = 'none';
        else {
            // show again
            quantumGui.domElement.style.display = '';
        }
    });


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
	
	function adjustGUIControls() {
  if (!params || !visualizationFolder || !spinFolder) return;

  const isVector =
    params.harmonicType === 'Magnetic Vector Harmonics' ||
    params.harmonicType === 'Electric Vector Harmonics';
  const isScalar = params.harmonicType === 'Scalar';
  const isSpinor =
    params.applySpinHalf ||
    params.harmonicType === 'Spinor Only' ||
    params.spinorMode !== 'Off';
  const isCharge =
    params.applyChargeTransformation ||
    params.harmonicType === 'Charge Only';

  visualizationFolder.__controllers.forEach((controller) => {
    if (
      controller.property === 'vectorScale' ||
      controller.property === 'skipPoints' ||
      controller.property === 'threshold'
    ) {
      controller.domElement.parentElement.style.display = isVector ? '' : 'none';
    }
  });

  params.overallScale = 1;
  if (isCharge && params.applyChargeTransformation == false) {
    params.sliceAxis = 'Z';
    params.slicePosition = 0;
    params.sliceWidth = 0.12;
    params.amplitudeScale = 0.1;
  } else if (isSpinor && params.applySpinHalf == false) {
    params.sliceAxis = 'Y';
    params.slicePosition = 0;
    params.sliceWidth = 0.12;
    params.amplitudeScale = 0.1;
    params.pointSize = 1.5;
    params.timeScale = 5;
  } else if (isScalar) {
    params.sliceAxis = 'None';
  }
  else if (isVector) {
    params.sliceAxis = 'Z';
    params.timeScale = 5;
    if (params.m ===0) {
        params.overallScale =2;
    }

  }

  visualizationFolder.__controllers.forEach((controller) => {
    if (controller.property === 'sliceAxis') {
      controller.updateDisplay();
    }
  });
  spinFolder.domElement.style.display = isSpinor ? '' : 'none';
  chargeFolder.domElement.style.display = isCharge ? '' : 'none';
  chargeFolder2.domElement.style.display = isCharge ? '' : 'none';

  sphericalFolder.domElement.style.display = isScalar ? '' : 'none';
  vectorFolder.domElement.style.display = isVector ? '' : 'none';
  experimentalFolder.domElement.style.display = params.showSecondWave ? 'none' : '';
  quantumFolder1.domElement.style.display = isScalar || isVector ? '' : 'none';
    quantumFolder2.domElement.style.display = (isScalar || isVector ) && (params.showSecondWave)? '' : 'none';

  if (isScalar || isVector) {
    quantumFolder1.open();
      quantumFolder2.open();
  } else {
    quantumFolder1.close();
      quantumFolder2.close();
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
