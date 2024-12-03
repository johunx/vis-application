/**
 * A Starting Template for Lab in Vis Applications course module in TNM093
 * -------------------------------------
 *
 * IMPORTANT:
 * - This is a basic template serving as a starting template and NOT intended to cover all requirements.
 * - You are encouraged to implement the lab in your own way.
 * - Feel free to ignore this template if you prefer to start from scratch.
 *
 */

// Main simulation logic

// Select the SVG container
const svg = d3.select("#simulation-area");
const width = svg.attr("width");
const height = svg.attr("height");

// examples of default settings that can be changed later
let rows = parseInt(document.getElementById("rows").value, 10);
let cols = parseInt(document.getElementById("cols").value, 10);
let restoreForce = parseFloat(document.getElementById("restore-force").value);
let damping = parseFloat(document.getElementById("damping").value);
const nodeRadius = 5;
const timeStep = 0.016;
const padding = 50;
const nodeMass = 0.2; // Mass of each node (kg)

// Structural spring parameters
const structuralSpringK = 20; // Structural spring stiffness (kg/s^2)
const structuralSpringB = 0.1; // Structural spring damping (kg/s)
const structuralRestLength = 100; // Structural spring rest length (pixels)

// Shear spring parameters
const shearSpringK = 7; // Shear spring stiffness (kg/s^2)
const shearSpringB = 0.05; // Shear spring damping (kg/s)
const shearRestLength = structuralRestLength * Math.sqrt(2); // Shear spring rest length (pixels)

// Arrays to hold positions, velocities, and forces
let positions = [];
let oldPositions = [];
let velocities = [];
let forces = [];
let isRunning = false;

let useAltForceMethod = false;

/**
 * Initialize the grid with nodes and reset their positions, velocities, and forces.
 */
function initializeGrid() {
  positions = [];
  oldPositions = []; // Initialize old positions array
  velocities = [];
  forces = [];
  const xStep = (width - 2 * padding) / (cols - 1);
  const yStep = (height - 2 * padding) / (rows - 1);

  for (let i = 0; i < rows; i++) {
    const positionRow = [];
    const oldPositionRow = [];
    const velocityRow = [];
    const forceRow = [];
    for (let j = 0; j < cols; j++) {
      const initialPos = [padding + j * xStep, padding + yStep * i];
      positionRow.push(initialPos); // ! TODO: think about how to calculate initial positions for the nodes
      oldPositionRow.push([...initialPos]); // Copy initial position for old positions
      velocityRow.push([0, 0]); // Initial velocity
      forceRow.push([0, 0]); // Initial force
    }
    positions.push(positionRow);
    oldPositions.push(oldPositionRow);
    velocities.push(velocityRow);
    forces.push(forceRow);
  }

  console.log(positions[0]);

  drawNodes();
  drawEdges();
}

/**
 * Draw the nodes (circles) on the SVG.
 */
function drawNodes() {
  // example of how to draw nodes on the svg
  const nodes = svg.selectAll("circle").data(positions.flat());
  nodes
    .enter()
    .append("circle")
    .attr("r", nodeRadius)
    .merge(nodes)
    .attr("cx", (d) => d[0])
    .attr("cy", (d) => d[1])
    .attr("fill", "blue")
    .attr("stroke", "white")
    .attr("stroke-width", 2);

  nodes.exit().remove();
}

/**
 * Draw the edges (lines) connecting the nodes.
 */
function drawEdges() {
  const edges = [];

  // Create horizontal edges
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols - 1; j++) {
      edges.push({
        x1: positions[i][j][0],
        y1: positions[i][j][1],
        x2: positions[i][j + 1][0],
        y2: positions[i][j + 1][1],
        isShear: false, // Add this flag for structural springs
      });
    }
  }

  // Create vertical edges
  for (let i = 0; i < rows - 1; i++) {
    for (let j = 0; j < cols; j++) {
      edges.push({
        x1: positions[i][j][0],
        y1: positions[i][j][1],
        x2: positions[i + 1][j][0],
        y2: positions[i + 1][j][1],
        isShear: false, // Add this flag for structural springs
      });
    }
  }

  // Add diagonal (shear) springs
  for (let i = 0; i < rows - 1; i++) {
    for (let j = 0; j < cols - 1; j++) {
      // Diagonal from top-left to bottom-right
      edges.push({
        x1: positions[i][j][0],
        y1: positions[i][j][1],
        x2: positions[i + 1][j + 1][0],
        y2: positions[i + 1][j + 1][1],
        isShear: true,
      });

      // Add this part for diagonal from top-right to bottom-left
      edges.push({
        x1: positions[i][j + 1][0],
        y1: positions[i][j + 1][1],
        x2: positions[i + 1][j][0],
        y2: positions[i + 1][j][1],
        isShear: true,
      });
    }
  }

  // Draw edges
  const edgeSelection = svg.selectAll("line").data(edges);
  edgeSelection
    .enter()
    .append("line")
    .merge(edgeSelection)
    .attr("x1", (d) => d.x1)
    .attr("y1", (d) => d.y1)
    .attr("x2", (d) => d.x2)
    .attr("y2", (d) => d.y2)
    .style("stroke", (d) => (d.isShear ? "blue" : "lightgreen"))
    .style("stroke-width", (d) => (d.isShear ? 1 : 2));

  edgeSelection.exit().remove();
}

function updatePositions() {
  if (useAltForceMethod) {
    updatePositionsVerlet();
  } else {
    updatePositionsEuler();
  }
}

function calculateForces() {
  // Reset forces for all particles
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      forces[i][j][0] = 0; // Force in x-direction
      forces[i][j][1] = 0; // Force in y-direction
    }
  }

  // Structural Spring Constants
  const k = structuralSpringK; // Spring stiffness coefficient
  const b = structuralSpringB; // Damping coefficient
  const ℓ0 = structuralRestLength; // Rest length of the spring

  // Calculate forces for horizontal structural springs
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols - 1; j++) {
      // Positions of particles p and q
      const r_p = positions[i][j]; // [x_p, y_p]
      const r_q = positions[i][j + 1]; // [x_q, y_q]

      // Velocities of particles p and q
      const v_p = velocities[i][j]; // [v_px, v_py]
      const v_q = velocities[i][j + 1]; // [v_qx, v_qy]

      // Displacement vector components: Δr = r_p - r_q
      const Δr_x = r_p[0] - r_q[0]; // Δr_x = x_p - x_q
      const Δr_y = r_p[1] - r_q[1]; // Δr_y = y_p - y_q

      // Distance between particles p and q
      const distance = Math.sqrt(Δr_x * Δr_x + Δr_y * Δr_y);

      // Avoid division by zero
      if (distance === 0) continue;

      // Unit vector along the displacement: n̂ = Δr / |Δr|
      const n_x = Δr_x / distance;
      const n_y = Δr_y / distance;

      // **Spring Force Calculation**

      // Calculate the spring force magnitude: F_s = -k (|Δr| - ℓ0)
      const F_s_magnitude = -k * (distance - ℓ0);

      // Spring force components: F_s = F_s_magnitude * n̂
      const F_sx = F_s_magnitude * n_x;
      const F_sy = F_s_magnitude * n_y;

      // **Damping Force Calculation**

      // Relative velocity components: Δv = v_p - v_q
      const Δv_x = v_p[0] - v_q[0]; // Δv_x = v_px - v_qx
      const Δv_y = v_p[1] - v_q[1]; // Δv_y = v_py - v_qy

      // Damping force components: F_d = -b * Δv
      const F_dx = -b * Δv_x;
      const F_dy = -b * Δv_y;

      // **Total Force Components**

      // Total force on particle p due to particle q: F_total = F_s + F_d
      const F_total_x = F_sx + F_dx;
      const F_total_y = F_sy + F_dy;

      // **Apply Forces to Particles**

      // Update forces on particle p
      forces[i][j][0] += F_total_x;
      forces[i][j][1] += F_total_y;

      // Update forces on particle q (Newton's Third Law)
      forces[i][j + 1][0] -= F_total_x;
      forces[i][j + 1][1] -= F_total_y;
    }
  }
}

function updatePositionsEuler() {
  // TODO: think about how to calculate positions and velocities. (e.g. Euler's method)
  calculateForces();

  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      // TODO: potentially implement position and velocity updates here.
      // Example:
      // velocities[i][j][0] += some calculation
      // velocities[i][j][1] += some calculation
      // positions[i][j][0] += some calculation;
      // positions[i][j][1] += some calculation;
      // Calculate acceleration (F = ma -> a = F/m)

      const ax = forces[i][j][0] / nodeMass;
      const ay = forces[i][j][1] / nodeMass;

      // Update velocity: v = v + at
      velocities[i][j][0] += ax * timeStep;
      velocities[i][j][1] += ay * timeStep;

      // Update position: x = x + vt
      positions[i][j][0] += velocities[i][j][0] * timeStep;
      positions[i][j][1] += velocities[i][j][1] * timeStep;
    }
  }

  // TODO: Think about how to redraw nodes and edges with updated positions
  drawNodes();
  drawEdges();
}

function updatePositionsVerlet() {
  calculateForces();

  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      // Calculate acceleration (F = ma)
      const ax = forces[i][j][0] / nodeMass;
      const ay = forces[i][j][1] / nodeMass;

      // Store current position
      const oldX = positions[i][j][0];
      const oldY = positions[i][j][1];

      // Verlet integration
      // new_position = 2 * current_position - old_position + acceleration * (timestep)^2
      positions[i][j][0] =
        2 * positions[i][j][0] -
        oldPositions[i][j][0] +
        ax * timeStep * timeStep;
      positions[i][j][1] =
        2 * positions[i][j][1] -
        oldPositions[i][j][1] +
        ay * timeStep * timeStep;

      // Equation 12: v_{n+1} = (r_{n+1} - r_{n-1})/(2h)
      velocities[i][j][0] =
        (positions[i][j][0] - oldPositions[i][j][0]) / (2 * timeStep);
      velocities[i][j][1] =
        (positions[i][j][1] - oldPositions[i][j][1]) / (2 * timeStep);

      // Update old positions for next frame
      oldPositions[i][j][0] = oldX;
      oldPositions[i][j][1] = oldY;
    }
  }
  drawNodes();
  drawEdges();
}

/**
 * Main simulation loop.
 * Continuously updates the simulation as long as `isRunning` is true.
 */
function simulationLoop() {
  if (!isRunning) return;

  // TODO: think about how to implement the simulation loop. below are some functions that you might find useful.
  updatePositions();
  requestAnimationFrame(simulationLoop);
}

// ********** Event listeners examples for controls **********

document.getElementById("force-method").addEventListener("change", (e) => {
  useAltForceMethod = e.target.checked;
});

// Start/Stop simulation
document.getElementById("toggle-simulation").addEventListener("click", () => {
  isRunning = !isRunning;
  document.getElementById("toggle-simulation").innerText = isRunning
    ? "Stop Simulation"
    : "Start Simulation";
  if (isRunning) simulationLoop();
});

// Update grid rows
document.getElementById("rows").addEventListener("input", (e) => {
  rows = parseInt(e.target.value, 10);
  initializeGrid();
});

// Update grid columns
document.getElementById("cols").addEventListener("input", (e) => {
  cols = parseInt(e.target.value, 10);
  initializeGrid();
});

// Update restore force
document.getElementById("restore-force").addEventListener("input", (e) => {
  restoreForce = parseFloat(e.target.value);
  document.getElementById("restore-force-value").textContent =
    restoreForce.toFixed(2);
});

// Update damping
document.getElementById("damping").addEventListener("input", (e) => {
  damping = parseFloat(e.target.value);
  document.getElementById("damping-value").textContent = damping.toFixed(2);
});

// Initialize the simulation
initializeGrid();
// additional functions
