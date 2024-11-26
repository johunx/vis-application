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
let velocities = [];
let forces = [];
let isRunning = false;

/**
 * Initialize the grid with nodes and reset their positions, velocities, and forces.
 */
function initializeGrid() {
  positions = [];
  velocities = [];
  forces = [];
  const xStep = (width - 2 * padding) / (cols - 1);
  const yStep = (height - 2 * padding) / (rows - 1);

  for (let i = 0; i < rows; i++) {
    const positionRow = [];
    const velocityRow = [];
    const forceRow = [];
    for (let j = 0; j < cols; j++) {
      positionRow.push([padding + j * xStep, padding + yStep * i]); // ! TODO: think about how to calculate initial positions for the nodes
      velocityRow.push([0, 0]); // Initial velocity
      forceRow.push([0, 0]); // Initial force
    }
    positions.push(positionRow);
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

/**
 * Calculate forces acting on each node.
 * This function is a placeholder for students to implement force calculations.
 */
function calculateForces() {
  // Reset forces
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      forces[i][j][0] = 0;
      forces[i][j][1] = 0;
    }
  }

  // TODO: add your implementation here.
  // Example:
  // - Calculate spring forces (horizontal, vertical, diagonal/sheer).
  // - Add restoring forces.
  // - Add damping forces.
  // Let's start with just structural springs (horizontal connections)
  // Calculate forces for horizontal springs
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols - 1; j++) {
      // Get positions and velocities
      const pos1 = positions[i][j];
      const pos2 = positions[i][j + 1];
      const vel1 = velocities[i][j];
      const vel2 = velocities[i][j + 1];

      // Calculate spring direction and length
      const dx = pos2[0] - pos1[0];
      const dy = pos2[1] - pos1[1];
      const length = Math.sqrt(dx * dx + dy * dy);

      if (length === 0) continue;

      // Calculate relative velocity
      const dvx = vel2[0] - vel1[0];
      const dvy = vel2[1] - vel1[1];

      // Spring force (as before)
      const springForce = structuralSpringK * (length - structuralRestLength);

      // Add damping force
      const dampingForce = structuralSpringB * ((dvx * dx + dvy * dy) / length);

      // Total force
      const totalForce = springForce + dampingForce;
      const forceMagnitude = totalForce / length;
      const forceX = forceMagnitude * dx;
      const forceY = forceMagnitude * dy;

      // Apply forces
      forces[i][j][0] += forceX;
      forces[i][j][1] += forceY;
      forces[i][j + 1][0] -= forceX;
      forces[i][j + 1][1] -= forceY;
    }
  }
}

function updatePositions() {
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
