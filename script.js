// Constants from the assignment
const h = 0.01; // Step size (s)
const m = 0.2; // Mass (kg)
const k = 20; // Spring stiffness (kg/s^2)
const b = 0.1; // Damping coefficient (kg/s)
const l0 = 1; // Rest length (m)

// Simulation state
let isRunning = false;
let masses = [];
let previousPositions = [];

// SVG setup
const width = 800;
const height = 400;
const svg = d3
  .select("#simulation")
  .append("svg")
  .attr("width", width)
  .attr("height", height);

// Initialize masses
function initializeMasses() {
  masses = [
    {
      x: width / 2 - 50,
      y: height / 2,
      vx: 0,
      vy: 0,
      fx: 0,
      fy: 0,
      fixed: true,
    },
    {
      x: width / 2 + 50,
      y: height / 2,
      vx: 0,
      vy: 0,
      fx: 0,
      fy: 0,
      fixed: false,
    },
  ];
  previousPositions = masses.map((m) => ({ x: m.x, y: m.y }));
}

// Calculate spring force between two masses
function calculateSpringForce(mass1, mass2) {
  const dx = mass2.x - mass1.x;
  const dy = mass2.y - mass1.y;
  const distance = Math.sqrt(dx * dx + dy * dy);
  const force = k * (distance - l0);

  return {
    fx: (force * dx) / distance,
    fy: (force * dy) / distance,
  };
}

// Calculate damping force between two masses
function calculateDampingForce(mass1, mass2) {
  const dvx = mass2.vx - mass1.vx;
  const dvy = mass2.vy - mass1.vy;

  return {
    fx: b * dvx,
    fy: b * dvy,
  };
}

// Euler integration step
function eulerStep() {
  // Calculate forces
  const springForce = calculateSpringForce(masses[0], masses[1]);
  const dampingForce = calculateDampingForce(masses[0], masses[1]);

  // Apply forces to non-fixed mass
  if (!masses[1].fixed) {
    masses[1].fx = -springForce.fx - dampingForce.fx;
    masses[1].fy = -springForce.fy - dampingForce.fy;

    // Update velocity and position using Euler method
    masses[1].vx += (masses[1].fx / m) * h;
    masses[1].vy += (masses[1].fy / m) * h;
    masses[1].x += masses[1].vx * h;
    masses[1].y += masses[1].vy * h;
  }
}

// Draw the simulation
function draw() {
  // Draw spring
  const spring = svg.selectAll(".spring").data([masses]);

  spring
    .enter()
    .append("line")
    .attr("class", "spring")
    .merge(spring)
    .attr("x1", (d) => d[0].x)
    .attr("y1", (d) => d[0].y)
    .attr("x2", (d) => d[1].x)
    .attr("y2", (d) => d[1].y);

  // Draw masses
  const circles = svg.selectAll(".mass").data(masses);

  circles
    .enter()
    .append("circle")
    .attr("class", "mass")
    .attr("r", 10)
    .merge(circles)
    .attr("cx", (d) => d.x)
    .attr("cy", (d) => d.y);
}

// Animation loop
function animate() {
  if (isRunning) {
    eulerStep();
    draw();
    requestAnimationFrame(animate);
  }
}

// Event listeners
d3.select("#startBtn").on("click", () => {
  isRunning = !isRunning;
  if (isRunning) animate();
});

d3.select("#resetBtn").on("click", () => {
  isRunning = false;
  initializeMasses();
  draw();
});

// Initialize and draw initial state
initializeMasses();
draw();
