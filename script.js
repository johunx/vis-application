const mass = 0.2; // kg
const springLength = 1; // meters
const stiffness = 20; // kg/s^2
const damping = 0.1; // kg/s
const timeStep = 0.01; // s

let positions = [0, springLength];
let velocities = [0, 0];
let accelerations = [0, 0];

function calculateSpringForce(position1, position2) {
  const displacement = position2 - position1;
  return -stiffness * (displacement - springLength);
}

function calculateDampingForce(velocity1, velocity2) {
  return -damping * (velocity1 - velocity2);
}

function updateSystem() {
  // Calculate forces
  let springForce = calculateSpringForce(positions[0], positions[1]);
  let dampingForce = calculateDampingForce(velocities[0], velocities[1]);
  let totalForce = springForce + dampingForce;

  // Calculate acceleration (Newton's Second Law)
  accelerations[0] = totalForce / mass;
  accelerations[1] = -totalForce / mass;

  // Update velocity and position using Euler method
  for (let i = 0; i < positions.length; i++) {
    velocities[i] += accelerations[i] * timeStep;
    positions[i] += velocities[i] * timeStep;
  }
}

// Initialize SVG before using it in the render function
const svg = d3
  .select("#visualization")
  .append("svg")
  .attr("width", 600)
  .attr("height", 200);

// Update loop
function simulate() {
  updateSystem();
  render(); // Update visualization
  requestAnimationFrame(simulate);
}

simulate(); // Start the simulation

function render() {
  svg.selectAll("*").remove();

  // Draw masses
  svg
    .append("circle")
    .attr("cx", 100 + positions[0] * 50)
    .attr("cy", 100)
    .attr("r", 10)
    .attr("fill", "blue");

  svg
    .append("circle")
    .attr("cx", 100 + positions[1] * 50)
    .attr("cy", 100)
    .attr("r", 10)
    .attr("fill", "blue");

  // Draw spring
  svg
    .append("line")
    .attr("x1", 100 + positions[0] * 50)
    .attr("y1", 100)
    .attr("x2", 100 + positions[1] * 50)
    .attr("y2", 100)
    .attr("stroke", "black");
}