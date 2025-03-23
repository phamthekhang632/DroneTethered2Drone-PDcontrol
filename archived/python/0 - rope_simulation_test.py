import numpy as np
import pygame
import matplotlib.pyplot as plt

# Particle class
class Particle:
    def __init__(self, position, mass=1.0):
        self.position = np.array(position, dtype=float)
        self.previous_position = np.array(position, dtype=float)
        self.acceleration = np.array([0.0, 0.0], dtype=float)
        self.mass = mass

    def get_acceleration(self):
        return self.acceleration

# Constraint class
class Constraint:
    def __init__(self, particle1, particle2, desired_distance):
        self.particle1 = particle1
        self.particle2 = particle2
        self.desired_distance = desired_distance

    def relax(self):
        direction = self.particle2.position - self.particle1.position
        distance = np.linalg.norm(direction)
        if distance > 0:
            direction /= distance
            delta_distance = distance - self.desired_distance
            correction = delta_distance * direction / 2.0
            self.particle1.position += correction
            self.particle2.position -= correction

    def relax_fixed_start(self):
        direction = self.particle2.position - self.particle1.position
        distance = np.linalg.norm(direction)
        if distance > 0:
            direction /= distance
            delta_distance = distance - self.desired_distance
            correction = delta_distance * direction / 2.0
            # self.particle1.position += correction
            self.particle2.position -= 2 * correction

    def relax_fixed_end(self):
        direction = self.particle2.position - self.particle1.position
        distance = np.linalg.norm(direction)
        if distance > 0:
            direction /= distance
            delta_distance = distance - self.desired_distance
            correction = delta_distance * direction / 2.0
            self.particle1.position += 2 * correction
            # self.particle2.position -= 2 *correction

# Verlet integration function
def verlet_integration(particles, dt, fixed_particle_index):
    for i, p in enumerate(particles):
        if i == fixed_particle_index:
            continue  # Skip updating the fixed particle
        temp_position = p.position.copy()
        p.position = 2 * p.position - p.previous_position + dt**2 * p.get_acceleration()
        p.previous_position = temp_position

# Jakobsen method to enforce constraints
def jakobsen(constraints, iterations, fixed_particle_index):
    for _ in range(iterations):
        for i, constraint in enumerate(constraints):
            if i == fixed_particle_index:
                constraint.relax_fixed_start()
                continue
            constraint.relax()

# Simulation parameters
dt = 0.01  # Time step
num_steps = 500  # Number of steps
num_particles = 10  # Number of particles
rope_length = 10.0  # Length of the rope
desired_distance = rope_length / (num_particles - 1)
gravity = np.array([0.0, -9.81])  # Gravity vector
iterations = 5  # Jakobsen relaxation iterations per time step

# Initialize particles
particles = [Particle(position=[i * desired_distance, 0.0]) for i in range(num_particles)]

# Fix the first particle
fixed_particle_index = 0
particles[fixed_particle_index].previous_position = particles[0].position  # Keep fixed
particles[fixed_particle_index].get_acceleration = lambda: np.array([0.0, 0.0])  # No forces

# Add constraints
constraints = []
for i in range(num_particles - 1):
    constraints.append(Constraint(particles[i], particles[i + 1], desired_distance))

# Add gravity force
for i, p in enumerate(particles):
    if i == fixed_particle_index:
        continue  # Skip updating the fixed particle
    p.get_acceleration = lambda p=p: gravity

# Initialize Pygame for visualization
pygame.init()
screen = pygame.display.set_mode((800, 600))
pygame.display.set_caption("Rope Simulation with Verlet Integration and Jakobsen Constraints")
clock = pygame.time.Clock()

# Pygame Scaling
scale = 50  # Scale factor to fit rope in screen space
offset = np.array([400, 100])  # Screen offset to center the rope

# End particle position tracking
end_positions = []

# Simulation loop
running = True
for step in range(num_steps):
    # Check for Pygame events to allow quitting
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
            break
    if not running:
        break

    # Verlet integration (ensure fixed particle stays fixed)
    verlet_integration(particles, dt, fixed_particle_index=0)

    # Enforce constraints using Jakobsen method
    jakobsen(constraints, iterations, fixed_particle_index=0)

    # Track end particle position
    end_positions.append(particles[-1].position.copy())

    # Pygame Visualization
    screen.fill((255, 255, 255))  # Clear screen with white background
    for i in range(num_particles - 1):
        # Draw rope as lines between particles
        p1_screen = particles[i].position * [scale, -scale] + offset
        p2_screen = particles[i + 1].position * [scale, -scale] + offset
        pygame.draw.line(screen, (0, 0, 0), p1_screen, p2_screen, 2)

    # Draw particles
    for p in particles:
        p_screen = p.position * [scale, -scale] + offset
        pygame.draw.circle(screen, (0, 0, 255), p_screen.astype(int), 5)

    # Highlight the fixed particle and end particle
    fixed_particle_screen = particles[0].position * [scale, -scale] + offset
    end_particle_screen = particles[-1].position * [scale, -scale] + offset
    pygame.draw.circle(screen, (0, 255, 0), fixed_particle_screen.astype(int), 8)  # Fixed
    pygame.draw.circle(screen, (255, 0, 0), end_particle_screen.astype(int), 8)  # End

    pygame.display.flip()
    clock.tick(60)  # Limit to 60 FPS

pygame.quit()

# Plot the trajectory of the end particle
end_positions = np.array(end_positions)
plt.figure(figsize=(8, 6))
plt.plot(end_positions[:, 0], end_positions[:, 1], label="End Particle Trajectory")
plt.scatter(end_positions[-1, 0], end_positions[-1, 1], color='red', label="Final Position")
plt.title("End Particle Trajectory")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
plt.grid()
plt.legend()
plt.show()
