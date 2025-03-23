import numpy as np
import pygame
import matplotlib.pyplot as plt

# Particle class
class Particle:
    def __init__(self, position, mass=0.1):
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
            correction = delta_distance * direction
            self.particle2.position -= 2 * correction

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
dt = 0.1  # Time step
num_steps = 500  # Number of steps
num_particles = 3  # Number of particles
rope_length = 100.0  # Length of the rope
desired_distance = rope_length / (num_particles - 1)
gravity = np.array([0.0, -9.81])  # Gravity vector
iterations = 5  # Jakobsen relaxation iterations per time step

# Mothership position (user-defined)
mothership_position = np.array([0.0, 120])  # Set the mothership's position

# Initialize particles with positions relaxed below the mothership
particles = [Particle(position=mothership_position - np.array([0.0, i * desired_distance])) for i in range(num_particles)]

# Change mass of the minion
particles[-1].mass = 1

# Fix the first particle (mothership)
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

pygame.font.init()
font = pygame.font.Font(None, 24)  # Use the default font with size 24

# Pygame Scaling
scale = 3  # Scale factor to fit rope in screen space
offset = np.array([200, 500])  # Screen offset to center the rope

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

    # Move the mothership to the right 1m/s
    particles[0].position += np.array([5.0 * dt, 0.0])

    # Verlet integration (ensure fixed particle stays fixed)
    verlet_integration(particles, dt, fixed_particle_index=0)

    # Enforce constraints using Jakobsen method
    jakobsen(constraints, iterations, fixed_particle_index=0)

    # Track end particle position
    end_positions.append(particles[-1].position.copy())

    # Pygame Visualization
    screen.fill((255, 255, 255))  # Clear screen with white background

    # Draw ground
    pygame.draw.line(screen, (0, 0, 0), [0, offset[1]], [800, offset[1]], 3)  # Ground line at y=0

    for i in range(num_particles - 1):
        # Draw rope as lines between particles
        p1_screen = particles[i].position * [scale, -scale] + offset
        p2_screen = particles[i + 1].position * [scale, -scale] + offset
        pygame.draw.line(screen, (0, 0, 0), p1_screen, p2_screen, 2)

    # Draw particles
    for i, p in enumerate(particles):
        p_screen = p.position * [scale, -scale] + offset
        color = (0, 0, 0) if i == 0 else (0, 0, 255)  # Mothership as black, others as blue
        pygame.draw.circle(screen, color, p_screen.astype(int), 5)

    # Draw mothership (black rectangle)
    mothership_screen = particles[0].position * [scale, -scale] + offset
    pygame.draw.rect(screen, (0, 0, 0), (mothership_screen[0] - 10, mothership_screen[1] - 5, 20, 10))

    # Highlight end particle
    end_particle_screen = particles[-1].position * [scale, -scale] + offset
    pygame.draw.circle(screen, (255, 0, 0), end_particle_screen.astype(int), 8)  # End

    fps = clock.get_fps()  # Get the FPS
    fps_text = font.render(f"FPS: {fps:.2f}", True, (0, 0, 0))  # Render FPS text
    screen.blit(fps_text, (10, 10))  # Draw FPS text at the top-left corner

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
