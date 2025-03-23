import numpy as np
import pygame
import matplotlib.pyplot as plt

# Particle class
class Particle:
    def __init__(self, position, mass=0.1):
        self.position = np.array(position, dtype=float)
        self.previous_position = np.array(position, dtype=float)
        self.velocity = np.array([0.0, 0.0], dtype=float)
        self.acceleration = np.array([0.0, 0.0], dtype=float)
        self.mass = mass

    def apply_force(self, force):
        force = np.array(force, dtype=float)
        self.acceleration += force / self.mass

# SpringDamper class
class SpringDamper:
    def __init__(self, particle1, particle2, rest_length, spring_constant, damping_constant):
        self.particle1 = particle1
        self.particle2 = particle2
        self.rest_length = rest_length
        self.spring_constant = spring_constant
        self.damping_constant = damping_constant

    def apply_force(self):
        # Compute the vector and distance between particles
        direction = self.particle2.position - self.particle1.position
        distance = np.linalg.norm(direction)
        if distance > 0:
            direction /= distance  # Normalize

            # Hooke's Law: F_spring = -k * (distance - rest_length)
            spring_force = -self.spring_constant * (distance - self.rest_length) * direction

            # Damping force: F_damp = -c * relative_velocity_along_direction
            relative_velocity = self.particle2.velocity - self.particle1.velocity
            damping_force = -self.damping_constant * np.dot(relative_velocity, direction) * direction

            # Total force
            force = spring_force + damping_force

            # Apply force to particles
            self.particle1.apply_force(force)
            self.particle2.apply_force(-force)

# Simulation parameters
dt = 0.01  # Time step
num_steps = 500  # Number of steps
num_particles = 3  # Number of particles
rope_length = 100.0  # Length of the rope
rest_length = rope_length / (num_particles - 1)
gravity = np.array([0.0, -9.81])  # Gravity vector
spring_constant = 100.0  # Spring stiffness
damping_constant = 0.5  # Damping coefficient

# Mothership position (user-defined)
mothership_position = np.array([0.0, 120])  # Set the mothership's position

# Initialize particles with positions relaxed below the mothership
particles = [Particle(position=mothership_position - np.array([i * rest_length, 0.0])) for i in range(num_particles)]

# Change mass of the drones
particles[0].mass = 10  # Mothership mass
particles[-1].mass = 1  # End particle mass

# Add spring-damper connections
springs = []
for i in range(num_particles - 1):
    springs.append(SpringDamper(particles[i], particles[i + 1], rest_length, spring_constant, damping_constant))

# Initialize Pygame for visualization
pygame.init()
screen = pygame.display.set_mode((800, 600))
pygame.display.set_caption("Rope Simulation with Spring-Damper Model")
clock = pygame.time.Clock()

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

    # Apply gravity and reset accelerations
    for p in particles:
        p.acceleration = gravity

    # Apply hover thrust to the mothership
    hover_thrust = np.array([0.0, -particles[0].mass * gravity[1]])  # Thrust to counteract gravity
    particles[0].apply_force(hover_thrust)

    # Apply spring-damper forces
    for spring in springs:
        spring.apply_force()

    # Update particle positions
    for p in particles:
        p.velocity += p.acceleration * dt
        p.position += p.velocity * dt

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
