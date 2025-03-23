import numpy as np
import matplotlib.pyplot as plt

# Parameters
num_points = 50
x = np.linspace(0, 2 * np.pi, num_points)
y = np.sin(x)

# Create a plot
plt.ion()  # Turn on interactive mode
fig, ax = plt.subplots()
line, = ax.plot(x, y, label="Real-time Sin Wave")  # Initial plot
ax.set_ylim(-2, 2)
ax.legend()

# Real-time update loop
for i in range(100):
    y = np.sin(x + i * 0.1)  # Update data
    line.set_ydata(y)  # Update the plot
    plt.pause(0.05)  # Pause for a brief moment (controls the refresh rate)

plt.ioff()  # Turn off interactive mode
plt.show()
