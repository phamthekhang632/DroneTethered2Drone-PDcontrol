clear
clc
close

% Parameters
N = 10;               % Number of masses
total_length = 100;                % Total length of the rope
mass = 0.05;           % Mass of each segment
k = 100;              % Spring stiffness
damping = 0.005;        % Damping coefficient
g = 9.81;             % Gravity (m/s^2)
dt = 0.01;            % Time step (s)
total_time = 30;                % Total simulation time (s)

% User-defined initial conditions
mother_ship_position = [total_length*0.75, total_length*1.2]; % [x, y] position of the mother ship
rope_angle = deg2rad(30);             % Angle of the rope with respect to the horizontal (in radians)

% Derived parameters
rest_length = total_length / (N-1);  % Rest length of each spring

% Initial conditions
positions = zeros(N, 2);   % [x, y] positions of the masses
velocities = zeros(N, 2);  % [vx, vy] velocities of the masses
for i = 1:N
    positions(i, 1) = mother_ship_position(1) + (i-1) * rest_length * cos(rope_angle);
    positions(i, 2) = mother_ship_position(2) - (i-1) * rest_length * sin(rope_angle);
end

% Fixed point (first mass)
fixedPoint = positions(1, :);

% Simulation loop
time = 0:dt:total_time;
trajectory = zeros(length(time), N, 2);

for t = 1:length(time)
    % Store current positions
    trajectory(t, :, :) = positions;
    
    % Compute forces on each mass
    forces = zeros(N, 2);
    
    for i = 2:N
        % Spring force from the left segment
        if i > 1
            vecLeft = positions(i, :) - positions(i-1, :);
            distLeft = norm(vecLeft);
            forceLeft = -k * (distLeft - rest_length) * (vecLeft / distLeft);
            forces(i, :) = forces(i, :) + forceLeft;
        end
        
        % Spring force from the right segment
        if i < N
            vecRight = positions(i+1, :) - positions(i, :);
            distRight = norm(vecRight);
            forceRight = k * (distRight - rest_length) * (vecRight / distRight);
            forces(i, :) = forces(i, :) + forceRight;
        end
        
        % Damping force
        forces(i, :) = forces(i, :) - damping * velocities(i, :);
        
        % Gravity force
        forces(i, 2) = forces(i, 2) - mass * g;
    end
    
    % Update positions and velocities using Euler integration
    for i = 2:N  % Skip the fixed mass
        accelerations = forces(i, :) / mass;
        velocities(i, :) = velocities(i, :) + accelerations * dt;
        positions(i, :) = positions(i, :) + velocities(i, :) * dt;
    end
    
    % Keep the first mass fixed
    positions(1, :) = fixedPoint;
    velocities(1, :) = [0, 0];
end

% Visualization
figure;
for t = 1:10:length(time)
    plot(squeeze(trajectory(t, :, 1)), squeeze(trajectory(t, :, 2)), '-o', 'LineWidth', 2);
    axis equal;
    xlim([0, total_length*1.5]);
    ylim([-10, mother_ship_position(2)*1.2]);
    title(sprintf('Time: %.2f s', time(t)));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    drawnow;
    pause(0.05);
end

% Extract end particle positions
end_particle_positions = squeeze(trajectory(:, end, :));

% Time vector
time = 0:dt:total_time;

% Plot end particle positions
figure;
plot(time, end_particle_positions(:, 1), 'b', 'DisplayName', 'X Position');
hold on;
plot(time, end_particle_positions(:, 2), 'r', 'DisplayName', 'Y Position');
title('End Particle Position vs Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;
hold off;
