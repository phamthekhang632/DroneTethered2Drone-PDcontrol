clear;
clc;
close all;

% Parameters
N = 10;               % Number of masses
total_length = 100;   % Total length of the rope
mass = 0.05;          % Mass of each segment
k = 100;              % Spring stiffness
damping = 0.05;      % Damping coefficient
g = 9.81;             % Gravity (m/s^2)
dt = 0.01;            % Time step (s)
total_time = 30;      % Total simulation time (s)

% User-defined initial conditions
mother_ship_position = [10, 110]; % [x, y] position of the mother ship
rope_angle = deg2rad(80);         % Angle of the rope with respect to the horizontal (in radians)

% Derived parameters
rest_length = total_length / (N-1);  % Rest length of each spring

% Initial conditions
positions = zeros(N, 2);   % [x, y] positions of the masses
velocities = zeros(N, 2);  % [vx, vy] velocities of the masses
for i = 1:N
    positions(i, 1) = mother_ship_position(1) + (i-1) * rest_length * cos(rope_angle);
    positions(i, 2) = mother_ship_position(2) - (i-1) * rest_length * sin(rope_angle);
end

% Simulation loop
time = 0:dt:total_time;
trajectory = zeros(length(time), N, 2);

% PD controller parameters
target_position = [mother_ship_position(1) + 75, mother_ship_position(2)]; % Target position for particle 1
Kp = 5;                    % Proportional gain
Kd = 5;                     % Derivative gain

% Parameters for velocity limit
max_velocity = 10; % Maximum allowable velocity for the mother ship (m/s)

% % Parameters for wind
wind_force = [-0.1, 0]; % Wind force vector (default: blowing from right to left, 2 N in magnitude)

for t = 1:length(time)
    % Store current positions
    trajectory(t, :, :) = positions;
    
    % Compute forces on each mass
    forces = zeros(N, 2);
    
    for i = 1:N
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
        
        % Add wind force to all particles
        % wind_force = [0.5*sin(t*dt), 0]; % Oscillating wind
        forces(i, :) = forces(i, :) + wind_force;
    end
    
    % Apply PD controller force to particle 1
    error = target_position - positions(1, :);            % Proportional term
    error_dot = -velocities(1, :);                        % Derivative term
    pd_force = Kp * error + Kd * error_dot;               % PD control force
    forces(1, :) = forces(1, :) + pd_force;               % Add PD control force to particle 1
    
    % Update positions and velocities using Euler integration
    for i = 1:N
        accelerations = forces(i, :) / mass;
        velocities(i, :) = velocities(i, :) + accelerations * dt;
        
        % Limit velocity of the mother ship (particle 1)
        if i == 1
            velocity_magnitude = norm(velocities(1, :));
            if velocity_magnitude > max_velocity
                velocities(1, :) = velocities(1, :) * (max_velocity / velocity_magnitude);
            end
        end
        
        positions(i, :) = positions(i, :) + velocities(i, :) * dt;
    end
end

% Calculate velocities from trajectory
mother_ship_velocity_traj = diff(squeeze(trajectory(:, 1, :))) / dt; % Velocity of the mother ship
minion_velocity_traj = diff(squeeze(trajectory(:, end, :))) / dt;    % Velocity of the minion

% Pad the velocities to match the time vector length (since diff reduces the size by 1)
mother_ship_velocity_traj = [mother_ship_velocity_traj; [0, 0]];
minion_velocity_traj = [minion_velocity_traj; [0, 0]];

% Visualization
figure;
for t = 1:10:length(time)
    % Plot the rope
    plot(squeeze(trajectory(t, :, 1)), squeeze(trajectory(t, :, 2)), '-o', 'LineWidth', 2);
    hold on;
    
    % Extract velocities for current frame
    mother_ship_velocity = mother_ship_velocity_traj(t, :);
    minion_velocity = minion_velocity_traj(t, :);
    
    % Print velocity vectors on the screen
    annotation_str = sprintf(['Time: %.2f s\n' ...
                              'Mother Ship Velocity: [%.2f, %.2f] m/s\n' ...
                              'Minion Velocity: [%.2f, %.2f] m/s\n'], ...
                              time(t), mother_ship_velocity(1), mother_ship_velocity(2), ...
                              minion_velocity(1), minion_velocity(2));
    text(5, mother_ship_position(2)+15, annotation_str, 'FontSize', 10, 'BackgroundColor', 'w');
    
    % Plot formatting
    axis equal;
    xlim([0, 150]);
    ylim([-10, mother_ship_position(2)+20]);
    title(sprintf('Time: %.2f s', time(t)));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    drawnow;
    pause(dt);
    hold off;
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
