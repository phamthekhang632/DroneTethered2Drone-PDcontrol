clear;
clc;
close all;

% Parameters
N = 10;               % Number of masses
total_length = 200;   % Total length of the rope
default_mass = 0.05;  % Default mass of each segment
masses = [6.47, repmat(default_mass, 1, N-2), 1]; % Masses: Mothership (10kg), Minion (1kg), Others (0.05kg)
k = 100;              % Spring stiffness
damping = 0.05;       % Damping coefficient
g = 9.81;             % Gravity (m/s^2)
dt = 0.01;            % Time step (s)
total_time = 30;      % Total simulation time (s)

% User-defined initial conditions
mother_ship_position = [10, 240]; % [x, y] position of the mother ship
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

% Mothership controller
target_position = [mother_ship_position(1) + 150, mother_ship_position(2)]; % Target position for particle 1
Kp = 10;                    % Proportional gain
Kd = 5;                    % Derivative gain
max_velocity_mothership = 12;         % Maximum allowable velocity for the mother ship (m/s)

% PD controller parameters for the minion
Kp_minion = 8;      % Proportional gain for minion
Kd_minion = 4;      % Derivative gain for minion
max_force_minion = 20;  % Maximum horizontal thrust (N)
minion_thrust = zeros(length(time), 1);

% Parameters for wind variation
wind_base = [-1, 0];              % Base wind force
wind_amplitude = [0.1, 0.05];       % Amplitude of sinusoidal variation
wind_frequency = [0.5, 0.3];        % Frequency of sinusoidal variation (Hz)
wind_noise_amplitude = [0.05, 0.02];% Amplitude of random noise
wind_record = zeros(length(time), 2);

dropped_package = false;
drop_down_distance = 20;

% Simulation loop
for t = 1:length(time)
    % Dynamic wind force (sinusoidal + random noise)
    wind_force(1) = wind_base(1) + wind_amplitude(1) * sin(2 * pi * wind_frequency(1) * time(t)) ...
                    + wind_noise_amplitude(1) * (2 * rand() - 1);
    wind_force(2) = wind_base(2) + wind_amplitude(2) * sin(2 * pi * wind_frequency(2) * time(t)) ...
                    + wind_noise_amplitude(2) * (2 * rand() - 1);
    wind_record(t, :) = wind_force;

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
        forces(i, 2) = forces(i, 2) - masses(i) * g;
        
        % Add wind force to all particles
        forces(i, :) = forces(i, :) + wind_force;
    end
    
    % Apply PD controller force to the mothership (particle 1)
    error = target_position - positions(1, :);            % Proportional term
    error_dot = -velocities(1, :);                        % Derivative term
    pd_force = Kp * error + Kd * error_dot + [0, masses(1) * g];
    forces(1, :) = forces(1, :) + pd_force;
    
    % Apply PD controller for the minion (last particle)
    target_position_minion = positions(1, 1);           % Target x is the mothership's current x
    error_minion_x = target_position_minion - positions(N, 1); % Proportional term
    error_dot_minion_x = -velocities(N, 1);               % Derivative term
    pd_force_minion_x = Kp_minion * error_minion_x + Kd_minion * error_dot_minion_x;
    pd_force_minion_x = min(max(pd_force_minion_x, -max_force_minion), max_force_minion);
    forces(N, 1) = forces(N, 1) + pd_force_minion_x;  % Apply the constrained horizontal thrust
    minion_thrust(t) = pd_force_minion_x;

    % Update positions and velocities using Euler integration
    for i = 1:N
        accelerations = forces(i, :) / masses(i);
        velocities(i, :) = velocities(i, :) + accelerations * dt;
        
        % Limit velocity of the mothership (particle 1)
        if i == 1
            velocity_magnitude = norm(velocities(1, :));
            if velocity_magnitude > max_velocity_mothership
                velocities(1, :) = velocities(1, :) * (max_velocity_mothership / velocity_magnitude);
            end
        end
        
        % Update positions
        positions(i, :) = positions(i, :) + velocities(i, :) * dt;
    end

    % Check if the mothership has reached the target
    if norm(positions(1, :) - target_position) < 1.0 && dropped_package == false
        % Move the target downward
        target_position(2) = target_position(2) - drop_down_distance;
        fprintf('New target position: [%.2f, %.2f]\n', target_position(1), target_position(2));
        dropped_package = true;
    end
end
%%
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
    plot(squeeze(trajectory(t, 1, 1)), squeeze(trajectory(t, 1, 2)), 'ko', 'LineWidth', masses(1), 'DisplayName', 'Mothership');
    plot(squeeze(trajectory(t, N, 1)), squeeze(trajectory(t, N, 2)), 'yo', 'LineWidth', masses(1), 'DisplayName', 'Minion');

    % Plot the target point for the mothership
    plot(target_position(1), target_position(2)+drop_down_distance, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Initial Target');
    plot(target_position(1), target_position(2), 'bx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Final Target');

    % Extract velocities for current frame
    mother_ship_velocity = mother_ship_velocity_traj(t, :);
    minion_velocity = minion_velocity_traj(t, :);
    wind_force = wind_record(t, :);
    
    % Print velocity vectors and wind force on the screen
    annotation_str = sprintf(['Time: %.2f s\n' ...
                              'Mother Ship Velocity: [%.2f, %.2f] m/s\n' ...
                              'Minion Velocity: [%.2f, %.2f] m/s\n' ...
                              'Wind Force: [%.2f, %.2f] N\n'], ...
                              time(t), mother_ship_velocity(1), mother_ship_velocity(2), ...
                              minion_velocity(1), minion_velocity(2), ...
                              wind_force(1), wind_force(2));
    text(5, mother_ship_position(2)+15, annotation_str, 'FontSize', 10, 'BackgroundColor', 'w');
    
    % Plot formatting
    axis equal;
    xlim([0, target_position(1)*1.5]);
    ylim([-10, mother_ship_position(2)+20]);
    title(sprintf('Time: %.2f s', time(t)));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    drawnow;
    pause(dt);
    hold off;
end


% End particle position plot
end_particle_positions = squeeze(trajectory(:, end, :));
time = 0:dt:total_time;
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
