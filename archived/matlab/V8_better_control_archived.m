clear;
clc;
close all;
%%
% Parameters
N = 10;               % Number of masses
total_length = 100;   % Total length of the rope
default_mass = 0.05;  % Default mass of each segment
masses = [6.47, repmat(default_mass, 1, N-2), 1]; % Masses
k = 100;              % Spring stiffness
damping = 0.01;       % Damping coefficient
g = 9.81;             % Gravity (m/s^2)
dt = 0.01;            % Time step (s)
total_time = 30;      % Total simulation time (s)

% User-defined initial conditions
mother_ship_position = [10, 110]; % [x, y] position of the mother ship
rope_angle = deg2rad(90);         % Angle of the rope with respect to the horizontal (in radians)

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

% Mothership control
mothership_target_position = [positions(1, 1) + 75, positions(1, 2)]; % Target for mothership
mothership_Kp = 10; mothership_Kd = 5; % PD gains for x control
mothership_max_velocity = 12; % Velocity limit

% Initial states
minion_l = 0.2;
minion_I = 0.5*masses(N)*minion_l^2;        % Moment of inertia (kg·m²)
minion_max_thrust = 20;      % Maximum thrust (N)
minion_theta = zeros(length(time), 1);  % Initial pitch angle (radians)
minion_theta(1) = deg2rad(-30);
minion_omega = zeros(length(time));
minion_Kp = diag([5, 20]);
minion_Kd = diag([10, 10]);
minion_thrust = zeros(length(0:dt:total_time), 2);
minion_error_record = zeros(length(time), 2);
minion_theta_desired_record = zeros(length(time), 1);

% Parameters for wind variation
wind_base = [-1, 0];              % Base wind force
wind_amplitude = [0.1, 0.05];       % Amplitude of sinusoidal variation
wind_frequency = [0.5, 0.3];        % Frequency of sinusoidal variation (Hz)
wind_noise_amplitude = [0.05, 0.02];% Amplitude of random noise
wind_record = zeros(length(time), 2);

for t = 1:length(time)
    % Store current positions
    trajectory(t, :, :) = positions;

    % Forces on masses
    forces = zeros(N, 2);
    for i = 1:N
        % Spring forces
        if i > 1
            vecLeft = positions(i, :) - positions(i-1, :);
            distLeft = norm(vecLeft);
            forceLeft = -k * (distLeft - rest_length) * (vecLeft / distLeft);
            forces(i, :) = forces(i, :) + forceLeft;
        end
        if i < N
            vecRight = positions(i+1, :) - positions(i, :);
            distRight = norm(vecRight);
            forceRight = k * (distRight - rest_length) * (vecRight / distRight);
            forces(i, :) = forces(i, :) + forceRight;
        end

        % Damping and gravity
        forces(i, :) = forces(i, :) - damping * velocities(i, :);
        forces(i, 2) = forces(i, 2) - masses(i) * g;
    end

    % Mothership PD control
    mothership_error = mothership_target_position - positions(1, :);
    mothership_error_dot = -velocities(1, :);
    pd_force = mothership_Kp * mothership_error + mothership_Kd * mothership_error_dot + [0, masses(1) * g];
    forces(1, :) = forces(1, :) + pd_force;

    % Minion PD control
    % Error (x, theta)
    minion_state_current = [positions(N, 1); minion_theta(t)];
    minion_state_desired = [positions(1, 1); 0];
    minion_error = minion_state_desired - minion_state_current;

    % Calculate average of the last (2 / dt) frames of x-error
    frames_to_average = round(0.5 / dt); % Number of frames in the last 2 seconds
    if t > frames_to_average
        avg_x_error = mean(minion_error_record(t-frames_to_average:t-1, 1));
    else
        avg_x_error = mean(minion_error_record(1:t-1, 1)); % For initial frames
    end
    
    % Threshold for x-error to recalculate angular error
    x_error_threshold = 2; % Example threshold value
    gain = 0.1;    % Scaling factor for angular adjustment
    
    if abs(avg_x_error) > x_error_threshold
        % Scale and clip the desired angular state
        desired_theta = gain * ( - avg_x_error);
        desired_theta = min(max(desired_theta, -pi/4), pi/4); % Clip to [-pi/4, pi/4]
    
        % Recalculate angular error (theta) based on scaled and clipped value
        minion_state_desired(2) = desired_theta;
        minion_error(2) = minion_state_desired(2) - minion_theta(t);
    end

    % Update error record
    minion_theta_desired_record(t) = minion_state_desired(2);
    minion_error_record(t, :) = minion_error;

    minion_state_dot_current = [velocities(N, 1); minion_omega(t)];
    minion_error_dot = -minion_state_dot_current;

    % PD control law
    control_input = (minion_Kp * minion_error) + minion_Kd * minion_error_dot;
    thrust = masses(N) * (control_input(1)+g) / cos(minion_theta(t));
    torque = minion_I * control_input(2);

    F2 = (thrust / 2) + (torque / minion_l);
    F1 = (thrust / 2) - (torque / minion_l);

    % Clamp forces within limits
    F1 = min(max(F1, 0), minion_max_thrust);
    F2 = min(max(F2, 0), minion_max_thrust);

    thrust_force = (F1+F2) * [cos(minion_theta(t)), sin(minion_theta(t))];
    forces(N, :) = forces(N, :) + thrust_force;
    minion_thrust(t, :) = [F1, F2];

    % Angular acceleration
    minion_alpha = torque / minion_I;    % Angular acceleration

    % Update angular velocity
    if t==1
        minion_omega(t) = minion_alpha * dt;
    else 
        minion_omega(t) = minion_omega(t-1) + minion_alpha * dt;
    end

    % Update angle
    if t==1
        minion_theta(t) = minion_omega(t) * dt;
    else
        minion_theta(t) = minion_theta(t-1) + minion_omega(t) * dt;
    end

    % Update positions and velocities (Euler integration)
    for i = 1:N
        acceleration = forces(i, :) / masses(i);
        velocities(i, :) = velocities(i, :) + acceleration * dt;
        if i == 1 && norm(velocities(1, :)) > mothership_max_velocity
            velocities(1, :) = velocities(1, :) * (mothership_max_velocity / norm(velocities(1, :)));
        end
        positions(i, :) = positions(i, :) + velocities(i, :) * dt;
    end
end
%%
% Calculate velocities from trajectory
mother_ship_velocity_traj = diff(squeeze(trajectory(:, 1, :))) / dt; % Velocity of the mother ship
minion_velocity_traj = diff(squeeze(trajectory(:, end, :))) / dt;    % Velocity of the minion

% Pad the velocities to match the time vector length (since diff reduces the size by 1)
mother_ship_velocity_traj = [mother_ship_velocity_traj; [0, 0]];
minion_velocity_traj = [minion_velocity_traj; [0, 0]];

% % Visualization
% figure;
% for t = 1:10:length(time)
%     % Plot the rope
%     plot(squeeze(trajectory(t, :, 1)), squeeze(trajectory(t, :, 2)), '-o', 'LineWidth', 2);
%     hold on;
%     plot(squeeze(trajectory(t, 1, 1)), squeeze(trajectory(t, 1, 2)), 'ko', 'LineWidth', masses(1), 'DisplayName', 'Mothership');
%     plot(squeeze(trajectory(t, N, 1)), squeeze(trajectory(t, N, 2)), 'yo', 'LineWidth', masses(1), 'DisplayName', 'Minion');
% 
%     % Plot the target point for the mothership
%     % plot(mothership_target_position(1), mothership_target_position(2)+drop_down_distance, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Initial Target');
%     plot(mothership_target_position(1), mothership_target_position(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Initial Target');
%     % plot(mothership_target_position(1), mothership_target_position(2), 'bx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Final Target');
% 
%     % Extract velocities for current frame
%     mother_ship_velocity = mother_ship_velocity_traj(t, :);
%     minion_velocity = minion_velocity_traj(t, :);
%     wind_force = wind_record(t, :);
% 
%     % Print velocity vectors and wind force on the screen
%     annotation_str = sprintf(['Time: %.2f s\n' ...
%                               'Mother Ship Velocity: [%.2f, %.2f] m/s\n' ...
%                               'Minion Velocity: [%.2f, %.2f] m/s\n' ...
%                               'Wind Force: [%.2f, %.2f] N\n'], ...
%                               time(t), mother_ship_velocity(1), mother_ship_velocity(2), ...
%                               minion_velocity(1), minion_velocity(2), ...
%                               wind_force(1), wind_force(2));
%     text(5, mother_ship_position(2)+15, annotation_str, 'FontSize', 10, 'BackgroundColor', 'w');
% 
%     % Plot formatting
%     axis equal;
%     xlim([0, 150]);
%     ylim([-10, mother_ship_position(2)+20]);
%     title(sprintf('Time: %.2f s', time(t)));
%     xlabel('X Position (m)');
%     ylabel('Y Position (m)');
%     drawnow;
%     pause(dt);
%     hold off;
% end
%%
close all;

% End particle position plot
minion_positions = squeeze(trajectory(:, end, :));
mothership_position = squeeze(trajectory(:, 1, :));
time = 0:dt:total_time;
figure;
plot(time, minion_positions(:, 1), 'y', 'DisplayName', 'Minion');
hold on;
plot(time, mothership_position(:, 1), 'k', 'DisplayName', 'Mothership');
title('End Particle Position vs Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;
hold off;

% Plot Minion Thrust Over Time
figure;
plot(time, minion_thrust(:, 1), 'b', 'LineWidth', 2, 'DisplayName', 'F1');
hold on;
plot(time, minion_thrust(:, 2), '--r', 'LineWidth', 2, 'DisplayName', 'F2');
title('Minion Thrust Over Time');
xlabel('Time (s)');
ylabel('Thrust (N)');
grid on;
legend();
hold off;

% Create a figure with two subplots
figure;

% Subplot 1: Position error (X)
subplot(2, 1, 1);
plot(time, minion_error_record(:, 1), 'b', 'LineWidth', 2, 'DisplayName', 'X Position Error');
hold on;
plot(time, 5*ones(size(time)), '--r', 'LineWidth', 2, 'DisplayName', 'Limit');
plot(time, -5*ones(size(time)), '--r', 'LineWidth', 2, 'DisplayName', 'Limit');
title('Minion Position Error Over Time');
xlabel('Time (s)');
ylabel('X Error');
grid on;
legend();
hold off;

% Subplot 2: Pitch error
subplot(2, 1, 2);
plot(time, minion_error_record(:, 2)/pi*180, 'b', 'LineWidth', 2, 'DisplayName', 'Pitch Error');
hold on;
plot(time, minion_theta_desired_record/pi*180, '--r', 'LineWidth', 2, 'DisplayName', 'Pitch Desired');
% plot(time, minion_theta/pi*180, 'r', 'LineWidth', 2, 'DisplayName', 'Pitch');
title('Minion Pitch Error Over Time');
xlabel('Time (s)');
ylabel('Pitch Error (degrees)');
grid on;
legend();
hold off;