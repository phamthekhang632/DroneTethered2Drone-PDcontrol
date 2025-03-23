clear;
clc;
close all;

% Parameters
N = 10;               % Number of masses
total_length = 100;   % Total length of the rope
default_mass = 0.05;  % Default mass of each segment(
masses = [6.47, repmat(default_mass, 1, N-2), 1]; % Masses: Mothership (10kg), Minion (1kg), Others (0.05kg)
k = 100;              % Spring stiffness
damping = 0.05;       % Damping coefficient
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
mothership_Kp = 10; 
mothership_Kd = 5;
mothership_max_velocity = 12; % Velocity limit

% Initial states
minion_l = 0.2;
minion_I = 0.5*masses(N)*minion_l^2;        % Moment of inertia (kg·m²)
minion_max_thrust = 10;      % Maximum thrust (N)
minion_theta_record = zeros(length(time), 1);  % Initial pitch angle (radians)
minion_theta = deg2rad(-30);
minion_omega = 0;
minion_Kp = diag([5, 2]);
minion_Kd = diag([10, 5]);
minion_thrust = zeros(length(0:dt:total_time), 1);
minion_error_record = zeros(length(time), 2);
minion_theta_desired_record = zeros(length(time), 1);
%-----------------------------------------------------------------------------------------------------------
% Simulation loop
for t = 1:length(time)

    % Store current positions and velocity
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
        
        % Damping + Gravity force
        forces(i, :) = forces(i, :) - damping * velocities(i, :);
        forces(i, 2) = forces(i, 2) - masses(i) * g;
    end
    
    % Mothership PD control
    mothership_error = mothership_target_position - positions(1, :);
    mothership_error_dot = -velocities(1, :);
    pd_force = mothership_Kp * mothership_error + mothership_Kd * mothership_error_dot + [0, masses(1) * g];
    forces(1, :) = forces(1, :) + pd_force;
    
    % % Apply PD controller for the minion (last particle)
    % target_position_minion = positions(1, 1);           % Target x is the mothership's current x
    % error_minion_x = target_position_minion - positions(N, 1); % Proportional term
    % error_dot_minion_x = -velocities(N, 1);               % Derivative term
    % pd_force_minion_x = Kp_minion * error_minion_x + Kd_minion * error_dot_minion_x;
    % pd_force_minion_x = min(max(pd_force_minion_x, -max_force_minion), max_force_minion);
    % forces(N, 1) = forces(N, 1) + pd_force_minion_x;  % Apply the constrained horizontal thrust
    % minion_thrust(t) = pd_force_minion_x;

    % Minion PD Control
    minion_state_current = [positions(N, 1); minion_theta];
    minion_state_desired = [positions(1, 1); 0];
    minion_error = minion_state_desired - minion_state_current;
    minion_theta_desired_record(t) = minion_state_desired(2);
    minion_error_record(t, :) = minion_error;

    minion_state_dot_current = [velocities(N, 1); minion_omega];
    % minion_state_dot_desired = [velocities(1, 1); 0];
    minion_state_dot_desired = [0; 0];
    minion_error_dot = minion_state_dot_desired - minion_state_dot_current;

    % Torque-ish
    control_input = (minion_Kp * minion_error) + (minion_Kd * minion_error_dot);
    % if (minion_theta < 0.01)
    %     minion_theta = 0;
    % end
    thrust = masses(N) * (control_input(1)) / sin(minion_theta);
    if thrust > 10000
        thrust = 0;
    end
    thrust = min(max(thrust, 0), minion_max_thrust);

    thrust_force = thrust * [sin(minion_theta), cos(minion_theta)];
    forces(N, :) = forces(N, :) + thrust_force;
    minion_thrust(t) = thrust;

    torque = minion_I * control_input(2);

    minion_alpha = torque / minion_I;
    minion_omega = minion_omega + minion_alpha * dt;
    minion_theta = minion_theta + minion_omega * dt;
    minion_theta_record(t) = minion_theta;

    % Update positions and velocities using Euler integration
    for i = 1:N
        accelerations = forces(i, :) / masses(i);
        velocities(i, :) = velocities(i, :) + accelerations * dt;
        
        % Limit velocity of the mothership (particle 1)
        if i == 1
            velocity_magnitude = norm(velocities(1, :));
            if velocity_magnitude > mothership_max_velocity
                velocities(1, :) = velocities(1, :) * (mothership_max_velocity / velocity_magnitude);
            end
        end
        
        % Update positions
        positions(i, :) = positions(i, :) + velocities(i, :) * dt;
    end
end
%-----------------------------------------------------------------------------------------------------------
% Position tracking
mothership_positions_x = squeeze(trajectory(:, 1, 1));
minion_positions_x = squeeze(trajectory(:, end, 1));
time = 0:dt:total_time;
figure;
plot(time, minion_positions_x, 'y', 'LineWidth', 2, 'DisplayName', 'Minion');
hold on;
plot(time, mothership_positions_x, 'k', 'LineWidth', 2, 'DisplayName', 'Mothership');
title('End Particle Position vs Time (no wind, no minion thrust');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;
hold off;
%-----------------------------------------------------------------------------------------------------------
% Plot thrust
figure;
plot(time, minion_thrust, 'b', 'LineWidth', 2);
title('Thrust Forces Over Time');
xlabel('Time (s)');
ylabel('Thrust Force (N)');
grid on;
hold off;
%----------------------------------------------------------------------------------------------
% Plot pitch angle over time
figure;
plot(time, rad2deg(minion_theta_record), 'g', 'LineWidth', 2, 'DisplayName', 'Pitch Angle');
title('Pitch Angle Over Time');
xlabel('Time (s)');
ylabel('Pitch Angle (degrees)');
legend;
grid on;
%-----------------------------------------------------------------------------------------------------------
% Create a figure with two subplots
figure;

% Subplot 1: Position error (X)
subplot(2, 1, 1);
plot(time, minion_error_record(:, 1), 'b', 'LineWidth', 2, 'DisplayName', 'X Position Error');
hold on;
% plot(time, ones(size(time)), '--r', 'LineWidth', 2, 'DisplayName', 'Limit');
% plot(time, ones(size(time)), '--r', 'LineWidth', 2, 'DisplayName', 'Limit');
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
plot(time, minion_theta_record/pi*180, 'r', 'LineWidth', 2, 'DisplayName', 'Pitch');
title('Minion Pitch Error Over Time');
xlabel('Time (s)');
ylabel('Pitch Error (degrees)');
grid on;
legend();
hold off;
%%
% Calculate velocities from trajectory
mother_ship_velocity_traj = diff(squeeze(trajectory(:, 1, :))) / dt; % Velocity of the mother ship
minion_velocity_traj = diff(squeeze(trajectory(:, end, :))) / dt;    % Velocity of the minion

% Pad the velocities to match the time vector length (since diff reduces the size by 1)
mother_ship_velocity_traj = [mother_ship_velocity_traj; [0, 0]];
minion_velocity_traj = [minion_velocity_traj; [0, 0]];

% % Initialize video writer
% video = VideoWriter('8-torque_thrust_10.mp4', 'MPEG-4');
% video.FrameRate = 20; % Set desired frame rate
% open(video);

% Visualization
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen figure
for t = 1:10:length(time)
    % Plot the rope
    plot(squeeze(trajectory(t, :, 1)), squeeze(trajectory(t, :, 2)), '-o', 'LineWidth', 2);
    hold on;
    plot(squeeze(trajectory(t, 1, 1)), squeeze(trajectory(t, 1, 2)), 'ko', 'LineWidth', masses(1), 'DisplayName', 'Mothership');
    plot(squeeze(trajectory(t, N, 1)), squeeze(trajectory(t, N, 2)), 'yo', 'LineWidth', masses(1), 'DisplayName', 'Minion');

    % Plot the target point for the mothership
    plot(mothership_target_position(1), mothership_target_position(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Initial Target');

    % Plot formatting
    axis equal;
    xlim([0, 150]);
    ylim([-10, mother_ship_position(2) + 20]);
    title(sprintf('Time: %.2f s', time(t)));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    drawnow;

    % % Capture the frame for video
    % frame = getframe(gcf);
    % writeVideo(video, frame); % Write the captured frame to the video

    hold off;
end

% Close the video file
close(video);
