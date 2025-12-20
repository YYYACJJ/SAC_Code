%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-500, 500]

% Parameter initialization
dt = 0.001;        % Simulation time step
t_end = 5;         % Simulation end time
t = 0:dt:t_end;    % Time vector

% Initial system state
x = [0; 0];        % Initial state vector
x_ref = @(t) 0.05 + (t > 2) * 0.02 * sin(t*pi);   % Reference position
dx_ref = @(t) (t > 2) * 0.02 * pi * cos(t*pi);   % Reference velocity

% Definition of nonlinear system dynamics and gain function
f = @(x, t) (t <= 1) * (-2*x(1) - cos(t)*x(2)^3) ...
          + (t > 1)  * 200 * (cos(t)^2*x(1) + 2*sin(t)*cos(t));  % Nonlinear function f(x,t)
g = @(x) 1 + x(1);                                 % Gain function g(x)
U = @(x,t) 0.5*sin(t)*x(1) + 0.6*cos(t)*x(2);      % Additional nonlinear term
integral_error = 0;                                % Integral error initialization

% Data storage
x_trajectory = zeros(2, length(t));   % State trajectory
control_input = zeros(1, length(t));  % Control input
s_trajectory = zeros(1, length(t));   % Sliding variable trajectory
previous_error = 0.05;                % Initial tracking error

% ====== 2) Pre-generate random variables for each 1 s interval (uniformly distributed in [-1, 1]) ======
N  = floor(t_end) + 2;                 % Number of intervals (extra one to avoid index overflow)
U1 = 2*rand(N,1) - 1;                  % Random term used in disturbance d4
U2 = 2*rand(N,1) - 1;                  % Random term used in disturbance d4
U3 = 2*rand(N,1) - 1;                  % Random term used in disturbance d4

% ====== 3) Utility function: return the random value of the current 1 s interval (piecewise constant) ======
stepRand = @(t,U) U(min(floor(t)+1, numel(U)));

% Persistent disturbance
d4 = @(x,t) -10*sin(t)*cos(t) ...
            + 10*stepRand(t,U1)*x(2) ...
            + 5*stepRand(t,U2)*x(1) ...
            + 2*stepRand(t,U3);

% ===================== Main simulation loop =====================
% PID controller parameters
Kp = 9000;        % Proportional gain
Ki = 200000;      % Integral gain
Kd = 100;         % Derivative gain

for k = 1:length(t)

    % Persistent disturbance
    d = d4(x, t(k));

    % Tracking error computation
    error = x_ref(t(k)) - x(1);
    e  = error;                     % Position error
    de = dx_ref(t(k)) - x(2);       % Velocity error

    % PID error update
    integral_error = integral_error + error * dt;
    derivative_error = (error - previous_error) / dt;

    % Auxiliary sliding-like variables
    s = Kp * error / Kd + Ki * integral_error / Kd + derivative_error;

    % PID control law
    u = Kp * error + Ki * integral_error + Kd * derivative_error;

    % Update previous error
    previous_error = error;

    % System dynamics
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d + U(x, t(k));

    % State update using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;

    % Data storage
    x_d(:, k) = x_ref(t(k));             % Reference trajectory
    x_trajectory_PID(:, k) = x;          % System state
    XW_e_PID(k) = e;                     % Position error
    XW_de_PID(k) = de - 2*sin(t(k));     % Modified velocity error
    x_error_PID(:, k) = error;           % Tracking error
    control_input_PID(:, k) = u;         % Control input
    error_s1_PID(k) = derivative_error + Kp * error / Kd;  % Sliding variable s1
    error_s0_PID(k) = error + Ki * integral_error / Kd;    % Sliding variable s0
end

% ===================== Plotting =====================

% Sliding variable s0 versus time
figure;
plot(t, error_s0_PID(1, :), 'r'); hold on;
xlabel('t', 'FontSize', 16);
ylabel('s0', 'FontSize', 16);
ylim([-0.1, 0.1]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Sliding variable s1 versus time
figure;
plot(t, error_s1_PID(1, :), 'r'); hold on;
xlabel('t', 'FontSize', 16);
ylabel('s1', 'FontSize', 16);
ylim([-0.1, 0.1]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Control input versus time
figure;
plot(t, control_input_PID(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Tracking error versus time
figure;
plot(t, x_error_PID(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Position trajectory versus time
figure;
plot(t, x_trajectory_PID(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
