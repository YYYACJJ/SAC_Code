
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-15, 15]

% Parameter initialization
dt = 0.001;       % Simulation time step
t_end = 5;        % Simulation end time
t = 0:dt:t_end;   % Time vector

% Initial system state
x = 0;            % Initial position

% Reference trajectory
x_ref  = @(t) 0.05 + (t > 2) * 0.02 * sin(t*pi);     % Reference position
dx_ref = @(t) (t > 2) * 0.02 * pi * cos(t*pi);      % Reference velocity

% Initialization of integral error and previous error
integral_error = 0;
previous_error = 0;

% Storage variables
x_trajectory = zeros(1, length(t));   % State trajectory
control_input = zeros(1, length(t));  % Control input
x_error = zeros(1, length(t));        % Tracking error

% Definition of nonlinear system function and gain function
f = @(x, t) 0.9*cos(x) + 0.3*sign(x);  % Nonlinear function f(x,t)
g = @(x) 1;                            % Gain function g(x)

% PID controller parameters
integral_error = 0;
previous_error = 0;
Kp = 180;      % Proportional gain (diagonal)
Ki = 2500;     % Integral gain (diagonal)
Kd = 0.01;     % Derivative gain (diagonal)

% Additional nonlinear term
U=@(x)0.2*sign(x)+0.1*sin(x);

% ====== 2) Pre-generate random variables for each 1 s interval (uniform distribution in [-1, 1]) ======
N  = floor(t_end) + 2;                 % Number of intervals (extra one to avoid index overflow)
U1 = 2*rand(N,1) - 1;                  % Random disturbance term
U2 = 2*rand(N,1) - 1;                  % Random disturbance term
U3 = 2*rand(N,1) - 1;                  % Random disturbance term
U4 = 2*rand(N,1) - 1;                  % Random disturbance term

% Utility function: return the random value of the current 1 s interval (piecewise constant)
stepRand = @(t,UT) UT(min(floor(t)+1, numel(UT)));

% ====== 4) Disturbance functions (random terms updated every 1 s) ======
d1 = @(x,t) 2*sin(0.5*t) + cos(x) + 2*sign(0.4*t) ...
            + stepRand(t,U1)*x + stepRand(t,U2);
d2 = @(x,t) 2*sin(0.5*t) + cos(x) + 2*sign(0.4*t) ...
            + stepRand(t,U3)*x + stepRand(t,U4);

for k = 1:length(t)

    % Default disturbance
    d = d2(x, t(k));

    % Time-varying disturbance switching
    if k > 1000 && k < 1500
        d = d1(x, t(k)) + d2(x, t(k));
    elseif k > 3500 && k < 4500
        d = d1(x, t(k)) + d2(x, t(k));
    end

    % Tracking error computation
    error = x_ref(t(k)) - x;

    % PID error update
    integral_error = integral_error + error * dt;
    derivative_error = (error - previous_error) / dt;

    % PID control law
    u = Kp * error + Ki * integral_error + Kd * derivative_error;

    % Update previous error
    previous_error = error;

    % System dynamics
    dx = f(x, t(k)) + g(x) * u +U(x) + d;

    % State update using Euler method
    x = x + dx * dt;

    % Data storage
    x_d(:, k) = x_ref(t(k));            % Reference trajectory
    x_trajectory_PID(:, k) = x;         % System state
    x_error_PID(:, k) = error;          % Tracking error
    control_input_PID(k) = u;           % Control input
    dx_PID(k) = dx;                     % State derivative
end

% Plot: control input
figure;
plot(t, control_input_PID(1, :), 'r', 'DisplayName', 'PID');
xlabel('Time (s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
legend('show');

% Plot: tracking error
figure;
plot(t, x_error_PID(1, :), 'r', 'DisplayName', 'PID');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.1, 0.1]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
legend('show');

% Plot: state trajectory
figure;
plot(t, x_trajectory_PID(1, :), 'r', 'DisplayName', 'PID');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.1, 0.3]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
legend('show');
