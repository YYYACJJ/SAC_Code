%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common pitfalls in SAF function design:
% The function is not designed according to the error range, which causes
% the SAC to fail to track the function.
%
% We hope this code is helpful and inspires new ideas.
% If it is useful to you, please consider citing our paper:
% "Attractive Surface-Based State-Attracted Control for Uncertain Nonlinear Systems"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-15, 15]


% Parameter initialization
dt = 0.001;       % Simulation time step
t_end = 5;        % End time of simulation
t = 0:dt:t_end;   % Time vector
x = 0;            % Initial position
x_ref = @(t)  0.05 + (t > 2) * 0.02 * sin(t*pi);   % Reference position
dx_ref = @(t) (t > 2) * 0.02 * pi * cos(t*pi);    % Reference velocity

% Initialization of integral error and previous error
integral_error = 0;
previous_error = 0;

% Storage of results
x_trajectory = zeros(1, length(t));   % State trajectory
control_input = zeros(1, length(t));  % Control input
x_error = zeros(1, length(t));        % Tracking error

% Definition of nonlinear function and gain function
f = @(x, t) 0.9*cos(x) + 0.3*sign(x);  % Nonlinear function f(x,t)
g = @(x) 1;                            % Gain function g(x)

U = 0.2*sign(x) + 0.1*sin(x);          % Additional nonlinear term

% ====== 2) Pre-generate random variables for each 1 s interval (uniformly distributed in [-1, 1]) ======
N  = floor(t_end) + 2;                 % Number of intervals (extra one to avoid index overflow)
U1 = 2*rand(N,1) - 1;                  % Random disturbance term
U2 = 2*rand(N,1) - 1;                  % Random disturbance term
U3 = 2*rand(N,1) - 1;                  % Random disturbance term
U4 = 2*rand(N,1) - 1;                  % Random disturbance term

% ====== 3) Utility function: return the random value of the current 1 s interval (piecewise constant) ======
stepRand = @(t,U) U(min(floor(t)+1, numel(U)));

% ====== 4) Disturbance functions (random terms updated every 1 s) ======
d1 = @(x,t) 2*sin(0.5*t) + cos(x) + 2*sign(0.4*t) + stepRand(t,U1)*x + stepRand(t,U2);
d2 = @(x,t) 2*sin(0.5*t) + cos(x) + 2*sign(0.4*t) + stepRand(t,U3)*x + stepRand(t,U4);

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
    integral_error = integral_error + error * dt;
    derivative_error = (error - previous_error) / dt;
    previous_error = error;

    % Controller parameters
    k_cf = 180;
    k_f1 = 170;

    % Nonlinear compensation term
    game = fcn(integral_error, k_f1);

    % Control law
    u = k_cf * (error - game);
    control_input(k) = u;

    % System dynamics
    dx = f(x, t(k)) + g(x) * u + U + d;

    % State update using Euler method
    x = x + dx * dt;

    % Store simulation data
    x_trajectory_FSAC(:, k) = x;
    x_error_FSAC(:, k) = error;
    control_input_FSAC(k) = u;
    dx_FSAC(k) = dx;
end

% Plot: control input
figure;
plot(t, control_input_FSAC(1, :), 'r', 'DisplayName', 'SAC');
xlabel('Time (s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
legend('show');

% Plot: tracking error
figure;
plot(t, x_error_FSAC(1, :), 'r', 'DisplayName', 'SAC');
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
plot(t, x_trajectory_FSAC(1, :), 'r', 'DisplayName', 'SAC');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.1, 0.3]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
legend('show');

% FSAF
function game = fcn(e, k_L)
x = e;
game1 = -k_L * (atan(abs(x)) + 1e-8) * atan(10000 * x);
game = game1;
end
