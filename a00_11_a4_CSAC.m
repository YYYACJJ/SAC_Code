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
%Torque limited to [-10, 10]

% Parameter initialization
dt = 0.001;        % Simulation time step
t_end = 5;         % Simulation end time
t = 0:dt:t_end;    % Time vector
close all;         % Close all figure windows

% Initial system state
x = [0; 0];        % Initial state vector [x1, x2]
x_ref = @(t) 0.05; % Reference position
dx_ref = @(t) 0;   % Reference velocity

% Definition of nonlinear system dynamics and gain function
f = @(x, t) -2*x(1) - cos(t)*x(2)^3;   % Nonlinear function f(x,t)
g = @(x) 1 + x(1);                     % Gain function g(x)
U = @(x,t) 0.5*sin(t)*x(1) + 0.6*cos(t)*x(2); % Additional nonlinear term

% Data storage
x_trajectory = zeros(2, length(t));    % State trajectory
control_input = zeros(1, length(t));   % Control input
s_trajectory = zeros(1, length(t));    % Sliding surface trajectory
previous_error = 0.04;                 % Initial error value

% ====== 2) Pre-generate random variables for each 1 s interval (uniformly distributed in [-1, 1]) ======
N  = floor(t_end) + 2;                 % Number of intervals (extra one to avoid index overflow)
U1 = 2*rand(N,1) - 1;                  % Random disturbance term
U2 = 2*rand(N,1) - 1;                  % Random disturbance term
U3 = 2*rand(N,1) - 1;                  % Random disturbance term
U4 = 2*rand(N,1) - 1;                  % Random disturbance term
U5 = 2*rand(N,1) - 1;                  % Random disturbance term
U6 = 2*rand(N,1) - 1;                  % Random disturbance term

% ====== 3) Utility function: return the random value of the current 1 s interval (piecewise constant) ======
stepRand = @(t,UT) UT(min(floor(t)+1, numel(UT)));

% Time-varying disturbance
d1 = @(x,t) 50*(sin(t)*cos(t)*x(1) + cos(t)*x(2)) ...
            + 5*stepRand(t,U1)*x(2) ...
            + 10*stepRand(t,U2)*x(1) ...
            + stepRand(t,U3);

% Persistent disturbance
d2 = @(x,t) -sin(t)*cos(t) ...
            + stepRand(t,U4)*x(2) ...
            + stepRand(t,U5)*x(1) ...
            + 0.2*stepRand(t,U6);

% Controller parameters
k_cc = 1000;       % Control gain
k_c1 = 0.015;      % Sliding surface parameter

for k = 1:length(t)

    % Persistent disturbance
    d = d2(x, t(k));
    
    % Time-varying disturbance activation
    if k > 500 && k < 1500
        d = d2(x, t(k)) + d1(x, t(k));
    end

    % Sliding surface construction
    e  = x_ref(t(k)) - x(1);        % Position error
    de = dx_ref(t(k)) - x(2);       % Velocity error

    game1 = fcn(e, k_c1);           % Nonlinear compensation term
    s = de - game1;                 % Sliding surface
    u = k_cc * s;                   % Control law

    % System dynamics
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d + U(x, t(k));

    % State update using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;

    % Data storage
    x_trajectory_CSAC(:, k) = x;    % State trajectory
    control_input_CSAC(k) = u;      % Control input
    x_error_CSAC(k) = e;            % Position error
    x_derror_CSAC(k) = de;          % Velocity error
end

% Phase portrait: e vs. de
figure;
plot(x_error_CSAC(1, :), x_derror_CSAC(1, :), 'r'); hold on;
xlabel('e', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Velocity error versus time
figure;
plot(t, x_derror_CSAC(1, :), 'r'); hold on;
xlabel('t', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
ylim([-0.1, 0.1]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Control input versus time
figure;
plot(t, control_input_CSAC(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Position error versus time
figure;
plot(t, x_error_CSAC(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Position trajectory versus time
figure;
plot(t, x_trajectory_CSAC(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.001, 0.08]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% CSAC
function game = fcn(e, k_L)
    game = -k_L * ((pi/2 - atan(100*abs(e))) * atan(abs(e)*1000) + 1e-8) ...
           * atan(10000*e);
end
