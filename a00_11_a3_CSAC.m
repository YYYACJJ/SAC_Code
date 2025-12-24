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
%Torque limited to [-2000, 2000]

% Parameter initialization
dt = 0.001;        % Simulation time step
t_end = 5;         % End time of simulation
t = 0:dt:t_end;    % Time vector

% Initial system state
x = [0; 0];        % Initial state [x1, x2]
x_ref  = @(t) 10 + (t > 2) * 4*sin(t*pi);     % Reference position
dx_ref = @(t)      (t > 2) * 4*pi*cos(t*pi);  % Reference velocity

% Definitions of nonlinear dynamics and input gain
f = @(x, t) ((t <= 1) * (-2*x(1) - cos(t)*x(2)*x(2)*x(2)) + ...
            (t > 1)  * 2200*(cos(t)*cos(t)*x(1) + 2*sin(t)*cos(t))); % Nonlinear term f(x,t)
g = @(x) 1 + x(1);                                                % Input gain g(x)
U = @(x,t) 0.5*sin(t)*x(1) + 0.6*cos(t)*x(2);                     % Additional known term

% Data storage
x_trajectory = zeros(2, length(t));   % State trajectory
control_input = zeros(1, length(t));  % Control input
s_trajectory = zeros(1, length(t));   % Sliding surface trajectory

% ====== 2) Pre-generate random variables for each 1-second interval (uniform in [-1,1]) ======
N  = floor(t_end) + 2;      % Number of intervals (one extra to avoid out-of-bound indexing)
U1 = 2*rand(N,1) - 1;       % Random sequence
U2 = 2*rand(N,1) - 1;       % Random sequence
U3 = 2*rand(N,1) - 1;       % Random sequence

% ====== 3) Utility function: given t, return the random value of the current 1-second interval (piecewise constant) ======
stepRand = @(t,U) U( min(floor(t)+1, numel(U)) );

d4 = @(x,t) (-100*sin(t)*cos(t) + 100*stepRand(t,U1)*x(2) + ...
             50*stepRand(t,U2)*x(1) + 20*stepRand(t,U3)); % Persistent disturbance

% Main simulation loop
k_cc = 60;
k_c1 = 30;

for k = 1:length(t)

    d = d4(x, t(k)); % Persistent disturbance

    % Construct the sliding surface
    e  = x_ref(t(k))  - x(1);
    de = dx_ref(t(k)) - x(2);

    game1 = fcn(e, k_c1);
    s = de - game1;
    u = k_cc * s;

    % System dynamics
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d + U(x, t(k));

    % Update states using the Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;

    % Store trajectories
    x_trajectory_CSAC(:, k) = x;
    control_input_CSAC(k) = u;
    x_error_CSAC(k) = e;
    x_derror_CSAC(k) = de;

end

% Plotting
figure;
plot(x_error_CSAC(1, :), x_derror_CSAC(1, :), 'r'); hold on;
xlabel('e', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, control_input_CSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_error_CSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_trajectory_CSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

function game = fcn(e, k_L)
    game = -k_L * ((pi/2 - atan(1*abs(e))) * atan(5*abs(e)) + 1e-8) * atan(10000*e);
end


