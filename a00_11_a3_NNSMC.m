%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-500, 500]

% Parameter initialization
dt = 0.001;        % Simulation time step
t_end = 5;         % Simulation end time
t = 0:dt:t_end;    % Time vector
n = 0;
Da = 2;
Dt = 0.5;

% ===================== RBF neural network parameters =====================
c_a = 70;                     % Sliding surface gain
HidLayNuma = 5;               % Number of hidden-layer neurons
OutLayNum = 1;                % Number of output neurons
gamma_a = 60000;              % Adaptive learning rate
eta_a = 160;                  % Switching gain
c_i_a = 8 * [-2 -1  0  1  1    % RBF center vectors
             -2 -1  0  1  1];
b_i_a = 0.1;                  % RBF width
hat_W1_a = 0.08 * ones(HidLayNuma, OutLayNum);  % Initial adaptive weights

% Initial system state
x = [0; 0];

% Reference trajectory
xr_a  = @(t) 0.05 + (t > 2) * 0.02 * sin(t*pi);      % Reference position
dxr_a = @(t) (t > 2) * 0.02 * pi * cos(t*pi);       % Reference velocity

% Definition of nonlinear system dynamics and gain function
f = @(x, t) (t <= 1) * (-2*x(1) - cos(t)*x(2)^3) ...
          + (t > 1)  * 200 * (cos(t)^2*x(1) + 2*sin(t)*cos(t));  % Nonlinear function f(x,t)
g = @(x) 1 + x(1);                                  % Gain function g(x)
U = @(x,t) 0.5*sin(t)*x(1) + 0.6*cos(t)*x(2);       % Additional nonlinear term

integral_error = 0;

% ====== 2) Pre-generate random variables for each 1 s interval (uniformly distributed in [-1, 1]) ======
N  = floor(t_end) + 2;          % Number of intervals (extra one to avoid index overflow)
U1 = 2*rand(N,1) - 1;           % Random disturbance term
U2 = 2*rand(N,1) - 1;           % Random disturbance term
U3 = 2*rand(N,1) - 1;           % Random disturbance term

% ====== 3) Utility function: return the random value of the current 1 s interval (piecewise constant) ======
stepRand = @(t,U) U(min(floor(t)+1, numel(U)));

% Persistent disturbance
d4 = @(x,t) -10*sin(t)*cos(t) ...
            + 10*stepRand(t,U1)*x(2) ...
            + 5*stepRand(t,U2)*x(1) ...
            + 2*stepRand(t,U3);

% ===================== Main simulation loop =====================
tic;   % Start timing
for k = 1:length(t)

    % Persistent disturbance
    d = d4(x, t(k));

    % Tracking error and sliding surface
    e  = xr_a(t(k))  - x(1);      % Position error
    de = dxr_a(t(k)) - x(2);      % Velocity error
    s_a = c_a * e + de;           % Sliding surface

    % RBF neural network activation
    for j = 1:HidLayNuma
        h1a(j,:) = exp( -norm([x(1); x(2)] - c_i_a(:,j))^2 / (2*b_i_a^2) );
    end

    % Neural network output estimation
    hat_fx1_a = hat_W1_a' * h1a;

    % Adaptive weight update law
    hat_W1_a = hat_W1_a + dt * (gamma_a * s_a(1) * h1a);

    % Control law
    u = eta_a * sign(s_a) + hat_fx1_a;

    % System dynamics
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d + U(x, t(k));

    % State update using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;

    % Data storage
    x_trajectory_NNSMC(:, k) = x;   % State trajectory
    control_input_NNSMC(k) = u;     % Control input
    x_error_NNSMC(k) = e;           % Position error
    x_derror_NNSMC(k) = de;         % Velocity error
end

elapsed_time = toc;  % Stop timing and return elapsed time
fprintf('Total simulation time: %.4f seconds\n', elapsed_time);

% ===================== Plotting =====================

% Control input versus time
figure;
plot(t, control_input_NNSMC(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
ylim([-500, 500]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Tracking error versus time
figure;
plot(t, x_error_NNSMC(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Position trajectory versus time
figure;
plot(t, x_trajectory_NNSMC(1, :), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.001, 0.08]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;

% Skew-symmetric matrix function
function S = skew(v)
% Input v is a 3×1 vector; output S is a 3×3 skew-symmetric matrix
    S = [   0    -v(3)   v(2);
           v(3)   0    -v(1);
          -v(2)  v(1)    0   ];
end
