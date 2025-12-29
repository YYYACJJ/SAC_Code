%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Due to the existence of random terms, the system may get out of control 
%and it may be necessary to try several times
%Torque limited to [-2000, 2000]


% Parameter initialization
dt = 0.001;       % Simulation time step
t_end = 5;        % Simulation end time
t = 0:dt:t_end;   % Time vector
n = 0;
Da = 2;
Dt = 0.5;

% RBF (Radial Basis Function) network parameters
c_a = 100;
HidLayNuma = 5;          % Number of hidden-layer neurons
OutLayNum = 1;           % Number of output-layer neurons
gamma_a = 80000;         % Learning rate (adaptation gain)
eta_a = 1950;            % Switching gain
c_i_a = 8*[-2 -1  0  1  2
           -2 -1  0  1  2]; % RBF centers
b_i_a = 0.1;             % RBF width
hat_W1_a = 0.08*ones(HidLayNuma, OutLayNum); % Initial NN weights

x = [0; 0];
xr_a  = @(t) 10 + (t > 2) * 4*sin(t*pi);     % Reference position
dxr_a = @(t)      (t > 2) * 4*pi*cos(t*pi);  % Reference velocity

% Definitions of nonlinear dynamics and input gain
f = @(x, t) ((t <= 1) * (-2*x(1) - cos(t)*x(2)*x(2)*x(2)) + ...
            (t > 1)  * 2200*(cos(t)*cos(t)*x(1) + 2*sin(t)*cos(t))); % Nonlinear term f(x,t)
g = @(x) 1 + x(1);                                                % Input gain g(x)
U = @(x,t) 0.5*sin(t)*x(1) + 0.6*cos(t)*x(2);                     % Additional known term
integral_error = 0;

% ====== 2) Pre-generate random variables for each 1-second interval (uniform in [-1,1]) ======
N  = floor(t_end) + 2;    % Number of intervals (one extra to avoid out-of-bound indexing)
U1 = 2*rand(N,1) - 1;     % Random term for disturbance component
U2 = 2*rand(N,1) - 1;     % Random term for disturbance component
U3 = 2*rand(N,1) - 1;     % Random term for disturbance component

% ====== 3) Utility function: given t, return the random value of the current 1-second interval (piecewise constant) ======
stepRand = @(t,U) U( min(floor(t)+1, numel(U)) );

d4 = @(x,t) -100*sin(t)*cos(t) + 100*stepRand(t,U1)*x(2) + ...
            50*stepRand(t,U2)*x(1) + 20*stepRand(t,U3); % Persistent disturbance

% Main simulation loop
tic  % Start timing
for k = 1:length(t)

    %%%%%%%%%%%% Disturbance %%%%%%%%%%%%%%%%%
    d = d4(x, t(k));   % Persistent disturbance

    e  = xr_a(t(k))  - x(1);
    de = dxr_a(t(k)) - x(2);
    s_a = c_a*e + de;  % Sliding surface

    % Compute RBF activations
    for j = 1:HidLayNuma
        h1a(j,:) = exp( -norm([x(1); x(2)] - c_i_a(:,j))^2/(2*b_i_a^2) );
    end

    % NN estimation and weight update
    hat_fx1_a = hat_W1_a' * h1a;
    hat_W1_a  = hat_W1_a + dt * (gamma_a * s_a(1) * h1a);

    % Control law
    u = eta_a * sign(s_a) + hat_fx1_a;

    % System dynamics
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d + U(x, t(k));

    % Update states using the Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;

    % Store trajectories
    x_trajectory_NNSMC(:, k) = x;
    control_input_NNSMC(k)  = u;
    x_error_NNSMC(k)        = e;
    x_derror_NNSMC(k)       = de;

end

elapsed_time = toc;  % Stop timing and return elapsed time
fprintf('Total simulation time: %.4f s\n', elapsed_time);

% Plotting
figure;
plot(t, control_input_NNSMC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_error_NNSMC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_trajectory_NNSMC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

function S = skew(v)
% Input: v is a 3×1 vector; Output: S is a 3×3 skew-symmetric matrix
    S = [   0    -v(3)   v(2);
           v(3)   0    -v(1);
          -v(2)  v(1)    0   ];
end
