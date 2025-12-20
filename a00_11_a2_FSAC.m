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
%Torque limited to [-100, 100]

% Parameter initialization
dt = 0.001;        % Simulation time step
t_end = 5;         % Simulation end time
t = 0:dt:t_end;    % Time vector
close all;         % Close all figure windows

% Initial system state
x = [0; 0];        % Initial state
x_ref = @(t) 0.05 + (t > 2) * 0.02 * sin(t*pi);   % Reference position
dx_ref = @(t) (t > 2) * 0.02 * pi * cos(t*pi);    % Reference velocity

% Definition of nonlinear function and gain function
f = @(x, t) -2*x(1) - cos(t)*x(2)*x(2)*x(2);   % Nonlinear function f(x,t)
g = @(x) 1 + x(1);                             % Gain function g(x)
U = @(x,t) 0.5*sin(t)*x(1) + 0.6*cos(t)*x(2);  % Known auxiliary input

% Data storage
x_trajectory = zeros(2, length(t));   % State trajectory
control_input = zeros(1, length(t));  % Control input
s_trajectory = zeros(1, length(t));   % Sliding surface trajectory
previous_error = 0.04;

% ====== 2) Pre-generate random variables for each 1-second interval (uniform in [-1,1]) ======
N  = floor(t_end) + 2;    % Number of intervals (extra to avoid index overflow)
U1 = 2*rand(N,1) - 1;     % Random term for disturbance component
U2 = 2*rand(N,1) - 1;     
U3 = 2*rand(N,1) - 1;    
U4 = 2*rand(N,1) - 1;    
U5 = 2*rand(N,1) - 1;     
U6 = 2*rand(N,1) - 1;     
U7 = 2*rand(N,1) - 1;     

% ====== 3) Utility function: return the random value of the current 1-second interval ======
stepRand = @(t,U) U(min(floor(t)+1, numel(U)));

d3 = @(x,t) 70 * stepRand(t,U7);   % Impulsive disturbance
d2 = @(x,t) 100*(sin(t)*cos(t)*x(1) + cos(t)*x(2)) ...
          + 2*stepRand(t,U1)*x(2) ...
          + 4*stepRand(t,U2)*x(1) ...
          + stepRand(t,U3);        % Time-varying disturbance

% Persistent disturbance (multiple candidate forms)
d1=@(x,t)-sin(t)*cos(t)+2*stepRand(t,U4)*x(2)+5*stepRand(t,U5)*x(1)+stepRand(t,U6);
% d1=@(x,t)-30*sin(t)+40*cos(t)+5*stepRand(t,U4)*x(2)+11*stepRand(t,U5)*x(1)+stepRand(t,U6);
% d1=@(x,t)-10*cos(t)+20*sin(x(1))*x(2)+10*stepRand(t,U4)*x(2)+5*stepRand(t,U5)*x(1)+2*stepRand(t,U6);
% d1=@(x,t)-100*sin(x(1))+20+50*cos(t)+8*stepRand(t,U4)*x(2)+4*stepRand(t,U5)*x(1)+4*stepRand(t,U6);
% d1=@(x,t)90*sin(x(1))+80*sin(t)*cos(t)+10*stepRand(t,U4)*x(2)+2*stepRand(t,U5)*x(1)+5*stepRand(t,U6);
% d1=@(x,t)-100*sin(t)*cos(t)*sin(x(1))+20*stepRand(t,U4)*x(2)+8*stepRand(t,U5)*x(1)+10*stepRand(t,U6);
% d1=@(x,t)-80*sin(x(1))-80*sin(t)*cos(t)+10*cos(t)+5*stepRand(t,U4)*x(2)+5*stepRand(t,U5)*x(1)+2*stepRand(t,U6);
% d1=@(x,t)5*cos(x(2))+10+40*sin(t)*cos(t)+7*stepRand(t,U4)*x(2)+4*stepRand(t,U5)*x(1)+6*stepRand(t,U6);
% d1=@(x,t)-10*sin(t)*cos(t)*sin(x(1))-20+8*stepRand(t,U4)*x(2)+3*stepRand(t,U5)*x(1)+3*stepRand(t,U6);

% Control parameters
k_cf = 240; % Ensures system stability and attraction strength
k_f1 = 0.17; % Adjusts SAF intensity

tic   % Start timing
for k = 1:length(t)

    d = d1(x, t(k));   % Persistent disturbance
    if k > 500 && k < 1500
        d = d1(x, t(k)) + d2(x, t(k));
    elseif k > 2500 && k < 3500
        d = d1(x, t(k)) + d3(x, t(k));
    end

    %Surface construction
    e  = x_ref(t(k)) - x(1);
    de = dx_ref(t(k)) - x(2);
    game1 = fcn(e, k_f1);
    s = de - game1;
    u = k_cf * s;

    % System dynamics
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d + U(x, t(k));

    % State update using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;

    % Store state trajectory and sliding surface data
    x_trajectory_FSAC(:, k) = x;
    control_input_FSAC(k) = u;
    x_error_FSAC(k) = e;
    x_derror_FSAC(k) = de;
    error_s1_FSAC(k) = de - game1;
end

elapsed_time = toc;   % End timing and return elapsed time
fprintf('Total simulation time: %.4f seconds\n', elapsed_time);

% Plot: phase portrait
figure;
plot(x_error_FSAC(1,:), x_derror_FSAC(1,:), 'r'); hold on;
xlabel('e', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Plot: control input
figure;
plot(t, control_input_FSAC(1,:), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Plot: tracking error
figure;
plot(t, x_error_FSAC(1,:), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Plot: position trajectory
figure;
plot(t, x_trajectory_FSAC(1,:), 'r');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.001, 0.08]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% FSAF
function game = fcn(e, k_L)
    x = e;
    game1 = -k_L * (atan(abs(3000*x)) + 1e-8) .* atan(10000*x);
    game = game1;
end
