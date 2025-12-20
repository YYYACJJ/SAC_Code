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


% Parameter Initialization
dt = 0.001;       % Simulation time step
t_end = 5;      % Simulation end time
t = 0:dt:t_end;  % Time vector
x = 0;  % Initial position
x_ref = @(t)  0.05 + (t > 2) *0.02*sin(t*pi);   % Target position
dx_ref =@(t)  (t > 2) *0.02*pi*cos(t*pi);


% Integral error and error initialization
integral_error = 0;
previous_error = 0;

% Result storage
x_trajectory = zeros(1, length(t));  % Record state
control_input = zeros(1, length(t)); % Record control input
x_error = zeros(1, length(t));

% Nonlinear function and gain function definitions
f = @(x, t) 0.9*cos(x)+0.3*sign(x);   % Nonlinear function f(x,t)
g = @(x) 1 ;               % Gain function g(x)

U=0.2*sign(x)+0.1*sin(x);
% ====== 2) Pre-generate random variables for each 1s interval (uniformly distributed in [-1,1]) ======
N  = floor(t_end) + 2;       % Number of intervals (extra one for safety)
U1 = 2*rand(N,1) - 1;           % Random term
U2 = 2*rand(N,1) - 1;           % Random term
U3 = 2*rand(N,1) - 1;           % Random term 
U4 = 2*rand(N,1) - 1;           % Random term
stepRand = @(t,U) U( min(floor(t)+1, numel(U)) );

% ====== 4) Disturbance functions (adding random terms that change every 1s to the original expression) ======
d1 = @(x,t) 2*sin(0.5*t) + cos(x) + 2*sign(0.4*t) + stepRand(t,U1)*x+stepRand(t,U2);
d2 = @(x,t) 2*sin(0.5*t) + cos(x) + 2*sign(0.4*t) + stepRand(t,U3)*x+stepRand(t,U4);
for k = 1:length(t)

    d=d2(x,t(k));
    if k>1000 && k<1500
        d=d1(x,t(k))+d2(x,t(k));
    elseif k>3500 && k<4500
        d=d1(x,t(k))+d2(x,t(k));     
    end

    error = x_ref(t(k))-x(1);
    integral_error = integral_error + error * dt;
    derivative_error = (error - previous_error)/dt ;
    previous_error = error;
    
    k_cc=180;
    k_c1=6000;
    game=fcn(integral_error,k_c1);
    u = k_cc*(error-game);
    control_input(k) = u;
    
    
    
    % System dynamic equation
    dx = f(x, t(k)) + g(x) * u +U+ d;
    
    % Update state using Euler method
    x = x + dx  * dt;
    
    % Store state trajectory
    x_trajectory_CSAC(:, k) = x;
    x_error_CSAC(:, k) = error;
    control_input_CSAC(k) = u;
    dx_CSAC(k) = dx;
end

% Plotting
figure;
plot(t, control_input_CSAC(1, :), 'r', 'DisplayName', 'SAC');
xlabel('Time(s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;
legend('show');


figure;
plot(t, x_error_CSAC(1, :), 'r', 'DisplayName', 'SAC');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.1, 0.1]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;
legend('show');

figure;
plot(t, x_trajectory_CSAC(1, :), 'r', 'DisplayName', 'SAC');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.1, 0.3]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;
legend('show');


% CSAF
function game=fcn(e,k_L)
x = e;
game=-k_L*(pi/2-atan(0.01*x*x))*(atan(1 * x*x) + 0.00000001) * atan(10000 * x);

end