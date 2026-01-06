%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-100, 100]

% Parameter Initialization
dt = 0.001;        % Simulation time step
t_end = 5;        % Simulation end time
t = 0:dt:t_end;   % Time vector

close all;     % Close all figure windows

% System initial state
x = [0; 0];   % Initial state 
x_ref = @(t)  0.05 + (t > 2) *0.02*sin(t*pi);   % Target position
dx_ref =@(t)  (t > 2) *0.02*pi*cos(t*pi);
% Nonlinear function and gain function definitions
f = @(x, t) -2*x(1)-cos(t)*x(2)*x(2)*x(2);  % Nonlinear function f(x,t)
g = @(x) 1+x(1);               % Gain function g(x)
U=@(x,t) 0.5*sin(t)*x(1)+0.6*cos(t)*x(2);
integral_error = 0;

% Result storage
x_trajectory = zeros(2, length(t));  % Record state
control_input = zeros(1, length(t)); % Record control input
s_trajectory = zeros(1, length(t));  % Record sliding surface trajectory
previous_error=0.05;


% ====== 2) Pre-generate random variables for each 1s interval (uniformly distributed in [-1,1]) ======
N  = floor(t_end) + 2;       % Number of intervals (extra one for safety)
U1 = 2*rand(N,1) - 1;           % Random term
U2 = 2*rand(N,1) - 1;           % Random term
U3 = 2*rand(N,1) - 1;           % Random term
U4 = 2*rand(N,1) - 1;           % Random term
U5 = 2*rand(N,1) - 1;           % Random term
U6 = 2*rand(N,1) - 1;           % Random term
U7 = 2*rand(N,1) - 1;           % Random term
% ====== 3) A small utility: given t, returns the random value for the "current 1s interval" (piecewise constant) ======
stepRand = @(t,UT) UT(min(floor(t)+1, numel(UT)));


d3 = @(x,t) 70*stepRand(t,U7);  % Impulse disturbance
d2 = @(x,t) 100*(sin(t)*cos(t)*x(1)+cos(t)*x(2))+2*stepRand(t,U1)*x(2)+4*stepRand(t,U2)*x(1)+stepRand(t,U3);  % Time-varying function
d1=@(x,t)-sin(t)*cos(t)+2*stepRand(t,U4)*x(2)+5*stepRand(t,U5)*x(1)+stepRand(t,U6);%Continuous disturbance
% d1=@(x,t)-30*sin(t)+40*cos(t)+5*stepRand(t,U4)*x(2)+11*stepRand(t,U5)*x(1)+stepRand(t,U6);%Continuous disturbance
% d1=@(x,t)-10*cos(t)+20*sin(x(1))*x(2)+10*stepRand(t,U4)*x(2)+5*stepRand(t,U5)*x(1)+2*stepRand(t,U6);%Continuous disturbance
% d1=@(x,t)-100*sin(x(1))+20+50*cos(t)+8*stepRand(t,U4)*x(2)+4*stepRand(t,U5)*x(1)+4*stepRand(t,U6);%Continuous disturbance
% d1=@(x,t)90*sin(x(1))+80*sin(t)*cos(t)+10*stepRand(t,U4)*x(2)+2*stepRand(t,U5)*x(1)+5*stepRand(t,U6);%Continuous disturbance
% d1=@(x,t)-100*sin(t)*cos(t)*sin(x(1))+20*stepRand(t,U4)*x(2)+8*stepRand(t,U5)*x(1)+10*stepRand(t,U6);%Continuous disturbance
% d1=@(x,t)-80*sin(x(1))-80*sin(t)*cos(t)+10*cos(t)+5*stepRand(t,U4)*x(2)+5*stepRand(t,U5)*x(1)+2*stepRand(t,U6);
% d1=@(x,t)5*cos(x(2))+10+40*sin(t)*cos(t)+7*stepRand(t,U4)*x(2)+4*stepRand(t,U5)*x(1)+6*stepRand(t,U6);
% d1=@(x,t)-10*sin(t)*cos(t)*sin(x(1))-20+8*stepRand(t,U4)*x(2)+3*stepRand(t,U5)*x(1)+3*stepRand(t,U6);


% PID controller parameters
Kp = 1900; % Proportional gain (diagonal matrix)
Ki = 40000;   % Integral gain (diagonal matrix)
Kd = 100;   % Derivative gain (diagonal matrix)

tic  % Start timer
for k = 1:length(t)
    d=d1(x,t(k));%Continuous disturbance
    if k>500 && k<1500
        d=d1(x,t(k))+d2(x,t(k));    
    elseif k>2500 && k<3500
        d=d1(x,t(k))+d3(x,t(k));     
    end
    
      error = x_ref(t(k))-x(1)  ;
      e=x_ref(t(k))-x(1);
      de=dx_ref(t(k))-x(2);
      integral_error = integral_error + error * dt;
      derivative_error = (error - previous_error) / dt;

      u = (Kp * error + Ki * integral_error + Kd * derivative_error);

    
    previous_error = error;
    
    % System dynamic equation
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d+U(x, t(k));
    
    
    % Update state using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;
    
    % Store state trajectory
    x_d(:, k) = x_ref(t(k));
    x_trajectory_PID(:, k) = x;
    control_input_PID(k) = u;
    x_error_PID(k) = e;
    x_derror_PID(k) = de;

end
% Plotting
elapsed_time = toc;  % End timer and return elapsed time
fprintf('Total simulation time is %.4f seconds\n', elapsed_time);


% Plotting
figure;
plot(x_error_PID(1, :),x_derror_PID(1, :), 'r');hold on;
xlabel('e', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;






figure;
plot(t, control_input_PID(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;


% Plotting
figure;
plot(t, x_error_PID(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_trajectory_PID(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;