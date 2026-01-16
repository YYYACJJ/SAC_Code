%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-100, 100]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common pitfalls in SAF function design:
% The function is not designed according to the error range, which causes
% the SAC to fail to track the function.
%
% We hope this code is helpful and inspires new ideas.
% If it is useful to you, please consider citing our paper:
% "Attractive Surface-Based State-Attracted Control for Uncertain Nonlinear Systems"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter Initialization
dt = 0.001;        % Simulation time step
t_end = 5;        % Simulation end time
t = 0:dt:t_end;   % Time vector



% System initial state
x = [0; 0];   % Initial state
x_ref = @(t)  0.05 + (t > 2) *0.02*sin(t*pi);   % Target position
dx_ref =@(t)  (t > 2) *0.02*pi*cos(t*pi);

% Nonlinear function and gain function definitions
f = @(x, t) -2*x(1)-cos(t)*x(2)*x(2)*x(2);  % Nonlinear function f(x,t)
g = @(x) 1+x(1);               % Gain function g(x)
U=@(x,t) 0.5*sin(t)*x(1)+0.6*cos(t)*x(2);



% Result storage
x_trajectory = zeros(2, length(t));  % Record state
control_input = zeros(1, length(t)); % Record control input
s_trajectory = zeros(1, length(t));  % Record sliding surface trajectory
previous_error=0.04;

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

% Main simulation loop

k_cs=320;%
K_f1=0.17;
K_c1=1;

b1=0.040;
b2=0.0401;
[a1,a2,a3,a4]=afcn(b1,b2,K_f1,K_c1);
tic  % Start timer
for k = 1:length(t)

    d=d1(x,t(k));%Continuous disturbance
    if k>500 && k<1500
        d=d1(x,t(k))+d2(x,t(k));    
    elseif k>2500 && k<3500
        d=d1(x,t(k))+d3(x,t(k));     
    end
    
    
    % Construct surface
    e=x_ref(t(k))-x(1);
    de=dx_ref(t(k))-x(2);
    game1=fcn(e,K_f1,K_c1,b1,b2,a1,a2,a3,a4);
    s=de-game1;
    u = k_cs * s;

    
    % System dynamic equation
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d+U(x, t(k));
    
    % Update state using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;
    
    % Store state trajectory and sliding surface trajectory
    x_trajectory_SSAC(:, k) = x;
    control_input_SSAC(k) = u;
    x_error_SSAC(k) = e;
    x_derror_SSAC(k) = de;

end
elapsed_time = toc;  % End timer and return elapsed time
fprintf('Total simulation time is %.4f seconds\n', elapsed_time);

% Plotting
figure;
plot(x_error_SSAC(1, :),x_derror_SSAC(1, :), 'r');hold on;
xlabel('e', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;






% Plotting
figure;
plot(t, control_input_SSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;





% Plotting
figure;
plot(t, x_error_SSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_trajectory_SSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.001, 0.08]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;





%SSAF
function [a1,a2,a3,a4] = afcn(b1,b2,K_lf,K_lc)
    C1=fcnf(b1,K_lf);
    C2=dfcnf(b1,K_lf);
    C3=fcnc(b2,K_lc);
    C4=dfcnc(b2,K_lc);

    h=b1-b2;
    numerator = 2*h*(C1 - C3 + C2*h) + (b1^2 - b2^2 - 2*b1*h)*(C2 - C4);
    denominator = 2*h*(b1^3 - b2^3 - 3*b1^2*h) - 3*h*(b1^2 - b2^2 - 2*b1*h)*(b1 + b2);
    a4 = numerator / denominator;
    % Calculate a3
    a3 = (3 * a4 * h * (b1 + b2) + C2 - C4) / (-2 * h);
    % Calculate a2
    a2 = C2 - 2 * a3 * b1 - 3 * a4 * b1^2;
    % Calculate a1
    a1 = C1 - a2 * b1 - a3 * b1^2 - a4 * b1^3;
end






function game=fcn(e,K_lf,K_lc,b1,b2,a1,a2,a3,a4)

game=fcnf(e,K_lf);
if abs(e) > b1 && abs(e) < b2
game=fcna(e,a1,a2,a3,a4);
elseif abs(e) <= b1
game=fcnf(e,K_lf);
elseif abs(e) >= b2
game=fcnc(e,K_lc);
end
end% 参数初始化


function game=fcnf(e,K_lf)
ke=3000;
game=-K_lf * (atan(ke * abs(e)) + 0.00000001) * atan(10000 * e);
end% 参数初始化

function game=fcnc(e,K_lc)
ke=150;
game=-K_lc*((pi/2-atan(450*abs(e)))*atan(ke*abs(e))+0.00000001)*atan(10000*e);
end% 参数初始化

function dgame = dfcnf(e, K_lf)
ke = 3000;
epsi = 1e-8;

dgame = -K_lf * ( ...
    (ke * sign(e) ./ (1 + (ke * abs(e)).^2)) .* atan(10000 * e) ...
  + (atan(ke * abs(e)) + epsi) .* (10000 ./ (1 + (10000 * e).^2)) ...
);
end

function dgame = dfcnc(e, K_lc)
ke = 150;
epsi = 1e-8;

A = (pi/2 - atan(450*abs(e)));
B = atan(ke*abs(e));
C = atan(10000*e);

% A' = -(450*sign(e)) / (1 + (450*|e|)^2)
dA = -(450 * sign(e)) ./ (1 + (450 * abs(e)).^2);

% B' = (ke*sign(e)) / (1 + (ke*|e|)^2)
dB = (ke * sign(e)) ./ (1 + (ke * abs(e)).^2);

% C' = 10000 / (1 + (10000*e)^2)
dC = 10000 ./ (1 + (10000 * e).^2);

dgame = -K_lc * ( ...
    (dA .* B + A .* dB) .* C + (A .* B + epsi) .* dC ...
);
end


function game = fcna(e,a1,a2,a3,a4)
game = sign(e) .*(a1 ...
  + a2 * abs(e) ...
  + a3 * abs(e).^2 ...
  + a4 * abs(e).^3);
end