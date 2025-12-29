%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Torque limited to [-100, 100]

% Parameter Initialization
dt = 0.001;       % Simulation time step
t_end = 5;      % Simulation end time
t = 0:dt:t_end;  % Time vector
n=0;
Da=2;
Dt=0.5;


% RBF 
c_a = 50;
HidLayNuma = 5;
OutLayNum = 1;
gamma_a = 300;
eta_a = 70;
c_i_a = 8*[-2 -1  0  1  2
       -2 -1  0  1  2];
b_i_a = 0.4;
hat_W1_a = 0.2*ones(HidLayNuma, OutLayNum);


x=[0;0];



xr_a = @(t) 0.05 + (t > 2) *0.02*sin(t*pi); 
dxr_a =@(t) (t > 2) *0.02*pi*cos(t*pi);


% Nonlinear function and gain function definitions
f = @(x, t) -2*x(1)-cos(t)*x(2)*x(2)*x(2);  % Nonlinear function f(x,t)
g = @(x) 1+x(1);               % Gain function g(x)
U=@(x,t) 0.5*sin(t)*x(1)+0.6*cos(t)*x(2);


% ====== 2) Pre-generate random variables for each 1s interval (uniformly distributed in [-1,1]) ======
N  = floor(t_end) + 2;       % Number of intervals (extra one for safety)
U1 = 2*rand(N,1) - 1;           % Random term for d1
U2 = 2*rand(N,1) - 1;           % Random term for d2
U3 = 2*rand(N,1) - 1;           % Random term for d1
U4 = 2*rand(N,1) - 1;           % Random term for d1
U5 = 2*rand(N,1) - 1;           % Random term for d2
U6 = 2*rand(N,1) - 1;           % Random term for d1
U7 = 2*rand(N,1) - 1;           % Random term for d1
% ====== 3) A small utility: given t, returns the random value for the "current 1s interval" (piecewise constant) ======
stepRand = @(t,U) U( min(floor(t)+1, numel(U)) );


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
tic  % Start timer
for k = 1:length(t)
%%%%%%%%%%%%Disturbance%%%%%%%%%%%%%%%%%    
    d=d1(x,t(k));%Continuous disturbance
    if k>500 && k<1500
        d=d1(x,t(k))+d2(x,t(k));    
    elseif k>2500 && k<3500
        d=d1(x,t(k))+d3(x,t(k));     
    end

    

    e = xr_a(t(k))-x(1);  
    de = dxr_a(t(k))-x(2);  
    s_a = c_a*e + de;
    for j = 1:HidLayNuma
        h1a(j,:) = exp( -norm([x(1); x(2)]-c_i_a(:,j))^2/(2*b_i_a^2) );
    end
    hat_fx1_a = hat_W1_a' * h1a;
    hat_W1_a = hat_W1_a + dt * (gamma_a * s_a(1) * h1a);
    u =   eta_a * sign(s_a) + hat_fx1_a;
%     % System dynamic equation

    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d+U(x, t(k));
    
    % Update state using Euler method
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;



    x_trajectory_NNSMC(:, k) = x;
    control_input_NNSMC(k) = u;
    x_error_NNSMC(k) = e;
    x_derror_NNSMC(k) = de;

    
end
elapsed_time = toc;  % End timer and return elapsed time
fprintf('Total simulation time is %.4f seconds\n', elapsed_time);

figure;
plot(x_error_NNSMC(1, :),x_derror_NNSMC(1, :), 'r');hold on;
xlabel('e', 'FontSize', 16);
ylabel('de', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

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



% Plotting
figure;
plot(t, x_error_NNSMC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_trajectory_NNSMC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.001, 0.08]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;


function S = skew(v)
% Input v is a 3×1 vector, output S is a 3×3 skew-symmetric matrix
    S = [   0    -v(3)   v(2);
           v(3)   0    -v(1);
          -v(2)  v(1)    0   ];
end