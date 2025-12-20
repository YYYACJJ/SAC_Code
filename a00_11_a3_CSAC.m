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
%Torque limited to [-500, 500]

% 参数初始化
dt = 0.001;        % 仿真时间步长
t_end = 5;        % 仿真结束时间
t = 0:dt:t_end;   % 时间向量



% 系统初始状态
x = [0; 0];   % 初始状态 [x1, x2, x3, x4]
x_ref = @(t)  0.05 + (t > 2) *0.02*sin(t*pi);   % 目标位置
dx_ref =@(t)  (t > 2) *0.02*pi*cos(t*pi);

% 非线性函数和增益函数定义
f = @(x, t) (t <= 1) *(-2*x(1)-cos(t)*x(2)*x(2)*x(2))+(t > 1) *200*(cos(t)*cos(t)*x(1)+2*sin(t)*cos(t));  % 非线性函数 f(x,t)
g = @(x) 1+x(1);               % 增益函数 g(x)
U=@(x,t) 0.5*sin(t)*x(1)+0.6*cos(t)*x(2);
integral_error = 0;


% 结果存储
x_trajectory = zeros(2, length(t));  % 记录状态
control_input = zeros(1, length(t)); % 记录控制输入
s_trajectory = zeros(1, length(t));  % 记录滑模面轨迹
previous_error=0.04;



% ====== 2) 预生成每个 1s 区间的随机变量（均匀分布在 [-1,1]）======
N  = floor(t_end) + 2;       % 区间数（多留一个防越界）
U1 = 2*rand(N,1) - 1;           % 随机项
U2 = 2*rand(N,1) - 1;           % 随机项
U3 = 2*rand(N,1) - 1;           % 随机项


% ====== 3) 一个小工具：给定 t，返回“当前 1s 区间”的随机值（分段常数）======
stepRand = @(t,U) U( min(floor(t)+1, numel(U)) );

d4=@(x,t)-10*sin(t)*cos(t)+10*stepRand(t,U1)*x(2)+5*stepRand(t,U2)*x(1)+2*stepRand(t,U3);%持续扰动
% 主仿真循环

k_cc=500;
k_c1=250;

for k = 1:length(t)

    d=d4(x,t(k));%持续扰动

    
    
    
    
    % 构造滑模面
    e=x_ref(t(k))-x(1);
    de=dx_ref(t(k))-x(2);
    integral_error = integral_error + e * dt;
    derivative_error = (e - previous_error) / dt;
    previous_error = e;
    

    game1=fcn(e,k_c1);
    s=de-game1;
    u = k_cc * s;

    
    % 系统动态方程
    dx1 = x(2);
    dx2 = f(x, t(k)) + g(x) * u + d+U(x, t(k));
    
    % 使用Euler方法更新状态
    x(1) = x(1) + dx1 * dt;
    x(2) = x(2) + dx2 * dt;
    
    % 存储状态轨迹和滑模面轨迹
    x_trajectory_CSAC(:, k) = x;
    control_input_CSAC(k) = u;
    x_error_CSAC(k) = e;
    x_derror_CSAC(k) = de;
    error_s1_CSAC(k)=derivative_error-game1;
end



% 绘图
figure;
plot(t,error_s1_CSAC(1, :), 'r');hold on;
xlabel('t', 'FontSize', 16);
ylabel('error_s1', 'FontSize', 16);
ylim([-0.1, 0.1]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;






% 绘图
figure;
plot(t, control_input_CSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('u', 'FontSize', 16);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

% 绘图
figure;
plot(t, x_error_CSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Error', 'FontSize', 16);
ylim([-0.01, 0.01]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

figure;
plot(t, x_trajectory_CSAC(1, :), 'r');
xlabel('Time(s)', 'FontSize', 16);
ylabel('Position', 'FontSize', 16);
ylim([-0.001, 0.08]);
grid on;
set(gca, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 16;
hold on;

function game=fcn(e,k_L)


    game=-k_L*((pi/2-atan(450*abs(e)))*atan(1*abs(e))+0.00000001)*atan(10000*e);

end

