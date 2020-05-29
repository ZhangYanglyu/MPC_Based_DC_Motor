% the single-link flexible-joint robots
% M \ddot{q} + mgl sin(q) = K (\theta - q) + \omegat_1
% J \ddot{\theta} + B \dot{\theta} + K (\theta - q) = \tau + \omega_2

% M : inertia of link          1 kg \cdot m^2
% J : joint flexibility        1 kg \cdot m^2
% mgl : none                   10 kg \cdot m^2 / s^2
% B : joint stiffness          0.9 Nm \cdot s/rad
% K : joint stiffness          100 Nm/rad
% g : gravitational contant    9.8 m/s^2

% choose state variables x = [q, \dot{q}, \theta, theta]^T
% the above equation can be rewriten as dx/dt = f(x,u)
% \dot{x}_1 = x_2
% \dot{x}_2 = 1/M (K(x_3 - x_1) - mglsin(x_1) + w_1)
% \dot{x}_3 = x_4
% \dot{x}_4 = 1/J(-Bx - K(x_3 - x_1) + \tau + w_2)

% parameters associated with dynamic model of the system
sys.M = 1;
sys.J = 1;
sys.mgl = 10;
sys.B = 0.9;
sys.K = 100;

% system dynamic model and in the form of dx/dt =f(x,u). 
model = @(x,u) nonlin_eq(x,u,sys);

% samping interval (0.01s)
Ts = 0.01; 

% horizon length
Np = 100; 
Nc = 100;

% total simulation time (10s)
Tfinal = 10;
N_sim = Tfinal/Ts;

% r ：set-point value
% u ：system input
% y ：system output
% xi：system state
% dxm ：the derivitive of system state
xm = [0;0;0;0];
dxm = [0;0;0;0];
u = 0;
r = ones(N_sim, 1);
y = 0;
x1 = 0;
x2 = 0;
x3 = 0;
x4 = 0;

% linearization of nonliear model 
% continuous-time model is discretized

% [A_c, B_c]  = lineaModel(xm, sys);
[A_c, B_c, K_c] = lineaModel(dxm,u,xm,sys);     % 得到连续的线性模型
[Ap, Bp, Kd] = discretize(A_c, B_c, K_c, Ts);   % 将连续的线性模型转化为离散模型
Cp = [1 0 0 0];
Dp = zeros(1,1);

% agumented state-space model
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np); 

% n：number of state variables
% n_in：number of input variables
[n, n_in] = size(B_e);

% Feedback gain
Xf = zeros(n,1);

for kk = 1:N_sim
    DeltaU = inv(Phi_Phi + 0.5*eye(Nc, Nc))*(Phi_R*r(kk) - Phi_F*Xf);
    deltau = DeltaU(1,1);
    u = u + deltau;
    % 保存系统的输入变量
    u1(kk) = u;
    % 保存系统的输出变量
    y1(kk) = y;
    xm_old = xm;                   
    
    [xm, dxm] = RK4(xm_old,u,Ts,model);         % 当前时刻的采样值
    [A_c, B_c, K_c] = lineaModel(dxm,u,xm,sys); % 
    [Ad, Bd, Kd] = discretize(A_c, B_c, K_c, Ts);
    
    xm = Ad * xm + Bd * u;                      % 下一个采样时刻估计值 = 当前时刻的采样值 + 当前时刻的输入
    
%     xm = Ap*xm + Bp*u;              
%     [A_c, B_c]  = lineaModel(xm, sys);
    % continuous-time model is discretized
    % [Ap, Bp, Cp, Dp] = c2dm(A_c, B_c, C_c, D_c, Ts);
    % agumented 状态空间模型
    [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ad, Bd, Cp, Nc, Np); 
    x1(kk) = xm(1);
    x2(kk) = xm(2);
    x3(kk) = xm(3);
    x4(kk) = xm(4);
    y = Cp*xm;
    Xf = [xm-xm_old; y];    
end

% 仿真采样时刻
k = 0 : (N_sim-1);
figure(1)
plot(k, r, '--k', k, y1', 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant (ms)');
ylabel('Position (rad)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
axis([0 N_sim-1  0 1.2]);
set(gca, 'YTick', [0:0.6:1.2]);
h=legend('$q_r$','$q$');
set(h,'Interpreter','latex'); grid;

figure(2)
plot(k, x2, 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant (ms)');
ylabel('Velocity (rad/s)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
axis([0 N_sim-1  -0.7 1.4]);
set(gca, 'YTick', [-0.7:0.7:1.4]);
h=legend('$\dot{q}$');
set(h,'Interpreter','latex'); grid;

figure(3)
plot(k, x3, 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant (ms)');
ylabel('Position (rad)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
axis([0 N_sim-1  0 1.4]);
set(gca, 'YTick', [0:0.7:1.4]);
h=legend('$\theta$');
set(h,'Interpreter','latex'); grid;

figure(4)
plot(k, x4, 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant (ms)');
ylabel('Velocity (rad/s)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
axis([0 N_sim-1  -0.7 1.4]);
set(gca, 'YTick', [-0.7:0.7:1.4]);
h=legend('$\dot{\theta}$');
set(h,'Interpreter','latex'); grid;

figure(5)
plot(k, u1, 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant (ms)');
ylabel('torque (N \cdot m)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
axis([0 N_sim-1  -5 15]);
set(gca, 'YTick', [-5:5:15]);
h=legend('$\tau$');
set(h,'Interpreter','latex'); grid;
