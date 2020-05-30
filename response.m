% dc motor system

% system parameters
% m = 3;
% k = 1;
% c = 0.5;
J = 0.01;
b = 0.1;
K = 0.01;
R = 1;
L = 0.5;
% choose state variables x = [x_1(t) x_2 (t)]^T
% which can rewrite as 
% \dot{x} = A_c x + b_c u
% y = C x

A_c = [-b/J K/J; -K/L -R/L];
B_c = [0; 1/L];
C_c = [1 0];
D_c = zeros(1,1);

% samping interval
Ts = 0.01;

% continuous-time model is discretized
[Ap, Bp, Cp, Dp] = c2dm(A_c, B_c, C_c, D_c, Ts);

% % 离散的状态空间模型
% Ap = [1 1; 0 1];
% Bp = [0.5; 1];
% Cp = [1 0];
% Dp = 0;

Np = 100;
Nc = 100;

% agumented 状态空间模型
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np); 

% n：状态变量的个数
% n_in：输入变量的个数
[n, n_in] = size(B_e);

% 初始条件
xm = [0;0];
Xf = zeros(n,1);

% 仿真时间：1000个采样时间间隔
N_sim = 1000;

% r：目标设定值
% u：系统的输入
% y：系统的输出
r = ones(N_sim, 1);
u = 0; % u(k-1) = 0
y = 0;

for kk = 1:N_sim
    DeltaU = inv(Phi_Phi + 0.3*eye(Nc, Nc))*(Phi_R*r(kk) - Phi_F*Xf);
    deltau = DeltaU(1,1);
    u = u + deltau;
    % 保存系统的输入变量
    u1(kk) = u;
    % 保存系统的输出变量
    y1(kk) = y;
    
    xm_old = xm;        % 当前时刻的采样值
    xm = Ap*xm + Bp*u;  % 下一个采样时刻估计值 = 当前时刻的采样值 + 当前时刻的输入
    y = Cp*xm;
    Xf = [xm-xm_old; y];    
end

% 仿真采样时刻
k = 0 : (N_sim-1);
figure(1)
plot(k, r, '--k', k, y1', 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant');
ylabel('Speed (rad/s)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
axis([0 N_sim-1  0 1.2]);
set(gca, 'YTick', [0:0.6:1.2]);
h=legend('$\theta_r$','$\theta$');
set(h,'Interpreter','latex'); grid;

figure(2)
plot(k, u1, 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant');
ylabel('Voltage (v)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); % 加入虚线
h=legend('$u$');
set(h,'Interpreter','latex'); grid;