% Maintainer: Yang Zhang
% Time: 30/5/2020

% dc motor system
% J \ddot{\theta} + b \dot{\theta} = K_t i
% L di/dt + R i = V - K_e \dot{\theta}

% system parameters
% J: modent of inertia of the rotor
% b: viscous friction constant
% K_e: force constant
% K_t: torque constant
% R: resistant
% L: inductance
J = 0.01;
b = 0.1;
K = 0.01;
R = 1;
L = 0.5;

% choose state variables x = [x_1(t) x_2 (t)]^T
% x_1 = \dot{\theta}
% x_2 = i
% which can rewrite as 
% \dot{x} = A_c x + b_c u
% y = C x

A_c = [-b/J K/J; -K/L -R/L];
B_c = [0; 1/L];
C_c = [1 0];
D_c = zeros(1,1);

% samping interval
Ts = 0.01;

% continuous-time plant model is discretized
[Ap, Bp, Cp, Dp] = c2dm(A_c, B_c, C_c, D_c, Ts);

% prediction horizon
% control horizon

Np = 100;
Nc = 100;

% agumented state-space of plant model
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np); 

% n£ºthe number of state in agument state-space
% n_in£ºthe number of input in agument state-space
[n, n_in] = size(B_e);

% state value of plant model 
% state value of agumented state-space of plant model
xm = [0;0];
Xf = zeros(n,1);

% simulation time£º1000 samping interval
N_sim = 1000;

% r£ºset-point Speed
% u£ºinput in plant model
% y£ºoutput in plant model
r = ones(N_sim, 1);
u = 0; % u(k-1) = 0
y = 0;

for kk = 1:N_sim
    DeltaU = inv(Phi_Phi + 0.3*eye(Nc, Nc))*(Phi_R*r(kk) - Phi_F*Xf);
    deltau = DeltaU(1,1);
    u = u + deltau;
    % store next input added to plant model
    u1(kk) = u;
    % store current output of plant model
    y1(kk) = y;
    
    % store current state of plant model
    xm_old = xm;        
    % calculate next state of plant model
    xm = Ap*xm + Bp*u; 
    y = Cp*xm;
    % calculate next state of agumented state-space of plant model
    Xf = [xm-xm_old; y];    
end

k = 0 : (N_sim-1);
figure(1)
plot(k, r, '--k', k, y1', 'b','linewidth',2);
set(gca,'FontName','Times New Roman','FontSize',14);
xlabel('Sampling Instant');
ylabel('Speed (rad/s)');
% set(gca,'FontName','Times New Roman','FontSize',14);
set(gcf, 'unit', 'centimeters', 'position', [10 5 16 6])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
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
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
h=legend('$u$');
set(h,'Interpreter','latex'); grid;