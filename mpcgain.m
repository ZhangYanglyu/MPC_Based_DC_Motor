function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np)

% m1: �����������
% n1: ״̬��������
% n_in: �����������
[m1, ~] = size(Cp);
[n1, n_in] = size(Bp);

% ��ɢ��ʽ��״̬�ռ�ģ��---->extend state-model
A_e = eye(n1+m1,n1+m1); 
A_e(1:n1,1:n1) = Ap; 
A_e(n1+1:n1+m1,1:n1) = Cp*Ap; 
B_e = zeros(n1+m1,n_in); 
B_e(1:n1,:) = Bp; 
B_e(n1+1:n1+m1,:) = Cp*Bp; 
C_e = zeros(m1,n1+m1); 
C_e(:,n1+1:n1+m1) = eye(m1,m1);

% F ����
% n = n1 + m1;
h(1, :) = C_e;
F(1, :) = C_e * A_e;
for kk = 2 : Np
    h(kk, :) = h(kk-1, :)*A_e;
    F(kk, :) = F(kk-1, :)*A_e;
end

% Phi ����
v = h*B_e;
Phi = zeros(Np, Nc);
Phi(:, 1) = v;   % ��һ��
for i = 2:Nc
    Phi(:, i) = [zeros(i-1, 1); v(1:Np-i+1, 1)];
end

BarRs = ones(Np, 1); % ���������ڲ���
Phi_Phi =   Phi' * Phi;
Phi_F = Phi' * F;
Phi_R = Phi' * BarRs;

end

% �������
% Ap = 0.8;
% Bp = 0.1;
% Cp = 1;
% Nc = 4;
% Np = 10;
% 
% [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np)
% 
% ��ʼ����
% x_ki = [0.1, 0.2]';   % ϵͳ�ĳ�ʼֵ
% r_ki = 1;             % �����ο�ֵ 
% r_w = 10;             % ���ڲ���
% BarR = r_w .* eye(Nc);
% R_s = ones(Np, 1).*r_ki;
% delta_u = inv(Phi_Phi + BarR)*(Phi_R .* r_ki - Phi_F*x_ki); 