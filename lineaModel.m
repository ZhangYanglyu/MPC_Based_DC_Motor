function [Ac, Bc, Kc]  = lineaModel(dx, u, x, sys)
M = sys.M;
J = sys.J;
mgl = sys.mgl;
B = sys.B;
K = sys.K;

x1 = x(1);

Ac = [0 1 0 0;
    -(K + mgl*cos(x1))/M 0 K/M 0;
      0 0 0 1;
      K/J 0 -(B+K)/J 0];
Bc = [0; 0; 0; 1/J];

Kc = dx-Ac*x-Bc*u;

end