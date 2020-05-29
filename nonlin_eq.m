function dx = nonlin_eq(x,u,sys)
M = sys.M;
J = sys.J;
mgl = sys.mgl;
B = sys.B;
K = sys.K;

dx = zeros(4, 1);
dx(1) = x(2);
dx(2) = 1/M * (K*(x(3) - x(1)) - mgl * sin(x(1)));
dx(3) = x(4);
dx(4) = 1/J * (-B*x(3) - K*(x(3) - x(1)) + u);
end