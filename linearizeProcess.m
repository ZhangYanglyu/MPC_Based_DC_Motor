syms x1 x2 x3 x4 M K J B mgl u A B

f1 = x2;
f2 = 1/M * (K*(x3-x1)-mgl*sin(x1));
f3 = x4;
f4 = 1/J * (-B*x3 - K*(x3-x1) + u);

A(1,1) = diff(f1, x1);
A(1,2) = diff(f1, x2);
A(1,3) = diff(f1, x3);
A(1,4) = diff(f1, x4);
B(1,1) = diff(f1, u);

A(2,1) = diff(f2, x1);
A(2,2) = diff(f2, x2);
A(2,3) = diff(f2, x3);
A(2,4) = diff(f2, x4);
B(2,1) = diff(f2, u);

A(3,1) = diff(f3, x1);
A(3,2) = diff(f3, x2);
A(3,3) = diff(f3, x3);
A(3,4) = diff(f3, x4);
B(3,1) = diff(f3, u);

A(4,1) = diff(f4, x1);
A(4,2) = diff(f4, x2);
A(4,3) = diff(f4, x3);
A(4,4) = diff(f4, x4);
B(4,1) = diff(f4, u);
disp('A = ');
disp(A);

disp('B = ');
disp(B);