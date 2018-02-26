clear variables;
close all;

syms mu p0 q r a y1 y2 real;

%% solution for the recursive KF
x0s = mu;
% propagated covariance at t1-
p1s = a^2*p0 + q;
k1 = p1s/(r+p1s);
x1s = a*x0s + k1*(y1-a*x0s);
% propagated covariance at t1+
p1sp = (1-k1)*p1s;
% propagated covariance at t2-
p2s = a^2*p1sp + q;
k2 = p2s/(r+p2s);
x2s = a*x1s + k2*(y2-a*x1s);

%% solution for the batch KF with 1 measurement
M = [1/p0+a^2/q -a/q;
    -a/q 1/r+1/q];
b = [mu/p0;
    y1/r];
xsol1 = M\b;

simplify(xsol1(2) - x1s)

%% solution for the batch KF with 2 measurements
M = [1/p0 + a^2/q -a/q 0;
    -a/q 1/r+1/q+a^2/q -a/q;
    0 -a/q 1/r+1/q];
b = [mu/p0;y1/r;y2/r];
xsol2 = M\b;

simplify(xsol2(3) - x2s)