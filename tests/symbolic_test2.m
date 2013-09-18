x = sym('x', [1 1]);
u = sym('u', [1 1]);
params.m = .5;
params.c = 4;

syms t m c
obj(t, x, u) = x^2 + c*u^2;
rhs(t, x, u) = x*(m-x) - u;
bounds = [0, 1];

prob = make_from_symbolic(obj, rhs, 1, 1, params, bounds);
% soln = bvp_solver(prob, .1, [0, 5])