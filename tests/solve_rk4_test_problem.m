%% Setup the problem
clear all
global c m umax;
T0 = 0;
TF = 10;
tspan = linspace(T0, TF, 1001); 
x0 = 1;
c = 1.5;
m = 3;
umax = 1;
prob = rk4_test_problem();  % all problem definitions are in the function make_test_problem

nControlPts = 101;

%% Solve the problem
soln = single_shooting_new(prob, x0, tspan, nControlPts, 'IntegratorType', 'RK4Infinite')
