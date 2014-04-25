%% Setup the problem
clear all
global c m umax;
T0 = 0;
TF = 10;
tspan = linspace(T0, TF, 101); 
x0 = [1; 0];
c = 1.5;
m = 3;
umax = 1;
prob = rk4_test_problem();  % all problem definitions are in the function make_test_problem


%% Solve the problem
soln = rk4_single_shooting2(prob, x0, tspan)
