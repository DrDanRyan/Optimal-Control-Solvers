%% Setup the problem
clear all
global c m delta umax;
T = 10;
nSteps = 1000;
tspan = linspace(0, T, nSteps+1); 
tspanExtra = tspan + T;
x0 = 1;
c = 1.5;
m = 3;
delta = .10;
umax = 1;
prob = rk4_test_problem();  % all problem definitions are in the function make_test_problem
nControlPts = 101;


%% Compute equilibrium solution 
xGuess = 

%% Build RK4InfiniteIntegrator
integrator = RK4InfiniteIntegrator(tspan, tspanExtra, 4);


%% Solve the problem
soln = single_shooting_new(prob, x0, tspan, nControlPts, ...
                           'Integrator', integrator)
