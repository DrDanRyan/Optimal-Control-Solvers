%% Setup the problem
close all
clear all
global c m r umax;
T = 5;
nSteps = 1000;
tspan = linspace(0, T, nSteps+1); 
tspanExtra = linspace(T, 2*T, nSteps+1);
x0 = 1;
c = 1.5;
m = 3;
r = .05;
umax = 1;
prob = rk4_test_problem();  % all problem definitions are in the function make_test_problem
nControlPts = 101;


%% Compute equilibrium solution 
xGuess = 2.7;
lamGuess = 2.2;
uGuess = .7;

[~, ~, uStar] = compute_equilibrium(prob, xGuess, lamGuess, uGuess, r);


%% Build RK4InfiniteIntegrator
integrator = RK4InfiniteIntegrator(tspan, tspanExtra, uStar);


%% Solve the problem
soln = single_shooting_new(prob, x0, tspan, nControlPts, ...
                           'Integrator', integrator, ...
                           'u0', uStar);
