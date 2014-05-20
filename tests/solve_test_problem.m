%% Setup the problem
close all
clear all

T = 5;
nSteps = 1000;
tspan = linspace(0, T, nSteps+1); 
tspanExtra = linspace(T, 2*T, nSteps+1);

x0 = 1;

p.c = 1.5;
p.m = 3;
p.r = .05;
ControlBounds = [0, 1];
prob = TestOCProblem(p, ControlBounds);

nControlPts = 101;


%% Compute equilibrium solution 
xGuess = 2.7;
lamGuess = 2.2;
uGuess = .7;

[~, ~, uStar] = compute_equilibrium(prob, xGuess, lamGuess, uGuess, p.r);


%% Build RK4InfiniteIntegrator
integrator = RK4InfiniteIntegrator(tspan, tspanExtra, uStar);


%% Solve the problem
soln = single_shooting_new(prob, x0, tspan, nControlPts, ...
                           'Integrator', integrator, ...
                           'u0', uStar);
