%% Setup the problem
close all
clear all

T = 10;
nSteps = 500;
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

lb = [0; -Inf; ControlBounds(1)];
ub = [Inf; Inf; ControlBounds(2)];
[~, ~, uStar] = ...
   compute_equilibrium(prob, xGuess, lamGuess, uGuess, lb, ub, p.r);


%% Build RK4InfiniteIntegrator
integrator = RK4InfiniteIntegrator(tspan, tspanExtra, uStar);


%% Solve the problem
soln = single_shooting(prob, x0, tspan, nControlPts, ...
                       'u0', uStar');%, ...
                       %Integrator', integrator);
