%% Setup the problem
global c m umax;
T0 = 0;
TF = 10;
time = linspace(T0, TF, 1001); 
x0 = 1;
c = 1.5;
m = 3;
umax = 1;
prob = make_test_problem();  % all problem definitions are in the function make_test_problem

%% Forward/Backward Sweep with odevr7
sweepSoln = fb_sweep(prob, x0, [T0, TF])

%% Single Shooting 
controlPtArray = linspace(T0, TF, 51);
nlpSoln = single_shooting(prob, x0, controlPtArray)

%% Take nlpSoln as initial guess and solve using bvp5c
options.u0 = nlpSoln.u;
bvpSoln = bvp_solver(prob, x0, [T0, TF], options)

%% Plot all three solution controls 
plot(time, sweepSoln.u(time), time, nlpSoln.u(time), time, bvpSoln.u(time))
legend('f/b sweep', 'single shooting', 'bvp solver')