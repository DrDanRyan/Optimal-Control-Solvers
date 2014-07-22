%% Setup the problem
close all
clear all

T = 10;
nSteps = 500;
tspan = linspace(0, T, nSteps+1); 

x0 = 1;
p.c = 1.5;
p.m = 3;
p.r = .05;
ControlBounds = [0, 1];
prob = TestOCProblem(p, ControlBounds);

integrator = RK4Integrator(tspan);

nControlPts = 21;
v = .5 + .5*rand(nControlPts, 1);
control = PWLinearControl(integrator.t, nControlPts, 1);

epsilon = 1e-4; % finite difference step-size

%% Compute sensitivities with backprop
u = control.compute_u(v);
[~, J] = integrator.compute_states(prob, x0, u);
[~, dJdu] = integrator.compute_adjoints(prob, u);    
dJdv = control.compute_dJdv(dJdu);


%% Compute finite difference sensitivities
dJdv_FD = zeros(nControlPts, 1);
for i=1:nControlPts
   v(i) = v(i) + epsilon;
   u = control.compute_u(v);
   [~, J_FD] = integrator.compute_states(prob, x0, u);
   dJdv_FD(i) = (J_FD - J)/epsilon;
   v(i) = v(i) - epsilon;
end

%% Compare Derivatives
fprintf('Max Norm: %s\n', max(dJdv_FD - dJdv));
fprintf('Mean AbsDev: %s\n', sum(abs(dJdv_FD - dJdv))/nControlPts);



