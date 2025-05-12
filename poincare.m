clear; clc; close all;

% System Parameters
params.m = 1.0;
params.k = 1.0;
params.c = 0.04;
params.alpha = -0.5;
params.D = 1.0;
params.A = 1.0;
params.beta = 0.95;
params.gamma = 0.05;
params.n = 2;
params.F0 = 0.25;         % choose F0 > or < F_CR
params.Omega = 1.0;        % Forcing frequency
params.zeta = 0.02;
params.omega = sqrt(params.k / params.m);
params.epsilon = 1;

% Initial Conditions
IC = [0.01; 0.0; 0.0];      % initial condition

% Time Span: long enough to see steady-state behavior
T_force = 2*pi / params.Omega;
t_transient = 0;               % non-zero for transients to die out
t_final = 10000;
tspan = [0 t_final];

% Integrate with stiff solver
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t, sol] = ode15s(@(t, state) BW_HysteresisFullSystem(t, state, params), tspan, IC, options);

% Extract times at multiples of forcing period
poincare_times = t_transient : T_force : t_final;
poincare_states = deval(ode15s(@(t, state) BW_HysteresisFullSystem(t, state, params), tspan, IC, options), poincare_times);

% Plot Poincaré Section
figure;
plot(poincare_states(1,:), poincare_states(2,:), 'k.', 'MarkerSize', 8);
xlabel('x'); ylabel('y');
title('Poincaré Section (x vs y at t = nT)');
axis equal;
grid on;
