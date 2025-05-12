clear; clc; close all;

% System parameters
params.m = 1.0;
params.k = 1.0;
params.c = 0.04;
params.alpha = -0.5;
params.D = 1.0;
params.A = 1.0;
params.beta = 0.95;
params.gamma = 0.05;
params.n = 2;
params.F0 = 0.25;
params.Omega = 1.0;
params.epsilon = 1.0;
params.zeta = 0.02;
params.omega = 1.0;

% Initial conditions: two nearby points
delta0 = 1e-5; 
state1 = [0.05; 0.0; 0.0]; 
state2 = state1 + [delta0; 0; 0]; 

% Integration time
tspan = linspace(0, 500, 5000);

% Integrate both trajectories
[~, sol1] = ode15s(@(t, state) BW_HysteresisFullSystem(t, state, params), tspan, state1);
[~, sol2] = ode15s(@(t, state) BW_HysteresisFullSystem(t, state, params), tspan, state2);

% Calculate distance growth over time
dist = sqrt(sum((sol1 - sol2).^2, 2));

% Avoid zero values and compute Lyapunov rate
dist(dist <= 0) = delta0;
log_dist = log(dist / delta0);

% Time vector (matched to log_dist)
%t = linspace(tspan(1), tspan(2), length(log_dist));

% Estimate slope (Lyapunov exponent) using linear fit over entire tspan
p = polyfit(tspan, log_dist, 1);
lambda = p(1);

% Plot
figure;
plot(tspan, log_dist, 'c', 'LineWidth', 2); hold on;
plot(tspan, polyval(p, tspan), 'r--', 'LineWidth', 2);
xlabel('Time');
ylabel('log(\delta(t)/\delta_0)');
title(['Lyapunov Exponent Estimation : \lambda = ' num2str(lambda, '%.4f')]);
legend('log(\deltat(t)/\delta_0)', ['Linear fit, slope = ' num2str(lambda, '%.4f')], 'Location', 'northwest');
grid on;


gamma_beta = params.gamma + params.epsilon*params.beta;

% Function to compute F_CR
F_CR = @(zeta, omega, alpha, A, gamma, epsilon, beta, Omega, tau_0) ...
    abs( (4 * zeta * omega^3 * (alpha + (1-alpha)*A)^2) ./ ...
    ((gamma + epsilon*beta) * (1-alpha) * A^2 * 2*pi) ) .* ...
    sinh(Omega * pi ./ (2 * omega * sqrt((alpha + (1-alpha)*A)/2)));

% Using existing parameters
tau_0 = pi/(2*params.omega);  % Assuming tau_0 = π/(2ω) based on common formulations

% Calculate F_CR
Fcr_value = F_CR(params.zeta, params.omega, params.alpha, params.A, ...
                params.gamma, params.epsilon, params.beta, params.Omega, tau_0);

% Display the result
disp(['F_CR = ', num2str(Fcr_value)]);

coeff1 = 0.5 * params.omega^2 * (params.alpha + (1-params.alpha)*params.A);
coeff2 = (1/12) * params.omega^2 * params.A^2 * (1-params.alpha)*gamma_beta;

x_equil = [0; sqrt(coeff1/(4*coeff2)); -sqrt(coeff1/(4*coeff2))];
disp('Equilibrium points (x_equil):');
disp(x_equil);
