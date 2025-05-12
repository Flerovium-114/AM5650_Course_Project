clear; clc; close all;

% System parameters from paper
params.alpha = -0.5;     % Negative stiffness ratio
params.beta = 0.95;
params.gamma = 0.05;
params.A = 1.0;
omega = 1.0;             % Natural frequency
zeta = 0.02;             % Damping ratio
F0 = 0.5;                % Excitation amplitude (unused for these figures)
Omega = 1.0;             % Forcing frequency (unused for these figures)
epsilon = 1;             % sgn(dx)*sgn(z) = 1

%% Figure 3b: Potential Energy Curve
gamma_beta_term = params.gamma + epsilon*params.beta;
coeff1 = 0.5 * omega^2 * (params.alpha + (1-params.alpha)*params.A);
coeff2 = (1/12) * omega^2 * params.A^2 * (1-params.alpha)*gamma_beta_term;

x = linspace(-2, 2, 500);
V = coeff1 * x.^2 - coeff2 * x.^4;

figure;
plot(x, V, 'b', 'LineWidth', 2);
xlabel('Displacement x');
ylabel('Potential Energy V_{ap}(x)');
title('(b) Potential Energy Curve');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

%% Figure 3a: Phase Space with Heteroclinic Orbits
% Duffing equation approximation
duffing_eq = @(t,X) [X(2); 
                    -coeff1*X(1) + 4*coeff2*X(1)^3];

% Find equilibrium points
x_equil = [0; 
           sqrt(coeff1/(4*coeff2)); 
           -sqrt(coeff1/(4*coeff2))];

% Simulation parameters
tspan = [0 100];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

% Initial conditions near saddle points
ics = [x_equil(2)*1.01 0;
       x_equil(3)*0.99 0;
       0 0.1];

figure; 
hold on;
colors = ['r', 'g', 'b'];

% Plot trajectories
for i = 1:size(ics,1)
    [~,X] = ode15s(duffing_eq, tspan, ics(i,:), options);
    plot(X(:,1), X(:,2), colors(i), 'LineWidth', 1.5);
    
end

% Plot equilibrium points
plot(x_equil, [0;0;0], 'ko', 'MarkerSize', 8,...
    'MarkerFaceColor','k', 'LineWidth', 1.5);

% Formatting
xlabel('Displacement x');
ylabel('Velocity y');
title('Phase Space with Heteroclinic Orbits');
grid on;
axis([-2 2 -2 2]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
hold off;

%% Figure 4: Heteroclinic Orbit
tau = linspace(-10, 10, 1000);
x_het = sqrt(3*(params.alpha + (1-params.alpha)*params.A)/...
          ((1-params.alpha)*gamma_beta_term*params.A^2))...
        * tanh(omega*sqrt((params.alpha + (1-params.alpha)*params.A)/2)*tau);

y_het = sqrt(3*(params.alpha + (1-params.alpha)*params.A)/...
          ((1-params.alpha)*gamma_beta_term*params.A^2))...
        * omega*sqrt((params.alpha + (1-params.alpha)*params.A)/2)...
        .* sech(omega*sqrt((params.alpha + (1-params.alpha)*params.A)/2)*tau).^2;

figure;
plot(x_het, y_het, 'b', 'LineWidth', 2);
hold on;
plot(-x_het, -y_het, 'r', 'LineWidth', 2);
xlabel('Displacement x');
ylabel('Velocity y');
title('Heteroclinic Orbit of Monostable Potential');
grid on;
axis([-3 3 -3 3]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
hold off;
