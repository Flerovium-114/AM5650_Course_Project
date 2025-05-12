% Main script to plot hysteresis loops from the paper
% Requires bouc_wen_ode.m - a simplified ode function, adjusted only for
% the initial computation

clear; close; clc;
% Common parameters
params.alpha = 0.4;     % Negative stiffness ratio
params.A = 1.0;          % Hysteresis amplitude control
omega = 1.0;             % Natural frequency
zeta = 0.02;             % Damping ratio
F0 = 0.7;         %0.5       % Excitation amplitude
Omega = 1.0;             % Forcing frequency
tspan = [0 14];        % Simulation time span

% Case 1: Softening hysteresis (Fig 1)
params1 = params;
params1.beta = 0.95;
params1.gamma = 0.05;

% Case 2: Hardening hysteresis (Fig 2)
params2 = params;
params2.beta = 0.35;
params2.gamma = -0.65;

% Solve both systems
[t1, X1] = solve_system(params1, omega, zeta, F0, Omega, tspan);
[t2, X2] = solve_system(params2, omega, zeta, F0, Omega, tspan);

% Plot results

plot_hysteresis(X1, 'Softening Hysteresis (β=0.95, γ=0.05)');
plot_hysteresis(X2, 'Hardening Hysteresis (β=0.35, γ=-0.65)');

function [t,X] = solve_system(params, omega, zeta, F0, Omega, tspan)
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t,X] = ode45(@(t,X) bouc_wen_ode(t,X,params,omega,zeta,F0,Omega),...
                  tspan, [0 0 0], options);
    
    % Remove transient (keep last 50% of data)
    idx = t > tspan(2)/2;
    t = t(idx);
    X = X(idx,:);
end

function plot_hysteresis(X, titleText)
    figure;
    plot(X(:,1), X(:,3), 'b', 'LineWidth', 1.5);
    xlabel('Displacement x');
    ylabel('Hysteretic displacement z');
    title(titleText);
    grid on;
    axis tight;
    set(gca, 'FontSize', 12);
end

