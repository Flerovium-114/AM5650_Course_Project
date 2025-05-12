clear; clc; close all;

% System parameters
params.m = 1.0;
params.k = 1.0;
params.c = 0.04;
params.alpha = -0.5;
params.D = 1.0;
params.A = 1;
params.beta = 0.95;
params.gamma = 0.05;
params.n = 2;
params.F0 = 0.25;
params.Omega = 1.0;
params.epsilon = 1.0;
params.zeta = 0.02;
params.omega = 1.0;

% Initial conditions
state0 = [0.5, -0.5, 0.5];

% Integration
tspan = [0, 1000]; % Long time for attractor to emerge
[t, sol] = ode15s(@(t, state) BW_HysteresisFullSystem(t, state, params), tspan, state0);

% Extract x, xdot, z
x = sol(:,1);
xdot = sol(:,2);
z = sol(:,3);

% Plot 3D phase space
figure;
scatter3(state0(1), state0(2), state0(3), 80, 'r', 'filled'); 
hold on;
plot3(x, xdot, z, 'k');
xlabel('x (Displacement)', 'FontSize', 24);
ylabel('v (Velocity)', 'FontSize', 24);
zlabel('z (Hysteresis State)', 'FontSize', 24);
title('3D Phase Space Trajectory for F_{0} = 0.25', 'FontSize', 26);
set(gca, 'FontSize', 24);
grid on;
view(45, 30); % good angle to spot attractor shapes

% % Find local maxima of z(t)
% [zn_peaks, peak_locs] = findpeaks(z, t);

% % Create return map: z_{n+1} vs z_n
% zn = zn_peaks(1:end-1);
% zn_plus1 = zn_peaks(2:end);
% 
% % Plot return map
% figure;
% plot(zn, zn_plus1, 'ko', 'MarkerFaceColor', 'k');
% hold on;
% plot([min(zn) max(zn)], [min(zn) max(zn)], 'r--'); % line z_{n+1} = z_n
% xlabel('$z_n$', 'Interpreter', 'latex', 'FontSize', 16);
% ylabel('$z_{n+1}$', 'Interpreter', 'latex', 'FontSize', 16);
% title('Return Map: $z_{n+1}$ vs $z_n$', 'Interpreter', 'latex', 'FontSize', 18);
% grid on;
% axis equal;

% % Compute square root term
% sqrt_term = sqrt(params.A / (params.gamma + params.epsilon * params.beta));
% 
% % Compute x component
% x_val = (1 - params.alpha) / params.alpha * sqrt_term;
% 
% % Fixed points
% fp1 = [-x_val, 0,  sqrt_term];
% fp2 = [ x_val, 0, -sqrt_term];
% 
% % Display results
% disp('Fixed point 1:');
% disp(fp1);
% disp('Fixed point 2:');
% disp(fp2);



