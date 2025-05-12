clear; clc; close all;

% System parameters (same as in Fig. 6 of the paper)
params.m = 1;
params.k = 1;
params.c = 0.02;
params.alpha = -0.4;
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

% Calculate F_CR (critical force threshold)
gamma_beta = params.gamma + params.epsilon*params.beta;
alpha_term = params.alpha + (1 - params.alpha)*params.A;
F_CR = abs((4 * params.zeta * params.omega^3 * alpha_term^2) / ...
    (gamma_beta * (1 - params.alpha) * params.A^2 * params.Omega * pi) * ...
    sinh(params.Omega * pi / (2 * params.omega * sqrt(alpha_term/2))));
fprintf('Critical force F_CR = %.4f (current F0 = %.4f)\n', F_CR, params.F0);

% Higher resolution grid (adjust based on computational resources)
n = 200; % Number of points in each dimension
x_vals = linspace(-2, 2, n);
v_vals = linspace(-3, 3, n);
basin_image = zeros(n, n);

% Integration parameters
tspan = [0, 15]; % Longer integration for better classification
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.1);

% Perform the basin computation
fprintf('Computing basin of attraction (%dx%d grid)...\n', n, n);
tic; % Start timer

% Process row by row with progress updates
for i = 1:n
    if mod(i, 10) == 0
        fprintf('Row %d of %d (%.1f%%) - elapsed time: %.1f seconds\n', ...
                i, n, 100*i/n, toc);
    end
    
    for j = 1:n
        x0 = x_vals(j);
        v0 = v_vals(i);
        
        % Initial condition
        init = [x0; v0; 0];
        
        % Simulate system
        [t, state] = ode15s(@(t, state) BW_HysteresisFullSystem(t, state, params), ...
                            tspan, init, options);
                        
        % Extract end state
        x_final = state(end, 1);
        
        % Calculate escape/convergence time for coloring
        if abs(x_final) < 1  % Bounded trajectory
            % Find when trajectory stabilizes
            for k = length(t):-1:2
                if abs(state(k,1) - state(k-1,1)) > 0.01
                    break;
                end
            end
            
            % Color based on stabilization time (normalized to [1-64])
            stability_time = t(k)/tspan(2);
            basin_image(i,j) = 1 + floor(min(63, 64*stability_time));
        else
            % Unbounded trajectory (different color scheme)
            basin_image(i,j) = 65 + floor(min(63, 64*abs(x_final)/5));
        end
    end
end

fprintf('Calculation complete (%.1f seconds)\n', toc);

% Create custom colormap to highlight fractal boundaries
% Blue shades for bounded orbits, red shades for unbounded
% blue_cmap = flipud(hot(64));  % Reversed hot colormap for bounded regions
% red_cmap = jet(64);           % Jet colormap for unbounded regions
% combined_cmap = [blue_cmap; red_cmap]; 
% 
% % Plot the basin image
% figure('Position', [100, 100, 800, 800]);
% imagesc(x_vals, v_vals, basin_image);
% colormap(combined_cmap);
% axis xy equal tight;
% xlabel('Initial Position (x)', 'FontSize', 14);
% ylabel('Initial Velocity (v)', 'FontSize', 14);
% title(sprintf('Basin of Attraction (F_0 = %.3f)', params.F0), 'FontSize', 16);
% set(gca, 'FontSize', 12);
% 
% % Add a zoomed inset to show fractal detail
% axes('Position', [0.55, 0.2, 0.3, 0.3]);
% x_zoom = [-0.4, 0.4];
% v_zoom = [0.8, 1.6];
% x_idx = find(x_vals >= x_zoom(1) & x_vals <= x_zoom(2));
% v_idx = find(v_vals >= v_zoom(1) & v_vals <= v_zoom(2));
% imagesc(x_vals(x_idx), v_vals(v_idx), basin_image(v_idx, x_idx));
% colormap(combined_cmap);
% axis xy equal tight;
% title('Zoomed Fractal Boundary', 'FontSize', 10);

% Normalize basin_image to [0, 1] for colormap mapping
basin_image_norm = (basin_image - min(basin_image(:))) / (max(basin_image(:)) - min(basin_image(:)));

% Choose a single colormap (e.g., 'turbo' or 'parula' for nice contrast)
cmap = turbo(256); % or parula(256), jet(256), etc.

% Plot the basin image with a single colormap
figure('Position', [100, 100, 800, 800]);
imagesc(x_vals, v_vals, basin_image_norm);
colormap(cmap);
colorbar;
axis xy equal tight;
xlabel('Initial Position (x)', 'FontSize', 20);
ylabel('Initial Velocity (v)', 'FontSize', 20);
title(sprintf('Basin of Attraction (F_0 = %.3f)', params.F0), 'FontSize', 22);
set(gca, 'FontSize', 18);
colormap(cmap);

% % Add a zoomed inset to show fractal detail (using same colormap)
% axes('Position', [0.55, 0.2, 0.3, 0.3]);
% x_zoom = [-0.4, 0.4];
% v_zoom = [0.8, 1.6];
% x_idx = find(x_vals >= x_zoom(1) & x_vals <= x_zoom(2));
% v_idx = find(v_vals >= v_zoom(1) & v_vals <= v_zoom(2));
% imagesc(x_vals(x_idx), v_vals(v_idx), basin_image_norm(v_idx, x_idx));
% axis xy equal tight;
% title('Zoomed Boundary', 'FontSize', 10);

