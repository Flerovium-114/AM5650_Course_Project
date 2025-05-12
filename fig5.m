% Figure 5: Critical Force F_CR as a function of various parameters
clear; clc; close all

% Fixed parameters
zeta = 0.02;
omega = 1.0;
Omega_range = linspace(0.1, 2, 200);
alpha_list = [-0.9, -0.5, -0.1];
A_list = [1.0, 1.5, 2.0];
gamma_beta_list = [0.5, 1.0, 1.5];

% Melnikov threshold formula (from Eq. 20 in the paper)
F_CR = @(Omega, alpha, A, gamma_beta) ...
    abs(4 * zeta * omega^3 * (alpha + (1 - alpha) * A).^2 ./ ...
    (gamma_beta .* (1 - alpha) .* A.^2 .* Omega * pi) .* ...
    sinh(Omega * pi ./ (2 * omega * sqrt((alpha + (1 - alpha) * A)/2))));

%% (a) F_CR vs Omega for different alpha
figure;
hold on;
for i = 1:length(alpha_list)
    alpha = alpha_list(i);
    plot(Omega_range, F_CR(Omega_range, alpha, 1.0, 1.0), 'LineWidth', 2);
end
xlabel('\Omega');
ylabel('Critical Force F_{CR}');
legend(arrayfun(@(a) sprintf('\\alpha=%.1f', a), alpha_list, 'UniformOutput', false));
title('F_{CR} vs \Omega for different \alpha');
grid on;

%% (b) F_CR vs alpha for fixed Omega
figure;
alpha_range = linspace(-1.0, 0.0, 200);
FCR_alpha = F_CR(1.0, alpha_range, 1.0, 1.0);
plot(alpha_range, FCR_alpha, 'LineWidth', 2);
xlabel('\alpha');
ylabel('Critical Force F_{CR}');
title('F_{CR} vs \alpha for \Omega=1.0');
grid on;

%% (c) F_CR vs Omega for different A
figure;
hold on;
for i = 1:length(A_list)
    A = A_list(i);
    plot(Omega_range, F_CR(Omega_range, -0.5, A, 1.0), 'LineWidth', 2);
end
xlabel('\Omega');
ylabel('Critical Force F_{CR}');
legend(arrayfun(@(A) sprintf('A=%.1f', A), A_list, 'UniformOutput', false));
title('F_{CR} vs \Omega for different A');
grid on;

%% (d) F_CR vs A for fixed Omega
figure;
A_range = linspace(0.5, 2.0, 200);
FCR_A = F_CR(1.0, -0.5, A_range, 1.0);
plot(A_range, FCR_A, 'LineWidth', 2);
xlabel('A');
ylabel('Critical Force F_{CR}');
title('F_{CR} vs A for \alpha=-0.5, \Omega=1.0');
grid on;

%% (e) F_CR vs Omega for different (gamma+beta)
figure;
hold on;
for i = 1:length(gamma_beta_list)
    gb = gamma_beta_list(i);
    plot(Omega_range, F_CR(Omega_range, -0.5, 1.0, gb), 'LineWidth', 2);
end
xlabel('\Omega');
ylabel('Critical Force F_{CR}');
legend(arrayfun(@(gb) sprintf('\\gamma+\\beta=%.1f', gb), gamma_beta_list, 'UniformOutput', false));
title('F_{CR} vs \Omega for different (\gamma+\beta)');
grid on;

%% (f) F_CR vs (gamma+beta) for fixed Omega
figure;
gb_range = linspace(0.5, 2.0, 200);
FCR_gb = F_CR(1.0, -0.5, 1.0, gb_range);
plot(gb_range, FCR_gb, 'LineWidth', 2);
xlabel('\gamma+\beta');
ylabel('Critical Force F_{CR}');
title('F_{CR} vs (\gamma+\beta) for \alpha=-0.5, \Omega=1.0');
grid on;
