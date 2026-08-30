% Activation_Inhibition_Delay.m
% Generates a 4-panel figure showing the destabilization of a 
% 2-node activation-inhibition network as transcriptional delay increases.

clear; clc; close all;

%% 1. Define Network Parameters
alpha_x = 10; 
alpha_y = 10;
gamma_x = 2; 
gamma_y = 2;
K_x = 2; 
K_y = 2;
n = 4;

% Setup the delay sweep and simulation timeframe
tau_values = [0.095, 0.4, 0.7, 1.2];
tspan = [0 30];
hist = [0; 0]; % Initial history for X and Y

% Initialize the figure
figure('Position', [100, 100, 900, 700]);

%% 2. Run the DDE Sweep
for i = 1:length(tau_values)
    tau = tau_values(i);
    
    % Define the DDE system
    % Y(1) = x(t), Y(2) = y(t)
    % Z(1,1) = x(t-tau), Z(2,1) = y(t-tau)
    % Repression of X has delayed Y. Activation of Y has instantaneous X.
    dde_fun = @(t, Y, Z) [ ...
        (alpha_x * K_y^n) / (K_y^n + Z(2,1)^n) - gamma_x * Y(1); ...
        (alpha_y * Z(1,1)^n) / (K_x^n + Z(1,1)^n) - gamma_y * Y(2) ...
    ];

    % Solve the DDE
    sol = dde23(dde_fun, tau, hist, tspan);
    
    % Plotting
    subplot(2, 2, i);
    plot(sol.x, sol.y(1,:), 'LineWidth', 4, 'Color', '#0072BD'); hold on;
    plot(sol.x, sol.y(2,:), 'LineWidth', 4, 'Color', '#D95319');
    
    title(sprintf('\\tau = %.1f', tau), 'FontSize', 14);
    xlabel('time', 'FontSize', 12); 
    ylabel('concentration', 'FontSize', 12);
    set(gca,'fontsize',20)
    
    % Add limits to keep the visual scale identical across subplots
    ylim([0 6]);
    grid on;
    
    if i == 1
        legend('X', 'Y', 'Location', 'best');
    end
end