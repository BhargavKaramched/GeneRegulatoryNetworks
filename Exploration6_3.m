% Ensemble_Repressilator.m
% Simulates an ensemble of 1,000 stochastic repressilators to demonstrate 
% phase diffusion and the damping of the population mean.

clear; clc; close all;

%% 1. Define Network Parameters
alpha = 150;
K = 10;
gamma = 1;
n = 4; % Above the Hopf bifurcation threshold (deterministic oscillator)

N_cells = 1000;
t_max = 250;
t_grid = 0:0.1:t_max; % Uniform time grid for sampling

% Preallocate a matrix to hold the interpolated X counts for all cells
X_ensemble = zeros(N_cells, length(t_grid));

fprintf('Simulating %d cells... This may take a moment.\n', N_cells);

%% 2. Run the Gillespie Algorithm for the Ensemble
for cell_idx = 1:N_cells
    
    % Initialize stochastic variables for this cell
    max_steps = 100000;
    T_stoch = zeros(max_steps, 1);
    X_stoch = zeros(max_steps, 1);
    
    t_current = 0;
    x = 150; y = 0; z = 0;
    
    T_stoch(1) = t_current;
    X_stoch(1) = x;
    step = 1;
    
    while t_current < t_max && step < max_steps
        % Propensities
        a1 = (alpha * K^n) / (K^n + z^n);
        a2 = (alpha * K^n) / (K^n + x^n);
        a3 = (alpha * K^n) / (K^n + y^n);
        a4 = gamma * x;
        a5 = gamma * y;
        a6 = gamma * z;
        
        a0 = a1 + a2 + a3 + a4 + a5 + a6;
        
        % Random draws
        r1 = rand();
        r2 = rand();
        
        % Time update
        tau = (1/a0) * log(1/r1);
        t_current = t_current + tau;
        
        % Reaction selection
        rxn_thresh = r2 * a0;
        if rxn_thresh < a1
            x = x + 1;
        elseif rxn_thresh < a1 + a2
            y = y + 1;
        elseif rxn_thresh < a1 + a2 + a3
            z = z + 1;
        elseif rxn_thresh < a1 + a2 + a3 + a4
            x = max(0, x - 1);
        elseif rxn_thresh < a1 + a2 + a3 + a4 + a5
            y = max(0, y - 1);
        else
            z = max(0, z - 1);
        end
        
        step = step + 1;
        T_stoch(step) = t_current;
        X_stoch(step) = x;
    end
    
    % Trim excess
    T_stoch = T_stoch(1:step);
    X_stoch = X_stoch(1:step);
    
    % Interpolate onto uniform time grid
    % 'previous' maintains the discrete, step-like nature of the Gillespie data
    X_ensemble(cell_idx, :) = interp1(T_stoch, X_stoch, t_grid, 'previous', 'extrap');
end

%% 3. Calculate Ensemble Mean and Plot
X_mean = mean(X_ensemble, 1);

figure('Position', [150, 150, 900, 400]);

% Plot a few individual trajectories in the background to show they are oscillating
plot(t_grid, X_ensemble(1,:), 'Color', [0.7 0.7 0.7 0.5]); hold on;
plot(t_grid, X_ensemble(2,:), 'Color', [0.7 0.7 0.7 0.5]);
plot(t_grid, X_ensemble(3,:), 'Color', [0.7 0.7 0.7 0.5]);

% Plot the population mean in bold
plot(t_grid, X_mean, 'k', 'LineWidth', 3);

%title('Phase Diffusion: Ensemble Mean of 1,000 Cells', 'FontSize', 14);
xlabel('time', 'FontSize', 12);
ylabel('protein X count', 'FontSize', 12);
set(gca,'fontsize',20)
legend('single cell 1', 'single cell 2', 'single cell 3', 'mean', 'Location', 'northeast');
grid on;