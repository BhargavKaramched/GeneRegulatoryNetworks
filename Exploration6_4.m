% Ergodicity_PositiveFeedback.m
% Demonstrates ergodicity in a stochastic positive feedback circuit.
% Compares a population snapshot (Part A) to a single long trajectory (Part B).

clear; clc; close all;

%% 1. Define Network Parameters
beta_basal = 0.1;
alpha = 10;
K = 2;
n = 4;
gamma = 1;

%% Part A: The Multiverse (Ensemble of 5,000 Cells)
N_cells = 5000;
t_end_multiverse = 500; % Increased from 30 to allow full relaxation to steady state
final_states = zeros(N_cells, 1);

fprintf('Simulating the Multiverse (5,000 cells to t = 500)...\n');

for i = 1:N_cells
    t = 0;
    y = 0; % Initial protein count
    
    while t < t_end_multiverse
        % Propensities
        a_prod = beta_basal + (alpha * y^n) / (K^n + y^n);
        a_deg = gamma * y;
        a0 = a_prod + a_deg;
        
        % Time step and reaction
        t = t + (1/a0) * log(1/rand());
        
        if t < t_end_multiverse
            if rand() < (a_prod / a0)
                y = y + 1;
            else
                y = max(0, y - 1);
            end
        end
    end
    final_states(i) = y;
end

%% Part B: The Immortal Cell (Single Cell over 150,000 time units)
t_end_immortal = 150000;
t_transient = 1000;
sample_interval = 1; % Sample the state every 1 time unit

fprintf('Simulating the Immortal Cell (t = 150,000)...\n');

% Preallocate for speed based on estimated reaction count
est_steps = round(t_end_immortal * (alpha + 5)); 
T_hist = zeros(est_steps, 1);
Y_hist = zeros(est_steps, 1);

t = 0;
y = 0;
step = 1;
T_hist(1) = t;
Y_hist(1) = y;

while t < t_end_immortal
    % Propensities
    a_prod = beta_basal + (alpha * y^n) / (K^n + y^n);
    a_deg = gamma * y;
    a0 = a_prod + a_deg;
    
    % Time step and reaction
    t = t + (1/a0) * log(1/rand());
    
    if rand() < (a_prod / a0)
        y = y + 1;
    else
        y = max(0, y - 1);
    end
    
    step = step + 1;
    % Dynamically expand array if needed (rare if preallocated well)
    if step > length(T_hist)
        T_hist = [T_hist; zeros(100000, 1)]; 
        Y_hist = [Y_hist; zeros(100000, 1)];
    end
    
    T_hist(step) = t;
    Y_hist(step) = y;
end

% Trim arrays
T_hist = T_hist(1:step);
Y_hist = Y_hist(1:step);

% Sample the immortal cell at uniform intervals, discarding the transient
sample_times = t_transient:sample_interval:t_end_immortal;
immortal_samples = interp1(T_hist, Y_hist, sample_times, 'previous');

%% 3. Plotting the Ergodic Equivalence
figure('Position', [200, 200, 700, 500]);

% Create normalized histograms (Probability Density)
edges = -0.5:1:max([max(final_states), max(immortal_samples)]) + 0.5;

histogram(final_states, edges, 'Normalization', 'pdf', ...
    'FaceColor', '#0072BD', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
histogram(immortal_samples, edges, 'Normalization', 'pdf', ...
    'FaceColor', '#D95319', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

title('Ergodicity: Multiverse vs. Immortal Cell', 'FontSize', 14);
xlabel('Protein Count (y)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
set(gca,'fontsize',20)
legend('Multiverse', ...
       'Immortal Cell', 'Location', 'best');
grid on;