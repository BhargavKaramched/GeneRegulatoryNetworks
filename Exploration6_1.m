% Stochastic_Repressilator.m
% Demonstrates the emergence of noise-induced quasi-cycles in a 
% genetic repressilator that is deterministically stable.

clear; clc; close all;

%% 1. Define Network Parameters
alpha = 150;
K = 10;
gamma = 1;
n = 1.8; % Below the critical Hopf bifurcation threshold

% Initial conditions (X, Y, Z) and simulation time
X0 = [10; 80; 120]; 
tspan = [0 50];

%% 2. Deterministic ODE Solution
% X is repressed by Z, Y is repressed by X, Z is repressed by Y
ode_func = @(t, Y) [
    (alpha * K^n) / (K^n + Y(3)^n) - gamma * Y(1);
    (alpha * K^n) / (K^n + Y(1)^n) - gamma * Y(2);
    (alpha * K^n) / (K^n + Y(2)^n) - gamma * Y(3)
];

[t_det, Y_det] = ode45(ode_func, tspan, X0);

%% 3. Stochastic Gillespie Algorithm (SSA)
% Initialize stochastic variables
max_steps = 150000; % Preallocate memory for speed
T_stoch = zeros(max_steps, 1);
Y_stoch = zeros(max_steps, 3);

t_current = 0;
Y_current = X0;
T_stoch(1) = t_current;
Y_stoch(1, :) = Y_current';
step = 1;

while t_current < tspan(2) && step < max_steps
    x = Y_current(1);
    y = Y_current(2);
    z = Y_current(3);
    
    % Calculate Propensities (Reaction rates)
    % a(1)-a(3): Transcription/Translation
    a = zeros(6,1);
    a(1) = (alpha * K^n) / (K^n + z^n);
    a(2) = (alpha * K^n) / (K^n + x^n);
    a(3) = (alpha * K^n) / (K^n + y^n);
    
    % a(4)-a(6): Degradation
    a(4) = gamma * x;
    a(5) = gamma * y;
    a(6) = gamma * z;
    
    a0 = sum(a);
    
    % Draw random numbers to determine time step and next reaction
    r1 = rand();
    r2 = rand();
    
    % Time until next reaction
    tau = (1/a0) * log(1/r1);
    t_current = t_current + tau;
    
    % Determine which specific reaction occurred
    cumulative_a = cumsum(a);
    rxn = find(cumulative_a >= r2 * a0, 1);
    
    % Update the protein counts based on the reaction
    if rxn == 1;     Y_current(1) = Y_current(1) + 1;
    elseif rxn == 2; Y_current(2) = Y_current(2) + 1;
    elseif rxn == 3; Y_current(3) = Y_current(3) + 1;
    elseif rxn == 4; Y_current(1) = max(0, Y_current(1) - 1);
    elseif rxn == 5; Y_current(2) = max(0, Y_current(2) - 1);
    elseif rxn == 6; Y_current(3) = max(0, Y_current(3) - 1);
    end
    
    step = step + 1;
    T_stoch(step) = t_current;
    Y_stoch(step, :) = Y_current';
end

% Trim excess preallocated array zeros
T_stoch = T_stoch(1:step);
Y_stoch = Y_stoch(1:step, :);

%% 4. Plotting the Results
figure('Position', [100, 100, 1000, 450]);


% Subplot 1: The xy Phase Plane
subplot(1, 2, 1);
plot(Y_stoch(:,1), Y_stoch(:,2), 'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',3); hold on;
plot(Y_det(:,1), Y_det(:,2), 'k', 'LineWidth', 4);
% Mark the deterministic equilibrium
plot(Y_det(end,1), Y_det(end,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

title('x-y phase plane', 'FontSize', 14);
xlabel('protein X', 'FontSize', 12);
ylabel('protein Y', 'FontSize', 12);
set(gca,'fontsize',20)
legend('stochastic (Gillespie)', 'deterministic (ODE)', 'Equilibrium', 'Location', 'best');
grid on;

% Subplot 2: Time Series overlay
subplot(1, 2, 2);
stairs(T_stoch, Y_stoch(:,1), 'Color', [0.8500 0.3250 0.0980 0.7],'LineWidth',3); hold on;
plot(t_det, Y_det(:,1), 'k', 'LineWidth', 2.5);

title('time series for X', 'FontSize', 14);
xlabel('time', 'FontSize', 12);
ylabel('concentration / count', 'FontSize', 12);
set(gca,'fontsize',20)
legend('stochastic (Gillespie)', 'deterministic (ODE)', 'Location', 'best');
grid on;