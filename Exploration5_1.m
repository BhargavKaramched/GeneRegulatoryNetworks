% DDE_Oscillator_Analysis.m
% Simulates a gene regulatory network with delayed negative feedback.
% Analyzes the Hopf bifurcation, amplitude saturation, and frequency tuning.

clear; clc; close all;

%% 1. Define Network Parameters
alpha = 10;
gamma = 1;
K = 2;
n = 4;

% Define the DDE: dy/dt = (alpha * K^n) / (K^n + y(t-tau)^n) - gamma * y(t)
% In MATLAB dde23, Z(:,1) represents the delayed variable y(t-tau)
dde_fun = @(t, y, Z) (alpha * K^n) / (K^n + Z(1,1)^n) - gamma * y(1);

%% Part A: Time Series & Identifying the Bifurcation
tau_examples = [0.6, 1.0, 1.5]; % Sweeping past the critical delay
tspan = [0 50];
hist = 0; % Initial history 

figure('Position', [800, 800, 800, 800]);
subplot(2,2,[1 2])
hold on;
colors = {'#0072BD', '#D95319', '#EDB120'};

for i = 1:length(tau_examples)
    tau = tau_examples(i);
    sol = dde23(dde_fun, tau, hist, tspan);
    plot(sol.x, sol.y, 'LineWidth', 4, 'Color', colors{i}, ...
        'DisplayName', sprintf('\\tau = %.1f', tau));
    set(gca,'fontsize',20)
end
title('Emergence of the Biological Clock');
xlabel('time'); ylabel('$y(t)$','Interpreter','Latex');
legend('Location', 'best'); grid on;

%% Part B & C: Sweeping Delays for Amplitude and Frequency
tau_sweep = [1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 8.0];
amplitudes = zeros(size(tau_sweep));
periods = zeros(size(tau_sweep));

% Run a longer time span to ensure we only measure stable, sustained oscillations
tspan_long = [0 200]; 
eval_time = linspace(100, 200, 5000); % Only evaluate the last 100 time units

for i = 1:length(tau_sweep)
    tau = tau_sweep(i);
    sol = dde23(dde_fun, tau, hist, tspan_long);
    y_eval = deval(sol, eval_time);
    
    % Find peaks (maximums)
    [peaks, locs] = findpeaks(y_eval, eval_time);
    
    if length(peaks) >= 2
        % Find troughs (minimums) by inverting the signal
        [troughs, ~] = findpeaks(-y_eval, eval_time);
        troughs = -troughs;
        
        % Calculate Amplitude and Period
        amplitudes(i) = mean(peaks) - mean(troughs);
        periods(i) = mean(diff(locs));
    else
        % If no peaks are found, the system is damped (pre-bifurcation)
        amplitudes(i) = 0;
        periods(i) = NaN; 
    end
end

%% Plotting Part B & C Results
%figure('Name', 'Part B & C: Oscillator Analytics', 'Position', [150, 150, 800, 400]);

% Plot Amplitude vs Tau
subplot(2,2,3);
plot(tau_sweep, amplitudes, '-o', 'LineWidth', 4, 'MarkerFaceColor', 'b');
%title('Part B: Amplitude Saturation');
xlabel('\tau');
ylabel('amplitude');
set(gca,'fontsize',20)
grid on;

% Plot Period vs Tau
subplot(2,2,4);
plot(tau_sweep, 1./periods, '-o', 'LineWidth', 4, 'MarkerFaceColor', 'r','MarkerSize',5);
%title('Part C: Period Tuning');
xlabel('\tau');
ylabel('frequency');
set(gca,'fontsize',20)
grid on;