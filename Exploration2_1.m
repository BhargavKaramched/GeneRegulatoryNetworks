clear; clc;

%parameters
alpha = 10;
gamma = 2;
K = 2;
n = 2;

y2 = (1/(2*gamma))*(alpha - sqrt(alpha^2-4*(K^2)*(gamma^2)));

y3 = (1/(2*gamma))*(alpha + sqrt(alpha^2-4*(K^2)*(gamma^2)));

y = 0:0.01:5;

H = alpha*(y.^n)./(K^n + y.^n);

