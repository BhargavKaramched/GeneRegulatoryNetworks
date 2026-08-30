clear; clc;

%parameters
alpha = 10;
gamma = 2;
K = 1;
n = 4;

x = 3;
y = 4;

x_dot = (alpha*(y^n)/(K^n + y^n))-gamma*x;
y_dot = (alpha*(x^n)/(K^n + x^n))-gamma*y;

