clear; clc;

%parameters
K = 2;
gamma =4;
n = 10;
alpha = 30;

syms y
h = (alpha*K^n/(K^n + y^n)) - gamma*y==0;

vpasolve(h,y)