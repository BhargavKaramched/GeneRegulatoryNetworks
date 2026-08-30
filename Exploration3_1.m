clear; clc;

global beta alpha K n gamma

beta = 5;
alpha = 10;
K = 1.189;
n = 4;
gamma = 8.2;

T = 10;

y0 = 0.7463;

sol = ode23s(@basal_positive_feedback,[0,T],y0);

sol.y(end)


function out = basal_positive_feedback(t,y)
global beta alpha K n gamma

out = beta + (alpha*(y(1)^n)/(K^n + y(1)^n)) - gamma*y(1);
end