clear; clc;

%parameters
global alpha K n gamma
alpha = 10;
K = 2;
n = 4;
gamma = 1;


%domain vector
y = 0:0.01:15;

%Hill function
H = alpha*(y.^n)./(K^n + y.^n);

figure(1)
subplot(2,2,1)
plot(y,H,y,gamma*y,'LineWidth',4)
h = legend('$h_+(y)$','$\gamma y$','Interpreter','Latex');
set(h,'box','off')
xlabel('$y$','Interpreter','latex')
ylabel('rate')
set(gca,'fontsize',20)

sol = ode23s(@positive_ode,[0,4],5);
sol1 = ode23s(@positive_ode,[0,4],1);

gamma = 5;

subplot(2,2,2)
plot(y,H,y,gamma*y,'LineWidth',4)
xlabel('$y$','Interpreter','latex')
ylabel('rate')
set(gca,'fontsize',20)
ylim([0 15])

sol2 = ode23s(@positive_ode,[0,2],5);


figure(1)
subplot(2,2,3)
plot(sol.x,sol.y,sol1.x,sol1.y,'LineWidth',4)
h = legend('$\gamma = 1$','$\gamma = 5$','Interpreter','Latex');
set(h,'box','off')
xlabel('time')
ylabel('$y(t)$','Interpreter','latex')
set(gca,'fontsize',20)

subplot(2,2,4)
plot(sol2.x,sol2.y,'LineWidth',4)
xlabel('time')
ylabel('$y(t)$','Interpreter','latex')
set(gca,'fontsize',20)



function out = positive_ode(t,y)

global alpha K n gamma

out = (alpha*y(1)^n/(K^n + y(1)^n)) - gamma*y(1);
end