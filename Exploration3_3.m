clear; clc;

%domain
y = 0:0.01:2;

%parameters
alpha = 10;
beta = 5.621;
gamma = 8.953;
K = 2;
n = 4;

h = beta+ (alpha*(y.^n)./(K + y.^n));
l = gamma*y;

figure(1)
plot(y,h,y,l,'LineWidth',4)
h = legend('$h_+(y)$','$\gamma y$','Interpreter','Latex');
set(h,'box','off')
xlabel('$y$','Interpreter','Latex')
ylabel('rate')
title('Rate Plot at the Cusp')
set(gca,'fontsize',20)