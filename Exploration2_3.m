clear; clc;

alpha1 = 10;
alpha2 = 10;
K1 = 1;
K2 = 1;
gamma1 = 2;
gamma2 = 2;
n = 4;

x = 0:0.01:6;
y = (alpha2/gamma2)*(x.^n)./(K2^n + x.^n);

y1 = 0:0.01:6;
x1 = (alpha1/gamma1)*(y1.^n)./(K1^n + y1.^n);

figure(1)
subplot(1,3,1)
plot(x,y,x1,y1,'LineWidth',4)
title('$\gamma_y = 2$','Interpreter','Latex')
h = legend('$y$ null cline','$x$ null cline','Interpreter','Latex');
set(h,'box','off')
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
set(gca,'fontsize',20)

gamma2 = 5;
y = (alpha2/gamma2)*(x.^n)./(K2^n + x.^n);

subplot(1,3,2)
plot(x,y,x1,y1,'LineWidth',4)
title('$\gamma_y = 5$','Interpreter','Latex')
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
set(gca,'fontsize',20)

gamma2 = 15;
y = (alpha2/gamma2)*(x.^n)./(K2^n + x.^n);

subplot(1,3,3)
plot(x,y,x1,y1,'LineWidth',4)
title('$\gamma_y = 15$','Interpreter','Latex')
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
set(gca,'fontsize',20)