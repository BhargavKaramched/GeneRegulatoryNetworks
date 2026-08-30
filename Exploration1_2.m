clear; clc;

%degradation parameter
G = [0.1, 0.5, 1.0, 2.0];
LG = length(G);

%domain
tt = 0:0.01:10;

for j = 1:LG
    gamma = G(j);

    figure(1)
    subplot(1,2,1)
    plot(tt,10*exp(-gamma.*tt),'LineWidth',4)
    h = legend('$\gamma = 0.1$','$\gamma = 0.5$','$\gamma = 1.0$','$\gamma = 2.0$','Interpreter','Latex');
    set(h,'box','off')
    xlabel('time')
    ylabel('$y(t)$','Interpreter','Latex')
    set(gca,'fontsize',20)
    hold on
end

T2 = log(2)./G;

gamma = 0.1:0.1:2;
t2 = log(2)./gamma;

figure(1)
subplot(1,2,2)
plot(G,T2,'bx',gamma,t2,'k--','MarkerSize',20,'LineWidth',4)
h = legend('data','$\frac{\log{2}}{\gamma}$','Interpreter','latex');
set(h,'box','off')
xlabel('$\gamma$','Interpreter','Latex')
ylabel('half life')
set(gca,'fontsize',20)

