clear; clc;

%parameters
alpha = 10;

%setting up the spatial vector
y = 0:0.01:4;

%loop 1 - varying k
b = 2;

for k = 1:4
    %Boltzmann function
    B = alpha./(1 + exp(-k.*(y-b)));

    %hyperbolic tangent
    T = (alpha/2)*(1+tanh(k.*(y-b)));

    figure(1)
    subplot(2,2,1)
    plot(y,B,'LineWidth',4)
    xlabel('$y$','Interpreter','latex')
    ylabel('$B(y)$','Interpreter','latex')
    h = legend('k = 1', 'k = 2', 'k =3', 'k = 4');
    set(h,'box','off')
    set(gca,'fontsize',20)
    hold on

    subplot(2,2,2)
    plot(y,T,'LineWidth',4)
    xlabel('$y$','Interpreter','latex')
    ylabel('$T(y)$','Interpreter','latex')
    set(gca,'fontsize',20)
    hold on
end

%loop 2 - varying b
k = 4;

for b = 1:0.5:2.5
     %Boltzmann function
    B = alpha./(1 + exp(-k.*(y-b)));

    %hyperbolic tangent
    T = (alpha/2)*(1+tanh(k.*(y-b)));

    figure(1)
    subplot(2,2,3)
    plot(y,B,'LineWidth',4)
    xlabel('$y$','Interpreter','latex')
    ylabel('$B(y)$','Interpreter','latex')
    h = legend('b = 1', 'b = 1.5', 'b =2', 'b = 2.5');
    set(h,'box','off')
    set(gca,'fontsize',20)
    hold on

    subplot(2,2,4)
    plot(y,T,'LineWidth',4)
    xlabel('$y$','Interpreter','latex')
    ylabel('$T(y)$','Interpreter','latex')
    set(gca,'fontsize',20)
    hold on
end