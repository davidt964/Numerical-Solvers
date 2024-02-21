%==========================================================================
% AUTHOR: David L. Tran
%
% 1st order, nonlinear IVP with numerical integration with trapezoidal
% method and integration error estimates that enable higher-order numerical
% integration schemes
%
% DESCRIPTION: The solution to the abstract 1st order IVP 
% y'(x) = (1 + x) / (1 + y(x)) in a given domain and IVP is sought with a
% given true solution. The trapezoidal rule is applied and a corrective 
% method is employed by adding the asymptotic error estimate to the
% integration estimate I_h(f). The first derivative of the integrand is 
% obtained. Due to the implicit nature of this new method, a rootfinding
% problem is written and solved using Newton's method. Heun's method is
% used to obtain an accurate initial guess.
%
%==========================================================================


%% Clear Cache
clc; close all; clearvars;

%% Variables

h = 1/2;    %step size

p = 2; %constant for power

EPS = 10^(-12); %newtons method convergence
IT_MAX = 150; %maximum number of iterations

x_low = -3; %x lower limit
x_high = 4; %x upper limit

n_n = (x_high - x_low) / h + 1;  %number of nodes
n_SI = n_n - 1; %number of subintervals

x_nodes = linspace(x_low,x_high,n_n); %x arr of nodes
y_nodes = zeros(n_n,1); %x arr of nodes
y_nodes(1) = 2; %initial condition


%% Anonymous Functions

y_true = @(x) sqrt(x.^2+2.*x+6)-1; %analytical soln
y_prime_true = @(x) (x + 1) ./ (sqrt(x.^2 + 2.*x + 6)); %true derivative

y_prime_ode = @(x,y) (1 + x) / (1 + y); %ode 
df_dy = @(x,y) -(x + 1) / (y + 1)^2;
d2f_dy2 = @(x,y) (2 * (x + 1)) / (y + 1)^3;
df_dx = @(x,y) ( (1 + y) - (1 + x)*y_prime_ode(x,y) ) / ( (1 + y)^2 ) + 1 / (y+1) * y_prime_ode(x,y);
dF_dy = @(h, x, y, n_SI) 1 - h/2 * (df_dy(x,y)) - h^2/(12 * n_SI^2) * d2f_dy2(x,y);
c = @(h, x, y, i) -h^2 / 12 * (df_dx(x(i+1), y(i+1)) - df_dx(x(i),y(i)));


%% Loop

for i = 1:n_SI
    
    %reset j
    j = 1;

    %forward euler for part of heun's method
    f_FE = y_nodes(i) + h * y_prime_ode(x_nodes(i), y_nodes(i));

    %initial guess using heun's method
    y_nodes(i+1) = y_nodes(i) + h/2 * (y_prime_ode(x_nodes(i), y_nodes(i)) + f_FE);

    %store initial integral approximation from heun's method
    I_h = h/2 * (y_prime_ode(x_nodes(i), y_nodes(i)) + f_FE);

    %calculate asymptotic error estimate
    E_T = c(h, x_nodes, y_nodes, i) / n_SI^p;

    %calculate residual
    r = -(y_nodes(i+1) - y_nodes(i) - I_h - E_T);

    %calculate F'(y_n+1)
    F_prime = dF_dy(h, x_nodes(i+1), y_nodes(i+1), n_SI);

    %calculate norm
    norm_2 = sqrt(r^2);

    while norm_2 > EPS && j < IT_MAX
        %calculate change (scalars)
        Delta_y = r/F_prime;
        
        %update counter
        j = j + 1;

        %update y_n+1
        y_nodes(i+1) = y_nodes(i+1) + Delta_y;

        %calculate integrand approximation
        I_h = h / 2 * ( y_prime_ode(x_nodes(i+1),y_nodes(i+1)) + y_prime_ode(x_nodes(i),y_nodes(i)));

        %calculate asymptotic error estimate
        E_T = c(h, x_nodes, y_nodes, i) / n_SI^p;

        %calculate residual
        r = -(y_nodes(i+1) - y_nodes(i) - I_h - E_T);
        
        %calculate norm
        norm_2 = sqrt(r^2);

    end
end

for i = 1:n_n
    if i == 1
        %use forward difference for first derivative
        y_prime_h(i) = (y_nodes(i+1) - y_nodes(i)) / h;
    elseif i == n_n
        %use backward difference for last derivative
        y_prime_h(i) = (y_nodes(i) - y_nodes(i-1)) / h;
    else
        %central difference
        y_prime_h(i) = (y_nodes(i+1) - y_nodes(i-1)) / (2 * h);
    end
end

%% Print
error = abs(y_true(x_high) - y_nodes(end));
fprintf('Error at x = %f    for    h = %f:\n', x_high, h);
fprintf('|y(x) - y_h(x)|  =  %e     \n',error);


%% Plots
figure(1);
set(gcf,'Position',[0 25 750 800]);    
plot(x_nodes,y_true(x_nodes), 'b', 'LineWidth',2);
hold on; 
plot(x_nodes,y_nodes, 'r', 'LineWidth',2);
title("Improved Numerical Method y(x)");                 
xlabel('$x$','Interpreter','LaTeX');
xlim([x_low x_high]);
ylabel('$y$','Interpreter','LaTeX');    
set(gca,'LineWidth',2,'FontSize',18);
legend({'$y(x)$', strcat('$y_h(x) \;\mathrm{ with }\; h = $ ',num2str(h))},'Interpreter','latex','location','best');

figure(2);
set(gcf,'Position',[0 25 750 800]);   
plot(x_nodes,y_prime_true(x_nodes), 'b', 'LineWidth',2);
hold on;  
plot(x_nodes,y_prime_h, 'r', 'LineWidth',2);
title("Improved Numerical Method y'(x)");                 
xlabel('$x$','Interpreter','LaTeX');
xlim([x_low x_high]);
ylabel('$y^\prime$','Interpreter','LaTeX');    
set(gca,'LineWidth',2,'FontSize',18);
legend({'$y''(x)$', strcat('$y_h''(x) \;\mathrm{ with }\; h = $ ',num2str(h))},'Interpreter','latex','location','best');
