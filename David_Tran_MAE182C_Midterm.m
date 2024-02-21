%==========================================================================
% AUTHOR: David L. Tran
%
% MAE 182C MIDTERM TAKE HOME PORTION
% DATE: NOVEMBER 06, 2023
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variable Initialization

x_i = [2, 5, 8, 11, 14, 17, 20]; %x data points
y_i = [-3, 4, 0, 5, 3, 3, -1]; %y data points
x_n = []; %x nodes

h = x_i(2) - x_i(1);  %constant


%% Anonymous Function
% First derivatives of the 3 Lagrangean basis functions
% Evaluate each derivative with inputs xEval and index i
dL1 = @(xEval, x, i)  (2*xEval - x(i+2) - x(i+1)) / ((x(i) - x(i+1) )* (x(i) - x(i+2)) );
dL2 = @(xEval, x, i)  (2*xEval - x(i) - x(i+2)) / ((x(i+1) - x(i)) * (x(i+1) - x(i+2)) ); 
dL3 = @(xEval, x, i)  (2*xEval - x(i+1) - x(i)) / ((x(i+2) - x(i)) * (x(i+2) - x(i+1)) );   

%% Loops/Iterative Process

% [x0, x2]
[xPlot, yPlot] = calcLagrNodes(x_i(1:3), y_i(1:3), 1);

% [x2, x3]
q0_prime = y_i(1) * dL1(x_i(3), x_i, 1) + y_i(2) * dL2(x_i(3), x_i, 1) + y_i(3) * dL3(x_i(3), x_i, 1);
a = [0.5; 0.5; 0.5];
b = [1; 2; 1];
f = [29/6; -7/3; 4/3];
[M_i] = calcThomas(a, b, f);
[xCubic, yCubic] = calcSpline(x_i, y_i, M_i, 3, h);

% [x3, x4]
q4_prime = y_i(5) * dL1(x_i(5), x_i, 5) + y_i(6) * dL2(x_i(5), x_i, 5) + y_i(7) * dL3(x_i(5), x_i, 5);
[xCubic2, yCubic2] = calcSpline(x_i, y_i, M_i, 4, h);

% [x4, x6]
[xPlot2, yPlot2] = calcLagrNodes(x_i, y_i, 5);


%% Integral Approximation
% [x0, x2]
S_1 = h/3 * (y_i(1) + 4*y_i(2) + y_i(3));

% [x2, x3]
s2_mid = ((h/2)^3*M_i(3) + (h/2)^3*M_i(4)) / 18 + (1.5*y_i(3) + 1.5*y_i(4)) / 3 - (1/2) * (1.5 * M_i(3) + 1.5*M_i(4));
S_2 = h/3 * (y_i(3) + 4 * s2_mid + y_i(4) );

% [x3, x4]
s3_mid = ((h/2)^3*M_i(4) + (h/2)^3*M_i(5)) / 18 + (1.5*y_i(4) + 1.5*y_i(5)) / 3 - (1/2) * (1.5 * M_i(4) + 1.5*M_i(5));
S_3 = h/3 * (y_i(4) + 4 * s3_mid + y_i(5) );

% [x4, x6]
S_4 = h/3 * (y_i(5) + 4*y_i(6) + y_i(7));

simp_approx = S_1 + S_2 + S_3 + S_4;

fprintf('Simpson''s Rule Result: %f\n', simp_approx);


%% Plots
figure;
scatter(x_i, y_i, 30, 'b', 'filled'); hold on;
plot(xPlot, yPlot, 'r', 'LineWidth',2);
plot(xCubic, yCubic, 'r', 'LineWidth',2);
plot(xCubic2, yCubic2, 'r', 'LineWidth',2);
plot(xPlot2, yPlot2, 'r', 'LineWidth',2);
xlabel('$x$','Interpreter','LaTeX');    
ylabel('$y$','Interpreter','LaTeX');    
set(gca,'LineWidth',2,'FontSize',18);   
xlim([x_i(1) x_i(end)]);


%% Functions
function [xPlot, yPlot] = calcLagrNodes(x, y, i)
% Calculates the plot nodes for the Lagrangean quadratic interpolants in
% subinterval [x_i, x_i+2].
xPlot = linspace(x(i), x(end), (x(end) - x(i)) * 100 );

yPlot = y(i) .* ((xPlot - x(i+1)) .* (xPlot - x(i+2))) ./ ((x(i) - x(i+1)) .* (x(i) - x(i+2))) ...
    + y(i+1) .* ((xPlot - x(i)) .* (xPlot - x(i+2))) ./ ((x(i+1) - x(i)) .* (x(i+1) - x(i+2))) ...
    + y(i+2) .* ((xPlot - x(i)) .* (xPlot - x(i+1))) ./ ((x(i+2) - x(i)) .* (x(i+2) - x(i+1)));

end

function [M_i] = calcThomas(a, b, f)
% Solves the 3 x 3 tridiagonal linear system AM = f to calculate and return
% the M_i weights of the cubic spline in subinterval [x_2, x_4] with
% clamped BCs.

m_i = [];
beta_i = [b(1)];
g_i = [f(1)];
M_i = [];

for i = 2:3
    if i == 2
        m_i(i) = a(i) / b(i-1);
        beta_i(i) = b(i) - m_i(i) * a(i);
        g_i(i) = f(i) - m_i(i) * f(i-1);
    else
        m_i(i) = a(i) / beta_i(i-1);
        beta_i(i) = b(i) - m_i(i) * a(i);
        g_i(i) = f(i) - m_i(i) * g_i(i-1);
    end
end


for i = 3:-1:1
    if i == 3
        M_i(i) = g_i(i) / beta_i(i);
    else
        M_i(i) = (g_i(i) - a(i) * M_i(i+1)) / beta_i(i);
    end

end

M_i = [0, 0, M_i, 0, 0];

end

function [xPlot, yPlot] = calcSpline(x, y, M, i, h)
xPlot = linspace(x(i), x(i+1), (x(i+1) - x(i)) * 100 );

yPlot = ( (x(i+1) - xPlot).^3 .* M(i) + (xPlot - x(i)).^3.*M(i+1))./ (6 .* h)  ...
    + ( y(i).*(x(i+1) - xPlot) + y(i+1).*(xPlot - x(i)) ) ./ (h)    ...
    - (1./2) .* ( (x(i+1) - xPlot).*M(i) + (xPlot - x(i)).*M(i+1) );

end