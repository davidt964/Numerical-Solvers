%==========================================================================
% AUTHOR: David L. Tran
%
% Newton's method for 1 or 2 nonlinear equations.
%
% DESCRIPTION: Performs Newton's method iteratively on either 1 or 2
% nonlinear equations that are to be specified by the user. The program
% continues to iterate until n reaches IT_MAX or the 2-norm of the residual
% is less than a thresshold epsilon value. The residual and Jacobian
% matrices are obtained and then used to update the next guess.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables

% Stopping Criteria
epsilon = 10^(-12);       %tolerance/residual stopping criterion
IT_MAX = 100;             %maximum number of iterations before program termination

% USER-SPECIFIED PROBLEM DIMENSIONS
DIM = 2;                  %number of dimensions (1 for one eqn and 2 for two eqns)

% Initial Guesses
if DIM == 1
    % one equation
    x_0 = 0;                     %initial guess (1 nonlinear eqn)
    x_n = [x_0];                 %array to store all iterative solutions
    r_n = [];                    %residual array
elseif DIM == 2
    % two equations
    x_0 = [4; -4];               %initial guess (2 nonlinear eqns)
    x_n = [x_0(1)];              %array to store all x iterative solutions
    y_n = [x_0(2)];              %array to store all y iterative solutions
end

n = 1;                    %counter
J_n = [];                 %Jacobian array

%% Iterative Procedure
%enter the loop; norm_2 is overwritten anyway in the loop.

if DIM == 1
    %Calculate residual
    r_n = calcResidual(x_n(n), DIM);
    
    %Calculate Jacobian
    J_n = calcJacobian(x_n(n), DIM);

    %Calculate the 2-norm of the residual
    norm_2 = sqrt(r_n^2);
elseif DIM == 2
    %Calculate residual
    r_n = calcResidual([x_n(n), y_n(n)], DIM);
    
    %Calculate Jacobian
    J_n = calcJacobian([x_n(n), y_n(n)], DIM);
    
    %Calculate the 2-norm of the residual
    norm_2 = sqrt(r_n(1)^2 + r_n(2)^2);
end

while norm_2 >= epsilon && n < IT_MAX
    if DIM == 1
        %Obtain Delta x_n
        Delta_xn(n) = J_n(n)\r_n(n);
    
        %Update x_n
        x_n(n+1) = x_n(n) + Delta_xn(n);

        %Increment counter
        n = n + 1;
        
        %Calculate residual
        r_n = calcResidual(x_n(n), DIM);
        
        %Calculate Jacobian
        J_n(n) = calcJacobian(x_n(n), DIM);

        %Calculate the 2-norm of the residual
        norm_2 = sqrt(r_n^2);
    
    elseif DIM == 2
        %Obtain Delta x_n
        Delta_xn = J_n\r_n;

        %Update x_n
        x_n(n+1) = x_n(n) + Delta_xn(1);
        y_n(n+1) = y_n(n) + Delta_xn(2);

        %Increment counter
        n = n + 1;
        
        %Calculate residual
        r_n = calcResidual([x_n(n), y_n(n)], DIM);
    
        %Calculate Jacobian
        J_n = calcJacobian([x_n(n), y_n(n)], DIM);

        %Calculate the 2-norm of the residual
        norm_2 = sqrt(r_n(1)^2 + r_n(2)^2);
    
    end

end


%% Command Window Display
disp('%======================Iterative Solutions (All n)======================');
if DIM == 1
    for i = 1:n
        fprintf('n = %d  ||  x_n = %.9f \n', i-1, x_n(i));
    end
elseif DIM == 2
    for i = 1:n
        fprintf('n = %d  ||  x_n = %.9f  ||  y_n = %.9f\n', i-1, x_n(i), y_n(i));
    end
end

%% Functions
function [r_xn] = calcResidual(x_n, DIM)
% Calculates the residual r(x_n) based on current iterative solution x_n.


% Anonymous Functions
if DIM == 1
    %add negative in front of function for residual
    resid = @(x) -(-x^5 + sin(5*x) - cos(8*x) + 3*exp(-x));

    r_xn = resid(x_n);
elseif DIM == 2
    resid_1 = @(x,y) -(4*x^2 + x^3*y-6);
    resid_2 = @(x,y) -(x^3*y^4 - y^2*cos(5*y) + 1);

    r_xn = [resid_1(x_n(1), x_n(2)); resid_2(x_n(1), x_n(2)) ];
end


end

function [J_xn] = calcJacobian(x_n, DIM)
% Calculates the Jacobian matrix J(x_n) based on current iterative solution
% x_n.


if DIM == 1
    %derivative of f(x)
    f_prime = @(x) -5*x^4 + 5*cos(5*x) + 8*sin(8*x) - 3*exp(-x);

    %the jacobian is simply the derivative of f(x) wrt x evaluated at x_n.
    J_xn = f_prime(x_n);
elseif DIM == 2
    %partial derivatives of f_1(x,y) and f_2(x,y) wrt x and y
    f_1_x = @(x,y) 8*x + 3*x^2*y;
    f_1_y = @(x,y) x^3;
    f_2_x = @(x,y) 3*x^2*y^4;
    f_2_y = @(x,y) 4*x^3*y^3-2*y*cos(5*y)+5*y^2*sin(5*y);
    
    %fill the Jacobian matrix
    J_xn = [f_1_x(x_n(1), x_n(2)), f_1_y(x_n(1), x_n(2)); f_2_x(x_n(1), x_n(2)), f_2_y(x_n(1), x_n(2))];
end


end
