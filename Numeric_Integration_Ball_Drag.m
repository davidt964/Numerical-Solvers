%==========================================================================
% AUTHOR: David L. Tran
%
% Performs numeric integration using 1 of 5 implemented methods (Forward
% Euler, Backward Euler, Trapezoidal Method, Heun's Method, and 4th-order
% Runge-Kutta method) on a ball descending in a fluid with drag and gravity.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables/Constants

h = 0.2;            %time step size

d = 0.2;            %ball diameter [m]
C_D = 0.47;         %drag coefficient
g = 9.81;           %gravitational acceleration [m/s^2]

% Stopping Criterion
EPS = 10^(-12);     %tolerance
IT_MAX = 1000;      %max iterations

rho_rel = 0.9;      %relative density (rho_ball/rho_water < 1)

n_SI = [10, 20, 40, 80];    %number of subintervals for several iterations


% User-Selected Method
whichMethod = 1;                            % 1 - Forward Euler method
                                            % 2 - Backward Euler method
                                            % 3 - Trapezoidal method
                                            % 4 - Heun's method
                                            % 5 - Runge-Kutta method

%Recast parameters
alpha = -(1 - 1/rho_rel)*g; 
beta = 3*C_D / (4*rho_rel*d);

%initial condition/guesses
v_0 = 5;            %initial speed [m/s]
t_0 = 1.8;          %initial time guess [s]

j_vals = [1, 2, 3, 4, 5];

%% Loops
%calculate b
n = 1;
r_n(n) = calcResidual(alpha, beta, v_0, t_0);
J_n(n) = calcJacobian(alpha, beta, v_0, t_0);
norm_2 = sqrt(r_n(n)^2);
t_n = [t_0];

while norm_2 > EPS && n < IT_MAX
        
        %Obtain Delta x_n
        Delta_tn(n) = J_n(n)\r_n(n);
    
        %Update x_n
        t_n(n+1) = t_n(n) + Delta_tn(n);

        %Increment counter
        n = n + 1;
        
        %Calculate residual
        r_n(n) = calcResidual(alpha, beta, v_0, t_n(n));
    
        %Calculate Jacobian
        J_n(n) = calcJacobian(alpha, beta, v_0, t_n(n));
        
        %Calculate the 2-norm of the residual
        norm_2 = sqrt(r_n(n)^2);
    
end

b = t_n(end);

%% True Solutions
t_MAX = b;         %max time [s]
t_analytical = linspace(0,t_MAX,t_MAX*1000); %analytical time array

%Anonymous Functions
%Position
y_true = @(t) (log(beta./alpha .* v_0^2 + 1) + 2 .* log(cos(sqrt(alpha.*beta).*t - atan(sqrt(beta./alpha).*v_0))  ) ) ./ (2.*beta);

%Velocity
v_true = @(t) sqrt(alpha./beta).*tan(atan(sqrt(beta./alpha).*v_0) - sqrt(alpha.*beta).*t);

%Calculate true solns
y_true_plot = y_true(t_analytical);
v_true_plot = v_true(t_analytical);

r_n = [];
J_n = [];

%% Anonymous Functions
f = @(v) -alpha - beta*(v)^2;                              % f(v)
dfdv = @(v) -2*beta*v;                                   % df(v)/dv

%% CMD Window Display
fprintf('Time b for ball B to reach surface:\n');
fprintf('-------------------------------------------------------------------------\n');
fprintf(' t_0:   %2g        b:   %.12f [s]\n', t_0, b);
fprintf('-------------------------------------------------------------------------\n');

%% Numerical Methods to Solve IVPs
results = zeros(length(j_vals),length(n_SI)); % Results array

for i = 1:length(n_SI)                           % Loop over nSI-array
    t = linspace(0,b,n_SI(i)+1);                                      % Nodes array
    y = [];                                      % y-array (approx. y-sol.)
    y(1) = 0;                                   % Initial condition
    h_i = b / n_SI(i);

    

    switch whichMethod                          % Numerical methods
        case 1                                  % Forward Euler method     
            [v_FEuler,method] = computeFEulerSol(v_0,f,n_SI(i),h_i);
            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,v_true_plot, 'b','LineWidth',3);hold on;
            plot(t,v_FEuler,'r','LineWidth',3);
            title("Velocity (" + method + ")");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$v(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(i)))},'Interpreter','latex','location','best');
            y = computePos(v_FEuler,n_SI(i),h_i);

        case 2                                  % Backward Euler method
            [v_BEuler,method] = computeBEulerSol(v_0, EPS, f,dfdv,n_SI(i),h_i);
            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,v_true_plot, 'b','LineWidth',3);hold on;
            plot(t,v_BEuler,'r','LineWidth',3);
            title("Velocity (" + method + ")");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$v(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(i)))},'Interpreter','latex','location','best');
            y = computePos(v_BEuler,n_SI(i),h_i);
        case 3                                  % Trapezoidal method
            [v_trap,method] = computeTrapezoidalSol(v_0, EPS, f,dfdv,n_SI(i),h_i);
            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,v_true_plot, 'b','LineWidth',3);hold on;
            plot(t,v_trap,'r','LineWidth',3);
            title("Velocity (" + method + ")");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$v(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(i)))},'Interpreter','latex','location','best');
            y = computePos(v_trap,n_SI(i),h_i);
        case 4                                  % Heun's method
            [v_heun,method] = computeHeunSol(v_0,f,n_SI(i),h_i);
            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,v_true_plot, 'b','LineWidth',3);hold on;
            plot(t,v_heun,'r','LineWidth',3);
            title("Velocity (" + method + ")");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$v(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(i)))},'Interpreter','latex','location','best');
            y = computePos(v_heun,n_SI(i),h_i);
        case 5                                  % Runge-Kutta method
            [v_RK,method] = computeRungeKuttaSol(v_0,f,n_SI(i),h_i);
            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,v_true_plot, 'b','LineWidth',3);hold on;
            plot(t,v_RK,'r','LineWidth',3);
            title("Velocity (" + method + ")");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$v(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(i)))},'Interpreter','latex','location','best');
            y = computePos(v_RK,n_SI(i),h_i);
        otherwise
            error('Unknown numerical method!');
    end

    tInt = 2 .* j_vals .* h_i;

    %Loop over t-values of interest and store results
    for j = 1:length(tInt)
        results(j,i) = f(tInt(j));

         %Compute (true) errors at t-values of interest
        errors(j,i) = v_true(tInt(j)) - results(j,i);
    end

   
    %Now that the speed v_h(t) is found, find y_h(t) from second IVP:
    if i == 1 || i == length(n_SI)               % First/last element of nSI

        % Store results for plots of first and last element in nSI-array
        
        if i == 1
            t1 = t;                             % Store t-array
            y1 = y;                             % Store y-array
            
        elseif i == length(n_SI)
            t2 = t;                             % Store t-array
            y2 = y;                             % Store y-array
            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,y_true_plot, 'b','LineWidth',3);hold on;
            plot(t2,y2,'r','LineWidth',3);
            title("Position (Trapezoidal Method)");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$y(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(i)))},'Interpreter','latex','location','best');
        end
    end
    


    % CMD Window Display
    fprintf('\nv_h-values at t_int:\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('    h      %2g     %2g      %2g      %2g      %2g\n',tInt');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf(' %.4f  %11.8f  %11.8f  %11.8f  %11.8f  %11.8f\n',...
        [h_i,results(j,:)]);
    
    fprintf('\nTrue errors at t_int:\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('    h      %2g     %2g      %2g      %2g      %2g\n',tInt');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf(' %.4f  %11.8f  %11.8f  %11.8f  %11.8f  %11.8f\n',...
        [h_i,errors(j,:)]);

end

            figure;
            set(gcf,'Position',[0 25 750 800]);     
            plot(t_analytical,y_true_plot, 'b','LineWidth',3);hold on;
            plot(t1,y1,'r','LineWidth',3);
            title("Position (Trapezoidal Method)");                 
            xlabel('$t$','Interpreter','LaTeX');
            xlim([0 b]);
            ylabel('$y(t)$','Interpreter','LaTeX');    
            set(gca,'LineWidth',2,'FontSize',18);
            legend({'Analytical', strcat('$n = $ ',num2str(n_SI(1)))},'Interpreter','latex','location','best');


%% Functions
function [r_tn] = calcResidual(alpha, beta, v_0, t_n)
% Calculates the residual r(x_n) based on current iterative solution t_n.


% Anonymous Functions
    %add negative in front of function for residual
    resid = @(alpha, beta, v_0, t) -(log(beta./alpha * v_0^2 + 1) + 2 .* log(cos(sqrt(alpha.*beta).*t - atan(sqrt(beta./alpha).*v_0))  ) ) ./ (2.*beta);

    r_tn = resid(alpha, beta, v_0, t_n);

end

function [J_tn] = calcJacobian(alpha, beta, v_0, t_n)
% Calculates the Jacobian matrix J(t_n) based on current iterative solution
% t_n.

    %derivative of f(x)
    f_prime = @(alpha, beta, v_0, t) -(sqrt(alpha*beta)*sin(sqrt(alpha*beta)*t-atan(sqrt(beta/alpha)*v_0))) ...
        / (beta*cos(sqrt(alpha*beta)*t-atan(sqrt(beta/alpha)*v_0)) );

    %the jacobian is simply the derivative of y(t) wrt t evaluated at t_n.
    J_tn = f_prime(alpha, beta, v_0, t_n);

end

%% Forward Euler Method Function
function [v,method] = computeFEulerSol(v_0,f,n_SI,h)

% Initialization
method = 'Forward Euler Method';    % Define method for plot

v = [v_0];

for i = 1:n_SI
    v(i+1) = v(i) + h*f(v(i));
end

end

% Backward Euler Method Function
function [v,method] = computeBEulerSol(v0, EPS, f,dfdv,n_SI,h)

% Initialization

v = [];                             % v-array (approx. solution at t_n)
v(1) = v0;                          % Initial condition
IT_MAX = 150;                        % Maximum number of iterations
whichMethod = 2;                    % 1 - Newton's method
                                    % 2 - Predictor-Corrector method

v = [];                             %store solutions
v(1) = v0;

for n = 1:n_SI                          % Loop over all subintervals
    j = 1;                          % Iteration index
    
    v_calc = [];
    if whichMethod == 1
        %Newton's Method
        method = 'Backward Euler//Newton';   % Define method for plot
        
        %initial guess from forward euler
        v_n_0 = v(n) + h*f(v(n));
        v_calc(1) = v_n_0;
        
        %calculate residual and derivative to see if enter loop
        
        r = -(v_calc(1) - v(n) - h * f(v_calc(1)));
        F_prime = 1 - h *dfdv(v_n_0);
        
        norm_2 = sqrt(r^2);
        while norm_2 > EPS && j < IT_MAX
            %calculate change
            Delta_vn = r / F_prime;

            %calculate updated step
            v_calc(j+1) = v_calc(j) + Delta_vn;

            %calculate new residual and derivative
            r = -(v_calc(j+1) - v(n) - h * f(v_calc(j+1)));
            F_prime = 1 - h*dfdv(v_calc(j+1));

            %calculate 2 norm
            norm_2 = sqrt(r^2);
            
            %update count
            j = j + 1;

        end
    else
        %Predictor-Corrector Method
        method = 'Backward Euler//Predictor-Corrector';   % Define method for plot
        
        %initial guess from forward euler
        v_n_0 = v(n) + h*f(v(n));
        v_calc(1) = v_n_0;

        %diff
        diff = abs(v_calc(1) - v(n));
        while diff > EPS && j < IT_MAX

            %calculate the corrector
            v_calc(j+1) = v(n) + h*f(v_calc(j));

            %calculate difference
            diff = abs(v_calc(j+1) - v_calc(j));
            
            %update counter
            j = j + 1;

        end
    end

    % v_n array to store Newton/Predictor-corrector iteration soln
    % for each n_SI.
    v(n+1) = v_calc(end);
end

end

% Trapezoidal Method Function
function [v,method] = computeTrapezoidalSol(v0, EPS, f,dfdv,n_SI,h)

% Initialization
v = [];                             % v-array (approx. solution at t_n)
v(1) = v0;                          % Initial condition
IT_MAX = 150;                        % Maximum number of iterations
whichMethod = 2;                    % 1 - Newton's method
                                    % 2 - Predictor-Corrector method

v = [];                             %store solutions
v(1) = v0;

for n = 1:n_SI                          % Loop over all subintervals
    j = 1;                          % Iteration index
    
    v_calc = [];
    if whichMethod == 1
        %Newton's Method
        method = 'Trapezoidal//Newton';   % Define method for plot
        
        %initial guess from forward euler
        v_n_0 = v(n) + h*f(v(n));

        %calculate trapezoidal method
        v_calc(1) = v_n_0 + h/2*(f(v_n_0) + f(v(n)));
        
        %calculate residual and derivative to see if enter loop
        
        r = -(v_calc(1) - v(n) - h/2 * (f(v_calc(1)) + f(v(n)) ));
        F_prime = 1 - h/2 *dfdv(v_calc(1));
        
        norm_2 = sqrt(r^2);
        while norm_2 > EPS && j < IT_MAX
            %calculate change
            Delta_vn = r / F_prime;

            %calculate updated step
            v_calc(j+1) = v_calc(j) + Delta_vn;

            %calculate new residual and derivative
            r = -(v_calc(j+1) - v(n) - h/2 * (f(v_calc(j+1)) + f(v(n))));
            F_prime = 1 - h/2*dfdv(v_calc(j+1));

            %calculate 2 norm
            norm_2 = sqrt(r^2);
            
            %update count
            j = j + 1;

        end
    else
         %Predictor-Corrector Method
        method = 'Trapezoidal//Predictor-Corrector';   % Define method for plot

        %initial guess from forward euler
        v_n_0 = v(n) + h*f(v(n));
        v_calc(1) = v_n_0;

        %diff
        diff = abs(v_calc(1) - v(n));
        while diff > EPS && j < IT_MAX

            %calculate the corrector
            v_calc(j+1) = v(n) + h/2*( f(v_calc(j)) + f(v(n)));

            %calculate difference
            diff = abs(v_calc(j+1) - v_calc(j));
            
            %update counter
            j = j + 1;

        end
    end

    % Define v_n+1 as last Newton/Predictor-corrector iteration
    v(n+1) = v_calc(end);
end

end


%  Heun's Method Function
function [v,method] = computeHeunSol(v_0,f,n_SI,h)

% Initialization
method = 'Heun''s Method';          % Define method for plot

v = [v_0];

for i = 1:n_SI
    v(i+1) = v(i) + h/2*( f(v(i)) + f(v(i) + h*f(v(i))) );
end

end


% 4th-Order Explicit Runge-Kutta Method Function
function [v,method] = computeRungeKuttaSol(v_0,f,n_SI,h)

% Initialization
method = '4th-Order Explicit Runge-Kutta Method';      % Define method for plot

v = [v_0];

for i = 1:n_SI
    k_1 = f(v(i));
    k_2 = f(v(i) + h/2*k_1);
    k_3 = f(v(i) + h/2*k_2);
    k_4 = f(v(i)+h*k_3);
    v(i+1) = v(i) + h/6 * (k_1 + 2*k_2 + 2*k_3 + k_4);
end


end

function [y] = computePos(v,n_SI,h)
y(1) = 0;

for i = 1:n_SI
    y(i+1) = y(i) + h/2*(v(i) + v(i+1));
end

end
