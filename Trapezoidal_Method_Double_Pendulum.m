%==========================================================================
% AUTHOR: David L. Tran
%
% TRAPEZOIDAL METHOD 
%
% DESCRIPTION: Performs the trapezoidal method on the double pendulum
% system of differential equations. Plots the approximate angular 
% positions and velocities for both the linear and nonlinear first-order
% initial value problems. An animation is produced that plots the motion of
% the compound pendulum. Finally, the approximate total energy is obtained
% and superimposed with the true total initial energy (a constant).
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables/Constants
userChoice = 1;                 % 1-NonLinear, 0-Linear 

%Pendulum Parameters
m = 3;                          %mass of point masses in [kg] and bt [2,5] kg
L = 3;                          %length of rod in [m] and bt [2,4] m
k = 30;                         %spring constant of rod in [N/m] and bt [20,50] N/m
b = 15;                         %terminating time in [s] and >= 10 [s]
g = 9.81;                       %gravitational constant of Earth in [m/s^2]

h = 0.005;                      %time step in [s]
t = 0:h:b;                      %time array in [s]

EPS = 10^(-12);                 %tolerance value
IT_MAX = 1000;                  %maximum number of iterations for while loop

%Initial Angular Positions
theta_01 = 90 * pi/180;             %initial angular position of the first mass in [deg]
                                    %---NOTE---one initial condition must be
                                    %>= 90 deg for NONLINEAR IVPs and both
                                    %absolute values must be <= 20 deg for
                                    %LINEAR IVPs
theta_02 = 20 * pi/180;             %initial angular position of the second mass in [deg]
                                    %---NOTE---one initial condition must be
                                    %>= 90 deg for NONLINEAR IVPs and both
                                    %absolute values must be <= 20 deg for
                                    %LINEAR IVPs

%Initial Angular Velocities
omega_1 = 0;                    %initial angular velocity of the first mass in [rad/s]
omega_2 = 0;                    %initial angular velocity of the second mass in [rad/s]

%Initial Energy Components
T_i = 1/2 * m * L^2 * (omega_1^2 + omega_2^2);                             %initial kinetic energy
V_ei = -k * L^2 * (cos(theta_01 - theta_02) - 1);                          %initial elastic potential energy
V_gi = -m * g * L * (cos(theta_01) + cos(theta_02));                       %initial gravitational potential energy

%Total Initial Energy
Etotal_th = T_i + V_ei + V_gi;

theta = [theta_01 zeros(numel(t)-1,1)'; theta_02 zeros(length(t)-1,1)'];           %theta array
E_true = Etotal_th .* ones(1,length(t));

T = zeros(1,length(t));
T(1) = T_i;

V_e = zeros(1,length(t));  
V_e(1) = V_ei;

V_g = zeros(1,length(t));  
V_g(1) = V_gi;

E_approx = zeros(1,length(t));                      %numerical energy
E_approx(1) = T(1) + V_e(1) + V_g(1);

%% Open Figure 1
dummyX = zeros(2,1);                        % Dummy array for plotting
dummyM = 0;                                 % Dummy scalar for plotting

figure(1);                                  % Open figure 1
plotR1 = plot(dummyX,dummyX,'LineWidth',4); % Plot rod 1
hold on;                                    % Put hold to on
plotR2 = plot(dummyX,dummyX,'LineWidth',4); % Plot rod 2
plotS = plot(dummyX,dummyX,'LineWidth',2);  % Plot rubber band / spring
plotM1 = plot(dummyM,dummyM,'.','MarkerSize',sqrt(m)*50);   % Plot mass 1
plotM2 = plot(dummyM,dummyM,'.','MarkerSize',sqrt(m)*50);   % Plot mass 2
plot(-2,0,'.','MarkerSize',30);             % Plot pin 1
plot(2,0,'.','MarkerSize',30);              % Plot pin 2
title('Compound Pendulum, Trapezoidal Method');             % Set title
xlabel('$x$', 'Interpreter','latex');                                % Set x-label
ylabel('$y$','Interpreter','latex');                                % Set y-label
xlim([-2-1.2*L 2+1.2*L]);                   % Set x-limits
ylim([-2-1.2*L 2+1.2*L]);                   % Set y-limits
set(gcf,'Position',[50 50 900 900]);        % Change position and size
set(gca,'LineWidth',3,'FontSize',18);       % Change linewidth of axes
axis square;                                % Use same units in x and y
grid on;                                    % Turn on grid


%% Anonymous Functions
%Differential Equations
alpha_1 = @(theta_1, theta_2) -g/L * sin(theta_1) - k/m * sin(theta_1 - theta_2);
alpha_2 = @(theta_1, theta_2) -g/L * sin(theta_2) + k/m * sin(theta_1 - theta_2);

alpha_1L = @(theta_1, theta_2) -g/L * theta_1 - k/m * (theta_1 - theta_2);
alpha_2L = @(theta_1, theta_2) -g/L * theta_2 + k/m * (theta_1 - theta_2);

omega_1t = [omega_1];
omega_2t = [omega_2];

omega = [omega_1t;omega_2t];


%% Compute Trapezoidal Solution and Total Energy and Create Animation on-the-fly
for n = 1:length(t)-1                                  % Loop over all subintervals
%     ..                                      % Trapezoidal method
%     % Note that it is OK to only implement the method for the nonlinear
%     % problem and then to use the same method for the linear problem
%     % by simply changing the f- and Jf-functions
%     
%     ..

    %restart counter
    j = 1;

    %initial guess from forward euler
    if userChoice == 1
        omega(1,n+1) = omega(1,n) + h * alpha_1(theta(1,n), theta(2,n));
        omega(2,n+1) = omega(2,n) + h * alpha_2(theta(1,n), theta(2,n));
    elseif userChoice == 0
        omega(1,n+1) = omega(1,n) + h * alpha_1L(theta(1,n), theta(2,n));
        omega(2,n+1) = omega(2,n) + h * alpha_2L(theta(1,n), theta(2,n));
    end

    theta(1,n+1) = theta(1,n) + h/2 * (omega(1,n) + omega(1,n+1));
    theta(2,n+1) = theta(2,n) + h/2 * (omega(2,n) + omega(2,n+1));

    omega_calc = [];
    theta_calc = [];

    omega_calc(1,1) = omega(1,n+1);
    omega_calc(2,1) = omega(2,n+1);

    theta_calc(1,1) = theta(1,n+1);
    theta_calc(2,1) = theta(2,n+1);

    %Calculate residual
    resid(1,1) = -(theta(1,n+1) - theta(1,n) - h/2 * (omega_calc(1,1) + omega(1,n)));
    resid(2,1) = -(theta(2,n+1) - theta(2,n) - h/2 * (omega_calc(2,1) + omega(2,n)));  


    %calculate jacobian
    J_mat = eye(2) - h/2*Jf(m, g, L, k, theta(1,n+1), theta(2,n+1), userChoice);

    %calculate 2-norm to see if enter loop
    norm_2 = sqrt(resid(1,1)^2 + resid(2,1)^2);

    while norm_2 > EPS && j < IT_MAX

        %calculate change
        Delta_omega = J_mat\resid;

        %update step
        omega_calc(1,j+1) = omega_calc(1,j) + Delta_omega(1);
        omega_calc(2,j+1) = omega_calc(1,j) + Delta_omega(2);

        theta_calc(1,j+1) = theta(1,n) + h/2 * (omega_calc(1,j+1) + omega(1,n));
        theta_calc(2,j+1) = theta(2,n) + h/2 * (omega_calc(2,j+1) + omega(2,n));

        %calculate updated residual with new step
        resid(1,1) = -(theta_calc(1,j+1) - theta(1,n) - h/2 * (omega_calc(1,j+1) + omega(1,n) ));
        resid(2,1) = -(theta_calc(2,j+1) - theta(2,n) - h/2 * (omega_calc(2,j+1) + omega(2,n) ));

        %calculate Jacobian with new step
        J_mat =  eye(2) - h/2*Jf(m, g, L, k, theta_calc(1,j+1), theta_calc(2,j+1), userChoice);

        %update counter
        j = j + 1;

        %calculate updated 2-norm
        norm_2 = sqrt(resid(1,1)^2 + resid(2,1)^2);
    end

    
    %Trapezoidal method for both omega and theta
    if userChoice == 1
        %Nonlinear
        omega(1,n+1) = omega(1,n) + h/2 * (alpha_1(theta(1,n), theta(2,n)) + alpha_1(theta_calc(1,end), theta_calc(2,end)));
        omega(2,n+1) = omega(2,n) + h/2 * (alpha_2(theta(1,n), theta(2,n)) + alpha_2(theta_calc(1,end), theta_calc(2,end)));
    elseif userChoice == 0
        %Linear
        omega(1,n+1) = omega(1,n) + h/2 * (alpha_1L(theta(1,n), theta(2,n)) + alpha_1L(theta_calc(1,end), theta_calc(2,end)));
        omega(2,n+1) = omega(2,n) + h/2 * (alpha_2L(theta(1,n), theta(2,n)) + alpha_2L(theta_calc(1,end), theta_calc(2,end)));
    end

    theta(1,n+1) = theta(1,n) + h/2 * (omega(1,n+1) + omega(1,n));
    theta(2,n+1) = theta(2,n) + h/2 * (omega(2,n+1) + omega(2,n));

    
    deltaX1 = L*sin(theta(1,n+1));          % Delta x, pendulum 1 at t_n+1
    deltaY1 = L*cos(theta(1,n+1));          % Delta y, pendulum 1 at t_n+1
    deltaX2 = L*sin(theta(2,n+1));          % Delta x, pendulum 2 at t_n+1
    deltaY2 = L*cos(theta(2,n+1));          % Delta y, pendulum 2 at t_n+1

    m1x = deltaX1;                          % x-posit. of mass m1 at t_n+1
    m1y = -deltaY1;                         % y-posit. of mass m1 at t_n+1
    m2x = deltaX2;                          % x-posit. of mass m2 at t_n+1
    m2y = -deltaY2;                         % y-posit. of mass m2 at t_n+1

    %Calculate energies
    T(n+1) = 1/2 * m * L^2 * (omega(1,n+1)^2 + omega(2,n+1)^2);
    V_e(n+1) = -k * L^2 * (cos(theta(1,n+1) - theta(2,n+1)) - 1);
    V_g(n+1) = -m * g * L * (cos(theta(1,n+1)) + cos(theta(2,n+1)));

    E_approx(n+1) = T(n+1) + V_e(n+1) + V_g(n+1);

    % Copy positions of m1 and m2 into plotR1, plotR2, plotM1, and plotM2
    set(plotR1,'xdata',[-2 -2+m1x],'ydata',[0 m1y]);
    set(plotR2,'xdata',[2 2+m2x],'ydata',[0 m2y]);
    set(plotS,'xdata',[-2+m1x 2+m2x],'ydata',[m1y m2y]);
    set(plotM1,'xdata',-2+m1x,'ydata',m1y);
    set(plotM2,'xdata',2+m2x,'ydata',m2y);

    if n == 1
        pause(1);                              % Pause animation before start
    else
        pause(h);
    end

   drawnow;                                 % Update plot
end




%% Final Plots
%Angular Position
figure;
set(gcf,'Position',[0 25 750 800]);
plot(t, theta(1,:),'b', 'LineWidth',3); hold on;
plot(t, theta(2,:), 'r','LineWidth',3);
hold on;
title('Trapezoidal Approximations $\theta_{h,1}(t)$ and $\theta_{h,2}(t)$',...
    'Interpreter','LaTeX','FontSize',18);   
xlabel('$t$ (s)','Interpreter','LaTeX','FontSize',18);  
ylabel('$\theta_{h,1}$ and $\theta_{h,2}$ (rad)','Interpreter','LaTeX',...
    'FontSize',18);                        
legend('$\theta_{h,1}(t)$','$\theta_{h,2}(t)$','Interpreter','LaTeX',...
    'Location','best');    

%Angular Velocity
figure;
set(gcf,'Position',[0 25 750 800]);
plot(t, omega(1,:),'b', 'LineWidth',3);hold on;
plot(t, omega(2,:), 'r','LineWidth',3);
hold on;  
title('Trapezoidal Approximations $\omega_{h,1}(t)$ and $\omega_{h,2}(t)$',...
    'Interpreter','LaTeX','FontSize',18);   
xlabel('$t$ (s)','Interpreter','LaTeX','FontSize',18); 
ylabel('$\omega_{h,1}$ and $\omega_{h,2}$ (rad/s)','Interpreter','LaTeX',...
    'FontSize',18);                         
legend('$\omega_{h,1}(t)$','$\omega_{h,2}(t)$','Interpreter','LaTeX',...
    'Location','best');                

%Total Energy
figure;
set(gcf,'Position',[0 25 750 800]); hold on;
plot(t, E_true, 'b','LineWidth',3);
plot(t, E_approx, 'r','LineWidth',3)
% plot(t, E_h, 'r','LineWidth',3);
title('Energy Approximation $E_{h}(t)$','Interpreter','LaTeX',...
    'FontSize',18);                         % Set title
xlabel('$t$ (s)','Interpreter','LaTeX','FontSize',18);      % Set x-label
ylabel('$E$ (J)','Interpreter','LaTeX','FontSize',18);  % Set y-label
set(gcf,'Position',[30 350 1000 600]);      % Change position and size
set(gca,'LineWidth',3,'FontSize',18);       % Change linewidth of axes
legend('$E(t)$','$E_h(t)$','Interpreter','LaTeX','Location',...
    'best');                           % Define legend  


%% Functions
% Function to Compute the Right-hand Side Function f
function fReturn = f(alpha_1L, alpha_2L, theta_1, theta_2)

fReturn(1) = alpha_1L(theta_1, theta_2);
fReturn(2) = alpha_2L(theta_1, theta_2);

end

% Function to Compute the Jacobi Matrix J_f
function JfReturn = Jf(m, g, L, k, theta_1, theta_2, userChoice)

if userChoice == 1
    % dw_dot_1/dtheta_1
    JfReturn(1,1) = -g/L * cos(theta_1) - k/m * cos(theta_1 - theta_2);
    
    % dw_dot_1/dtheta_2
    JfReturn(1,2) = k/m * cos(theta_1 - theta_2);
    
    % dw_dot_2/dtheta_1
    JfReturn(2,1) = k/m * cos(theta_1 - theta_2);
    
    % dw_dot_2/dtheta_2
    JfReturn(2,2) = -g/L * cos(theta_2) - k/m * cos(theta_1 - theta_2);

elseif userChoice == 0
    JfReturn(1,1) = -(g/L) - (k/m);
    
    % dw_dot_1/dtheta_2
    JfReturn(1,2) = k/m;
    
    % dw_dot_2/dtheta_1
    JfReturn(2,1) = k/m;
    
    % dw_dot_2/dtheta_2
    JfReturn(2,2) = -(g/L) - (k/m);

end


end
