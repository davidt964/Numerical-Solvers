%==========================================================================
% AUTHOR: David L. Tran
%
% MAE 182C, HW 7, POINT COLLOCATION METHOD 
%
% DESCRIPTION: Performs the point collocation method on the double pendulum
% system of differential equations. Plots the approximate angular 
% positions and velocities for both the linear and nonlinear first-order
% initial value problems. An animation is produced that plots the motion of
% the compound pendulum. Finally, the approximate total energy is obtained
% and superimposed with the true total initial energy (a constant). The
% radial basis functions and shape parameter are obtained as well.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%Also, no point in running the script; it doesn't work.

%% Variables/Constants
userChoice = 1;                 % 1-Linear, 0-Nonlinear 

%Pendulum Parameters
m = 3;                          %mass of point masses in [kg] and bt [2,5] kg
L = 3;                          %length of rod in [m] and bt [2,4] m
k = 10;                         %spring constant of rod in [N/m] and bt [20,50] N/m
b = 11;                         %terminating time in [s] and >= 10 [s]
g = 9.81;                       %gravitational constant of Earth in [m/s^2]

x_min = 0;
x_max = 2;

h = 0.01;                       %time step in [s]
t = 0:h:b;                      %time array in [s]

EPS = 0.5;                  %shape parameter
TOL = 10^(-9);                  %Newton's method tolerance
IT_MAX = 100;               %maximum number of Newton iterations

nN = 10;                                    % Number of nodes (= N+1)
nPlot = 1000;                               % Number of x-values for plot.
x = linspace(x_min,x_max,nN);                       % Node-array
yhAtXPlot = zeros(nPlot,1);                 % Radial basis interpolant
Nbar = zeros(4*(nN+1));                           % Matrix of basis functions phi
Bbar = zeros(4*(nN+1));                           % Matrix of derivatives of phi
dfdw = zeros(4*(nN+1));                           % Matrix of derivatives of f

w = zeros(4,nN);                 %weight vector

inc = floor(nPlot/(x_max-x_min)/30);                % Increment for animation
                                            % (assuming 30 fps)

%Initial Angular Positions
theta_1i = 10 * pi/180;            %initial angular position of the first mass in [RAD]
                                    %---NOTE---one initial condition must be
                                    %>= 90 deg for NONLINEAR IVPs and both
                                    %absolute values must be <= 20 deg for
                                    %LINEAR IVPs
theta_2i = 10 * pi/180;            %initial angular position of the second mass in [RAD]
                                    %---NOTE---one initial condition must be
                                    %>= 90 deg for NONLINEAR IVPs and both
                                    %absolute values must be <= 20 deg for
                                    %LINEAR IVPs

%Initial Angular Velocities
omega_1i = 0;                    %initial angular velocity of the first mass in [rad/s]
omega_2i = 0;                    %initial angular velocity of the second mass in [rad/s]


%Initial Energy Components
T_i = 1/2 * m * L^2 * (omega_1i^2 + omega_2i^2);                             %initial kinetic energy
V_ei = -k * L^2 * (cos(theta_2i - theta_2i) - 1);                          %initial elastic potential energy
V_gi = -m * g * L * (cos(theta_2i) + cos(theta_2i));                       %initial gravitational potential energy

%Total Initial Energy
Etotal_th = T_i + V_ei + V_gi;

theta = [theta_1i zeros(numel(t)-1,1)'; theta_2i zeros(length(t)-1,1)'];           %theta array
E_true = Etotal_th .* ones(1,length(t));

T = zeros(1,length(t));
T(1) = T_i;

V_e = zeros(1,length(t));  
V_e(1) = V_ei;

V_g = zeros(1,length(t));  
V_g(1) = V_gi;

E_approx = zeros(1,length(t));                      %numerical energy
E_approx(1) = T(1) + V_e(1) + V_g(1);


%% Anonymous Functions
%Interpolating radial basis functions
phiJAtN = @(n,j) sqrt(1+(EPS*(x(n)-x(j)))^2);           % phi_j (t_n)
dPhiJdtAtN= @(n,j) EPS^2*(x(n)-x(j))/phiJAtN(n,j);   % phi_j'(t_n)

alpha_1 = @(theta_1, theta_2) -g/L * sin(theta_1) - k/m * sin(theta_1 - theta_2);
alpha_2 = @(theta_1, theta_2) -g/L * sin(theta_2) + k/m * sin(theta_1 - theta_2);

alpha_1L = @(theta_1, theta_2) -g/L * theta_1 - k/m * (theta_1 - theta_2);
alpha_2L = @(theta_1, theta_2) -g/L * theta_2 + k/m * (theta_1 - theta_2);


% Define Nbar-, Bbar-, f-, and all dfdw-arrays
for n = 1:nN                                % Loop over all rows n
    for j = 1:nN                            % Loop over all columns j
        Nbar(n,j) = phiJAtN(n,j);           % Store phi_j (x_n)
        Bbar(n,j) = dPhiJdtAtN(n,j);      % Store phi_j'(x_n)
    end

f(n) = Bbar(n,:) .* w(1,:) - Nbar(n,:) .* w(3,:);
f(nN+1) = Bbar(n*(N+1),:) .* w(2,:) - Nbar(n*(N+1),:) .* w(4,:);
f(2*(nN+1)) = g/L .* sin(Nbar(n,:) .* w(1,:)) + k/m .* sin(Nbar(n,:) .* w(1,:) - Nbar(n,:) .* w(2,:));
f(3*(nN+1)) = g/L .* sin(Nbar(n,:) .* w(2,:)) - k/m .* sin(Nbar(n,:) .* w(1,:) - Nbar(n,:) .* w(2,:));
dfdw(n,:) = Bbar(n,:) - Nbar(n,:);
dfdw(n*(N+1),:) = Bbar(n*(N+1),:) - Nbar(n,:);
dfdw(2*(n*N+1),:) = g/L * 1;
dfdw(3*(n*N+1),:) = 1;
end

% Newton's method
r = computeResidual(theta_1i, theta_2i, omega_1i, omega_2i, Bbar, Nbar, w, f, nN);                 % Compute initial residual

norm_2 = norm(r);

while norm_2 > EPS && i < IT_MAX          % Loop u.||r||<=tol or i>=maxIt
    J = computeJacobiMatrix(Bbar, Nbar, w, g, L, k, m);% Compute Jacobi matrix
    DeltaW_i = J\r;                         % Solve lin. sys. for DeltaW_i
    w = w + DeltaW_i;                       % Update w_i+1 (overwrite w)
    i = i+1;                                % Increase iteration counter

    % Update all arrays that depend on w
    for n = 1:nN                            % Loop over all rows n
        f(n) = Bbar(n,:) .* w(1) - Nbar(n,:) .* w(3);
        f(1*(nN+1)) = Bbar(n*(N+1),:) .* w(2) - Nbar(n*(N+1),:) .* w(4);
        f(2*(nN+1)) = g/L .* sin(Nbar(n,:) .* w(1)) + k/m .* sin(Nbar(n,:) .* w(1) - Nbar(n,:) .* w(2));
        f(3*(nN+1)) = g/L .* sin(Nbar(n,:) .* w(2)) - k/m .* sin(Nbar(n,:) .* w(1) - Nbar(n,:) .* w(2));
        dfdw(n,:) = Bbar(n,:) - Nbar(n,:);
        dfdw(n*(N+1),:) = Bbar(n*(N+1),:) - Nbar(n,:);
        dfdw(2*(n*N+1),:) = g/L * 1;
        dfdw(3*(n*N+1),:) = 1;
    end

    r = computeResidual(theta_1i, theta_2i, omega_1i, omega_2i, Bbar, Nbar, w, f, nN);  % Compute next residual

    norm_2 = norm(r);
end

% Note that it is OK to only implement the method for the nonlinear
% problem and then to use the same method for the linear problem
% by simply changing the f- and J-functions

% Construct radial basis interpolants
for i = 1:nPlot                             % Loop ov. all points for plot.
    for j = 1:nN                            % Loop over all nodes
        yhAtXPlot(i) = yhAtXPlot(i) + ...
            w(j) * sqrt(1+(eps*(xPlot(i)-x(j)))^2); % Calculate interpolant
    end
end

% Compute approximate energy E_h(t)



%% Create Animation
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
title('Compound Pendulum, Point Collocation Method');       % Set title
xlabel('x');                                % Set x-label
ylabel('y');                                % Set y-label
xlim([-2-1.2*L 2+1.2*L]);                   % Set x-limits
ylim([-2-1.2*L 2+1.2*L]);                   % Set y-limits
set(gcf,'Position',[50 50 900 900]);        % Change position and size
set(gca,'LineWidth',3,'FontSize',18);       % Change linewidth of axes
axis square;                                % Use same units in x and y
grid on;                                    % Turn on grid

for n = 1:inc:nPlot                         % Loop over all movie frames
    deltaX1 = L*sin(theta1hAtTPlot(n));     % Delta x, pendulum 1 at t_n
    deltaY1 = L*cos(theta1hAtTPlot(n));     % Delta y, pendulum 1 at t_n
    deltaX2 = L*sin(theta2hAtTPlot(n));     % Delta x, pendulum 2 at t_n
    deltaY2 = L*cos(theta2hAtTPlot(n));     % Delta y, pendulum 2 at t_n

    m1x = deltaX1;                          % x-posit. of mass m1 at t_n
    m1y = -deltaY1;                         % y-posit. of mass m1 at t_n
    m2x = deltaX2;                          % x-posit. of mass m2 at t_n
    m2y = -deltaY2;                         % y-posit. of mass m2 at t_n

    % Copy positions of m1 and m2 into plotR1, plotR2, plotM1, and plotM2
    set(plotR1,'xdata',[-2 -2+m1x],'ydata',[0 m1y]);
    set(plotR2,'xdata',[2 2+m2x],'ydata',[0 m2y]);
    set(plotS,'xdata',[-2+m1x 2+m2x],'ydata',[m1y m2y]);
    set(plotM1,'xdata',-2+m1x,'ydata',m1y);
    set(plotM2,'xdata',2+m2x,'ydata',m2y);

    if n == 1
        pause(1);                               % Pause animation before start
    else
        pause(0.005);
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
% Function computeResidual to Compute the Residual
function r = computeResidual(theta_1i, theta_2i, omega_1i, omega_2i, Bbar, Nbar, w, f, nN)


r = -(Bbar.*w - f);                          % Rows 2-nN (1-nN for conven.)

%need to overwrite row 1
r(1) = -(Nbar(1,:).*w(1,:) - theta_1i);                 % Row 1 will be overwritten

%need to overwrite row nN+1
r(nN+1) = -(Nbar(nN+1,:).*w(2,:) - theta_2i);

%need to overwrite row 2nN+1
r(2*(nN+1)) = -(Nbar(2*(nN+1),:).*w(3,:) - omega_1i);

%need to overwrite row 3nN+1
r(3*(nN+1)) = -(Nbar(3*(nN+1),:).*w(4,:) - omega_2i);

end

% Function computeJacobiMatrix to Compute the Jacobi Matrix
function J = computeJacobiMatrix(Bbar, Nbar, w, g, L, k, m)

J = 0;                            % Rows 2-nN (1-nN for conven.)

%need to overwrite row 1
J(1,1) = Bbar(1,:);                        % Row 1 will be overwritten

J(1,3) = -Nbar;                        % Row 1 will be overwritten

%need to overwrite row nN+1
J(nN+1, 2) = Bbar(2,:);

J(nN+1, 4) = -Nbar;


%need to overwrite row 2nN+1
r(2*(nN+1), 1) = g/L * cos(Nbar*w(1,:)) + k/m * sin(Nbar*w(1,:) - Nbar*w(2,:));
r(2*(nN+1), 2) = - k/m * cos(Nbar*w(1,:) - Nbar*w(2,:)) * Nbar;

%need to overwrite row 3nN+1
r(3*(nN+1), 1) = - k/m * cos(Nbar*w(2,:) - Nbar*w(2,:)) * Nbar;
r(3*(nN+1), 2) = g/L * cos(Nbar*w(2,:)) + k/m * sin(Nbar*w(1,:) - Nbar*w(2,:));



end