%==========================================================================
% AUTHOR: David L. Tran
%
% Compound Double Pendulum, Finite Element Method/Discrete Galerkin (0) & 
% (1) Methods
%
% DESCRIPTION: Solves the compound double pendulum problem using either
% a linear governing system (for angles small enough) or a nonlinear 
% governing system of differential equations (for larger angles). The
% results are computed using the discrete Galerkin (0) or (1) method.
% An animation is generated of the motion of the pendulums along with the
% computed energy compared with the theoretical constant conserved energy,
% along with the angular positions and velocities wrt t.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables
lin = 1;                                    % 1-Linear pr., 0-Nonlinear pr.
dG = 1;                                     % 0 - dG(0), 1 - dG(1)
m = 3;                          %mass of point masses in [kg] and bt [2,5] kg
L = 3;                          %length of rod in [m] and bt [2,4] m
k = 30;                         %spring constant of rod in [N/m] and bt [20,50] N/m
g = 9.81;                       %gravitational constant of Earth in [m/s^2]

TOL = 10^(-12);                 %predictor-corrector tolerance
MAX_IT = 150;                   % Maximum number of iterations
eps = 1.0e+00;                  %shape parameter

nN = 1000;                                   % Number of nodes
nE = nN-1;                                  % Number of elements (subint.)
a = 0;                          %starting time in [s]
b = 11;                         %terminating time in [s] and >= 10 [s]
h = (b-a)/nE;

theta_10 = 20 * pi / 180;                                % Initial condition for theta_1 in [rad]
theta_20 = 15 * pi / 180;                                % Initial condition for theta_2 in [rad]
omega_10 = 0;                                % Initial condition for omega_1
omega_20 = 0;                                % Initial condition for omega_2

t = linspace(0, b, nN);                           % Node-array
thetaPlus = zeros(4,nE);                    % theta^+_h, omega^+_h at t_n
thetaMinus = zeros(4,nE);                   % theta^-_h, omega^-_h at t_n
thetaMinus(1,1) = theta_10;                  % Initial theta^-_h(t_0) value
thetaMinus(2,1) = theta_20;                  % Initial theta^-_h(t_0) value
thetaMinus(3,1) = omega_10;                  % Initial omega^-_h(t_0) value
thetaMinus(4,1) = omega_20;                  % Initial omega^-_h(t_0) value


%Initial Energy Components
T_i = 1/2 * m * L^2 * (omega_10^2 + omega_20^2);                             %initial kinetic energy
V_ei = -k * L^2 * (cos(theta_10 - theta_20) - 1);                          %initial elastic potential energy
V_gi = -m * g * L * (cos(theta_10) + cos(theta_20));                       %initial gravitational potential energy

%Total Initial Energy
Etotal_th = T_i + V_ei + V_gi;

theta = [theta_10 zeros(numel(t)-1,1)'; theta_20 zeros(length(t)-1,1)'];           %theta array
E_true = Etotal_th .* ones(1,length(t));

T = zeros(1,length(t));
T(1) = T_i;

V_e = zeros(1,length(t));  
V_e(1) = V_ei;

V_g = zeros(1,length(t));  
V_g(1) = V_gi;

E_approx = zeros(1,length(t));                      %numerical energy
E_approx(1) = T(1) + V_e(1) + V_g(1);


I = eye(4);                                 % 4x4 identity matrix
K_n = 1/2.*[I, I; -I, I];                             % (Element) stima for dG(1)
color = lines(6);                           % Default Matlab colors
dummyX = zeros(2,1);                        % Dummy array for plotting
dummyM = 0;                                 % Dummy scalar for plotting

%% Open Figure 1
figure(1);                                  % Open figure 1
plotR1 = plot(dummyX,dummyX,'LineWidth',4); % Plot rod 1
hold on;                                    % Put hold to on
plotR2 = plot(dummyX,dummyX,'LineWidth',4); % Plot rod 2
plotS = plot(dummyX,dummyX,'LineWidth',2);  % Plot rubber band / spring
plotM1 = plot(dummyM,dummyM,'.','MarkerSize',sqrt(m)*50);   % Plot mass 1
plotM2 = plot(dummyM,dummyM,'.','MarkerSize',sqrt(m)*50);   % Plot mass 2
plot(-2,0,'.','MarkerSize',30);             % Plot pin 1
plot(2,0,'.','MarkerSize',30);              % Plot pin 2
title('Compound Pendulum, Finite Element Method');          % Set title
xlabel('$x$','Interpreter','latex');                                % Set x-label
ylabel('$y$','Interpreter','latex');                                % Set y-label
xlim([-2-1.2*L 2+1.2*L]);                   % Set x-limits
ylim([-2-1.2*L 2+1.2*L]);                   % Set y-limits
set(gcf,'Position',[50 50 900 900]);        % Change position and size
set(gca,'LineWidth',3,'FontSize',18);       % Change linewidth of axes
axis square;                                % Use same units in x and y
grid on;


%% Compute dG(0)- or dG(1)-solution and Create Animation
for n = 1:nE                                  % Loop over all elements
    j = 1;                                  % Iteration counter
    if dG == 0                              % dG(0)-method

        thetaMinusIter = [];       % Initialize iterations array
        thetaMinusIter(:,j) = thetaMinus(:,n) +  h.*f(thetaMinus(:,n),g,L,k,m,lin);  % Compute pred. (Fwd Euler)
        thetaMinusIter(:,j+1) = thetaMinus(:,n) + h.*f(thetaMinusIter(:,j),g,L,k,m,lin);          % Compute init. correc.

        % Loop until || theta(j+1) - theta(j) || < eps or j >= maxIt
        while norm(thetaMinusIter(:,j+1) - thetaMinusIter(:,j)) > TOL && j < MAX_IT
            j = j + 1;                     % Increase iteration counter j
            thetaMinusIter(:,j+1) = thetaMinus(:,n) + h/2*(f(thetaMinusIter(:,j),g,L,k,m,lin) + f(thetaMinus(:,n),g,L,k,m,lin));       % Comp. next corr.
        end

        if j == MAX_IT                     % Error if maxIter is reached
            error('Maximum number of iterations reached!');
        end

        thetaMinus(:,n+1) = thetaMinusIter(:,j);              % Store theta^- results
        thetaPlus(:,n) = thetaMinusIter(:,j);                 % Store theta^+ results

    elseif dG == 1                                    % dG(1)-method
        j = 1;

        thetaPlusIter = zeros(4,1);        % Initialize plus iter. array
        thetaMinusIter = zeros(4,1);       % Initialize minus iter. array

        % Compute predictor on plus side at t_n from minus side at t_n
        thetaPlusIter(:,j) = thetaPlus(:,n);
        % Compute predictor on minus side at t_n+1 from forward Euler
        thetaMinusIter(:,j) = thetaMinus(:,n) + ...
            h*f(thetaPlusIter(:,j),g,L,k,m,lin);

        % Define (element) load vector f_n for initial corrector
        f_n = [thetaMinus(:,n); zeros(4,1)] + ...
           h/2.*[f(thetaPlusIter(:,j),g,L,k,m,lin); f(thetaMinusIter(:,j),g,L,k,m,lin)];

        c = K_n \ f_n;                               % Solve lin. s. for init. corr.
        thetaPlusIter(:,j+1) = c(1:4,:);           % Define theta^+_n(j+1)
        thetaMinusIter(:,j+1) = c(5:8,:);          % Define theta^-_n+1(j+1)
    
        % Loop until || theta(j+1) - theta(j) || < eps or j >= maxIt
        while norm(c - [thetaPlusIter(:,j) ; thetaMinusIter(:,j)]) > TOL && j < MAX_IT
            j = j + 1;                      % Increase iteration counter j

            % Redefine (element) load vector f_n for next corrector
            f_n = [thetaMinus(:,n); zeros(4,1)] + ...
                 h/2.*[f(thetaPlusIter(:,j),g,L,k,m,lin);
                 f(thetaMinusIter(:,j),g,L,k,m,lin)];
%             
           
            c = K_n\f_n;                         % Solve lin. s. for next corr.

            thetaPlusIter(:,j+1) =  c(1:4,:);       % Redefine theta^+_n(j+1)
            thetaMinusIter(:,j+1) = c(5:8,:);      % Redefine theta^-_n+1(j+1)

        end

        if j == MAX_IT                     % Error if maxIter is reached
            error('Maximum number of iterations reached!');
        end

    thetaPlus(:,n) = thetaPlusIter(:,j+1);                 % Store theta^+ results
    thetaMinus(:,n+1) = thetaMinusIter(:,j+1);              % Store theta^- results

    if n == 1
        pause(1);                               % Pause animation before start
    else
        pause(h/2);
    end
    end 
       %Calculate energy components
    T(n+1) = 1/2 * m * L^2 * (thetaMinus(3,n+1)^2 + thetaMinus(4,n+1)^2);
    V_e(n+1) = -k .* L^2 .* (cos(thetaMinus(1,n+1) - thetaMinus(2,n+1)) - 1);
    V_g(n+1) = -m .* g .* L .* (cos(thetaMinus(1,n+1)) + cos(thetaMinus(2,n+1)));

    deltaX1 = L*sin(thetaMinus(1,n));       % Delta x, pendulum 1 at t_n
    deltaY1 = L*cos(thetaMinus(1,n));       % Delta y, pendulum 1 at t_n
    deltaX2 = L*sin(thetaMinus(2,n));       % Delta x, pendulum 2 at t_n
    deltaY2 = L*cos(thetaMinus(2,n));       % Delta y, pendulum 2 at t_n

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
   drawnow;                                 % Update plot
end

% Compute approximate (smoothed) energy
Eh = T + V_e + V_g;

%% Plot Results
figure(2);                                  % Open figure 2
subplot1 = subplot(2,1,1);                  % Create subplot 1
hold on;                                    % Put hold to on
for n = 1:nE                                % Loop over all elements
    plot([t(n) t(n+1)],[thetaPlus(1,n) thetaMinus(1,n+1)],'LineWidth',3,...
        'Color',color(1,:));                % Plot theta_{h,1}
    plot([t(n) t(n+1)],[thetaPlus(2,n) thetaMinus(2,n+1)],'LineWidth',3,...
        'Color',color(2,:));                % Plot theta_{h,2}
end
title('Finite Element Approximations $\theta_{h,1}(t)$ and $\theta_{h,2}(t)$',...
    'Interpreter','LaTeX','FontSize',18);   % Set title
xlabel('$t$ [s]','Interpreter','LaTeX','FontSize',18);  % Set x-label
ylabel('$\theta_{h,1}$ and $\theta_{h,2}$ [rad]','Interpreter','LaTeX',...
    'FontSize',18);                         % Set y-label
legend('$\theta_{h,1}(t)$','$\theta_{h,2}(t)$','Interpreter','LaTeX',...
    'Location','NorthEast');                % Define legend

subplot2 = subplot(2,1,2);                  % Create subplot 2
hold on;                                    % Put hold to on
for n = 1:nE                                % Loop over all elements
    plot([t(n) t(n+1)],[thetaPlus(3,n) thetaMinus(3,n+1)],'LineWidth',3,...
        'Color',color(1,:));                % Plot omega_{h,1}
    plot([t(n) t(n+1)],[thetaPlus(4,n) thetaMinus(4,n+1)],'LineWidth',3,...
        'Color',color(2,:));                % Plot omega_{h,2}
end
title('Finite Element Approximations $\omega_{h,1}(t)$ and $\omega_{h,2}(t)$',...
    'Interpreter','LaTeX','FontSize',18);   % Set title
xlabel('$t$ [s]','Interpreter','LaTeX','FontSize',18);  % Set x-label
ylabel('$\omega_{h,1}$ and $\omega_{h,2}$ [rad/s]','Interpreter','LaTeX',...
    'FontSize',18);                         % Set y-label
legend('$\omega_{h,1}(t)$','$\omega_{h,2}(t)$','Interpreter','LaTeX',...
    'Location','NorthEast');                % Define legend
set(gcf,'Position',[50 50 900 900]);        % Change position and size

figure(3);                                  % Open figure 3
plot(t,E_true,'LineWidth',2);  % Plot true energy E
hold on;                                    % Put hold to on
plot(t,Eh,'LineWidth',2);                   % Plot approximate energy E_h
title('Energy Approximation $E_{h}(t)$','Interpreter','LaTeX',...
    'FontSize',18);                         % Set title
xlabel('$t$ [s]','Interpreter','LaTeX','FontSize',18);      % Set x-label
ylabel('$E_{h}$ [J]','Interpreter','LaTeX','FontSize',18);  % Set y-label
set(gcf,'Position',[30 350 1000 600]);      % Change position and size
set(gca,'LineWidth',3,'FontSize',18);       % Change linewidth of axes
legend('$E(t)$','$E_h(t)$','Interpreter','LaTeX','Location',...
    'NorthEast');                           % Define legend
grid on;                                    % Turn on grid

%% Function to Compute the Right-hand Side Function f
function fReturn = f(theta,g,L,k,m,lin)

% Initialization
theta1 = theta(1,:);                         % Define theta_1
theta2 = theta(2,:);                         % Define theta_2
omega1 = theta(3,:);                         % Define omega_1
omega2 = theta(4,:);                         % Define omega_2

fReturn = [omega1;                              % Define f_1
          omega2];                             % Define f_2

if lin == 1                                 % Linear problem
    % Define f_3 and f_4 based on linear compound pendulum problem
    fReturn(3:4) = [-g./L .* theta1 - k./m .* (theta1 - theta2);
                    -g./L .* theta2 + k./m .* (theta1 - theta2)];
elseif lin == 0                                        % Nonlinear problem
    % Define f_3 and f_4 based on nonlinear compound pendulum problem
    fReturn(3:4) = [-g./L .* sin(theta1) - k./m .* sin(theta1 - theta2);
                    -g./L .* sin(theta2) + k./m .* sin(theta1 - theta2)];
    
end

end
