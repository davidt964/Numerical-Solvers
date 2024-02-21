%==========================================================================
% AUTHOR: David L. Tran
%
% 2D Membrane Boundary Value Problem with enforced Dirichlet/essential
% boundary condition with a given load f(x) and linear deflection 
% function u_bar(x).
%
% DESCRIPTION: Solves the 2D membrane BVP such that a PDE of the form 
% -Theta * Delta u(x) = f(x) within a rectangular boundary Omega. The
% Dirichlet/essential boundary condition u(x) = u_bar(x) on Gamma.
% Similarly, a constant tension Theta is assigned to the membrane which is
% a material parameter and f(x) is the vertical load on the membrane and
% u_bar(x) is the given linear deflection function on the boundary Gamma.
% The point collocation method using radial basis functions is employed 
% and n_n collocation points/nodes are used in Omega and on Gamma. We
% seek to find the solution that -Theta * Delta u_h(x_i) = f(x_i) in Omega
% and u_h(x_i) = u_bar(x_i) on Gamma. Visualizations of the contour plot of
% u_h(x) with the nodes x_i along with a surface plot of u_h(x) and the 
% nodes x_i.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables


L = 10;                                  % Length [m]
H = 6;                                   % Height [m]
theta = 17;                              % Tension of the membrane [kN/m]

h = 0.5;                                  % Node distance [m]
nodesX = 0:h:L;                        % Nodes in the x-direction
nodesY = 0:h:H;                        % Nodes in the y-direction
nX = length(nodesX);                        % Number of nodes in the x-direc.
nY = length(nodesY);                        % Number of nodes in the y-direc.
nN = nX*nY;                                   % Number of nodes
nE = nN-1;                              %number of elements (subintervals)

nNOmega = 0;                            % Counter for numb. of n. in Omega
nNGamma = 0;                            % Counter for numb. of n. on Gamma
xOmega = zeros(nN,2);                   % Nodes in Omega
xGamma = zeros(nN,2);                   % Nodes on Gamma

Nbar = zeros(nN,nN);                       % Matrix of basis funct. phi_j(x_i)
BbarXX = zeros(nN,nN);                     % Matrix of deriv. phi_j,xx(x_i)
BbarYY = zeros(nN,nN);                     % Matrix of deriv. phi_j,yy(x_i)
A = zeros(nN,nN);                          % Coefficient matrix A

eps = 1e+0;                             % Shape parameter epsilon

nPlot = 100;                            % Number of nodes for plotting
xPlot = linspace(0,L,nPlot);            % x-values of plot nodes
yPlot = linspace(0,H,nPlot);            % y-values of plot nodes
uhAtXPlot = zeros(nPlot);               % u_h at plot nodes (xPlot,yPlot)
uhAtX = zeros(nN,1);                    % u_h at nodes (x_i,y_i)



% Anonymous functions
ubar = @(x,y) x.*y./60;                        % Prescr. displacement on Gamma [m]
f = @(x,y) 25.*sin(x).*cos(y);                           % Given load in Omega [kN/m^2]
phiJAtI = @(i,j,x) sqrt(1 + (eps * norm(x(i,:)-x(j,:)) )^2 );                  % phi_j(x_i)
phiJAtX = @(j,x,X) sqrt(1 + (eps * norm(x(j,:) - X))^2);                % phi_j(X)

% Second derivatives of phi_j(x) w.r.t x and y evaluated at the nodes x_i
phiJAtIxx = @(i,j,x) eps^2/phiJAtI(i,j,x)^(3/2);
phiJAtIyy = @(i,j,x) eps^2/phiJAtI(i,j,x)^(3/2);


%% Generate Domain and Boundary Nodes
[X,Y] = meshgrid(nodesX,nodesY);        % Generate mesh grid of nodes

for i = 1:nY                            % Loop over all nodes in the y-dir.
    for j = 1:nX                        % Loop over all nodes in the x-dir.
        % Check if node is on the boundary Gamma or in the domain Omega
        if X(i,j) == 0 || X(i,j) == L || Y(i,j) == 0 || Y(i,j) == H
            nNGamma = nNGamma + 1;                  % Increase counter
            xGamma(nNGamma,:) = [X(i,j) Y(i,j)];    % Insert node into arr.
        else
            nNOmega = nNOmega + 1;                  % Increase counter
            xOmega(nNOmega,:) = [X(i,j) Y(i,j)];               % Insert node into arr.
        end
   end
end

xGamma = xGamma(1:nNGamma,:);           % Shrink array of nodes on Gamma
xOmega = xOmega(1:nNOmega,:);           % Shrink array of nodes in Omega
x = [xOmega ; xGamma];                  % Create array of all nodes

%% Construct and Solve Linear System for the Weights
for i = 1:nN                            % Loop over all rows i
    for j = 1:nN                        % Loop over all columns j
        Nbar(i,j) = phiJAtI(i,j,x);         % Store phi_j(x_i) in Nbar
        BbarXX(i,j) = phiJAtIxx(i,j,x);     % Store phi_j,xx(x_i) in BbarXX
        BbarYY(i,j) = phiJAtIyy(i,j,x);     % Store phi_j,yy(x_i) in BbarYY
    end
end

for i = 1:nNOmega                              % Loop over all domain nodes x_i
    A(i,:) = -theta.*(BbarXX(i,:) + BbarYY(i,:));                        % Constr. columns of A
end

for i = (nNOmega+1):nN                  % Loop over all boundary nodes x_i
    A(i,:) = Nbar(i,:);                        % Construct columns of A
end

% Construct the right-hand side vector b = [ f(x_i) ; ubar(x_i) ]
b = [f(x(1:nNOmega,1), x(1:nNOmega ,2)) ; ubar(x((nNOmega+1):nN ,1),x((nNOmega+1):nN ,2))];

w = A\b;                                  % Solve linear system for weights w

%% Construct Interpolant (Numerical Solution)
[XPlot,YPlot] = meshgrid(xPlot,yPlot);  % Generate mesh grid of plot nodes

for j = 1:nN                            % Loop over all nodes
    for i = 1:nPlot                     % Loop ov. all x-val. of plot nodes
        for m = 1:nPlot                 % Loop ov. all y-val. of plot nodes
            % Determine u_h(x)-values at plot nodes x = (xPlot,yPlot)
            uhAtXPlot(m,i) = uhAtXPlot(m,i)+ w(j) * ...
                phiJAtX(j, x,[xPlot(i),yPlot(m)]);
        end
    end

    for i = 1:nN                        % Loop over all nodes
        % Determine u_h(x)-values at nodes x = (x_i,y_i)
        uhAtX(i) = uhAtX(i) + w(j) * phiJAtI(i,j,x);
    end
end

%% Plot Results
figure(1)                               % Open figure 1
contour(XPlot,YPlot,uhAtXPlot,20,'LineWidth',1);    % Contour plot of u_h
hold on;                                % Put hold to on
scatter(x(:,1),x(:,2),20,'filled');     % Plot nodes (x_i,y_i)
title('Approximate Poisson Solution','FontSize',18);    % Set title
xlabel('$x$','Interpreter','latex');                            % Set x-label
ylabel('$y$','Interpreter','latex');                            % Set y-label
set(gcf,'Position',[30 350 850 450]);   % Change position and size
set(gca,'LineWidth',2,'FontSize',18);   % Change linewidth of axes
axis equal;                             % Same units in all directions
a = colorbar;                           % Create color bar
a.Label.String = 'u_h(x,y)';            % Set color bar label

figure(2)                               % Open figure 2
surf(XPlot,YPlot,uhAtXPlot);            % Surface plot of u_h
hold on;                                % Put hold to on
scatter3(x(:,1),x(:,2),uhAtX,50,'filled');  % Plot data p. (x_i,u_h(x_i))
title('Approximate Poisson Solution','FontSize',18);    % Set title
xlabel('$x$','Interpreter','latex');                            % Set x-label
ylabel('$y$','Interpreter','latex');   
set(gcf,'Position',[30 350 850 450]);   % Change position and size
set(gca,'LineWidth',2,'FontSize',18);   % Change linewidth of axes
axis equal;                             % Same units in all directions
a = colorbar;                           % Create color bar
a.Label.String = 'u_h(x,y)';            % Set color bar label

