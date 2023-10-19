% Ex Parcial 2019-20 1Q 
% Problem 1

clearvars
close all

format short e
format compact

clc

a = 0; b = pi;
uA = 0; uB = 8;

nodes = [0; pi/4; pi/2; pi];
elem = [1, 2, 3; 3, 4, 0];

numNodes = size(nodes,1);
numElem = size(elem, 1);

%1st Element
K1 = [7, -8, 1; -8, 16, -8; 1, -8, 7];
F1 = [1; 4; 1];

sys1 = [K1, F1];
fprintf("Local system for element 1:\n")
fprintf("%6.1f%6.1f%6.1f | %4.1f\n", ...
    sys1')

%2nd Element
K2 = [1, -1; -1, 1];
F2 = [3; 3];

sys2 = [K2, F2];
fprintf("\nLocal system for element 2:\n")
fprintf("%6.1f%6.1f | %4.1f\n", ...
    sys2')

%Coupled system
K = zeros(4);
F = zeros(4,1);
Q = zeros(4,1);

e = 1; %couple element 1
%rows = [1;2;3]; cols = rows;
rows = [elem(e,1); elem(e,2); elem(e,3)];
cols = rows;
K(rows, cols) = K(rows,cols) + K1;
F(rows) = F(rows) + F1;

e = 2; %couple element 2
%rows = [3;4]; cols = rows;
rows = [elem(e,1); elem(e,2)];
cols = rows;
K(rows,cols) = K(rows,cols) + K2;
F(rows) = F(rows) + F2;

coupledSys = [K, F];
fprintf("\nCoupled system:\n")
fprintf("%6.1f%6.1f%6.1f%6.1f | %4.1f\n", ...
    coupledSys')

%Boundary conditions
u = zeros(4,1);

fixedNodes = [1,4];
freeNodes = [2,3];

%Natural B.C.
Q(freeNodes) = 0.0;

%Essential B.C.
u(1) = uA; u(4) = uB;

%Reduced system
Fm = F(freeNodes) + Q(freeNodes) - ...
    K(freeNodes, fixedNodes) * u(fixedNodes);
Km = K(freeNodes,freeNodes);
 
reducedSys = [Km Fm];
fprintf("\nReduced system:\n")
fprintf("%6.1f%6.1f | %4.1f\n", ...
    reducedSys')

%Solve the reduced system
um = Km\Fm;

%Add the solution at the free nodes to 
%the solution of the global system
u(freeNodes) = um;

%Post Process
Q = K*u - F;

sols = [(1:numNodes)',nodes,u,Q];
fprintf("\nSolution:\n")
fprintf("%4s%14s%14s%14s\n",...
    'Nod.','x', 'U', 'Q')
fprintf("%4d%14.4e%14.4e%14.4e\n",sols')

%Approximation for the solution u at x = 3*pi/4 
% x = 3pi/4 belongs to the second element, so
xp = 3*pi/4;

%Using the shape functions
nods = elem(2,1:2); %num of nodes of element 2
X = nodes(nods);    %positions of nodes of elem 2

Psi21 = @(t) (t-X(2))/(X(1)-X(2));
Psi22 = @(t) (t-X(1))/(X(2)-X(1));
interpU = u(nods(1))*Psi21(xp) + ...
    u(nods(2))*Psi22(xp);
fprintf('\nInterpolated value of u at x = 3 pi /4:\n')
fprintf('Using shape functions: U_interp = %10.4e\n',interpU)

%Using polyfit
U = u(nods);
p = polyfit(X,U,1);
interpU = polyval(p, xp);

fprintf('Using polyfit: U_interp = %10.4e\n',interpU)

%Approximation for the solution u at x = pi/8 
% x = pi/8 belongs to the first element, so
xp = pi/8;
nods = elem(1,:);
X = nodes(nods);

Psi11 = @(t) ((t-X(2))*(t-X(3)))/((X(1)-X(2))*(X(1)-X(3)));
Psi12 = @(t) ((t-X(1))*(t-X(3)))/((X(2)-X(1))*(X(2)-X(3)));
Psi13 = @(t) ((t-X(1))*(t-X(2)))/((X(3)-X(1))*(X(3)-X(2)));
interpU = u(nods(1))*Psi11(xp) + ...
    u(nods(2))*Psi12(xp) + ...
    u(nods(3))*Psi13(xp);
fprintf('\nInterpolated value of u at x = pi /8:\n')
fprintf('Using shape functions U_interp = %10.4e\n',interpU)

%Using polyfit
U = u(nods);
p = polyfit(X,U,2);
interpU = polyval(p, xp);

fprintf('Using polyfit: U_interp = %10.4e\n',interpU)