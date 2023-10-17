clearvars
close all

format short e
format compact

clc

a = 0; b = pi;
uA = 0; uB = 8;

nodes = [0; pi/4; pi/2; pi];
elem = [1,2,3; 3,4,0];

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