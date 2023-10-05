clearvars
close all

a = 0.0; b = 1.0;
a1 = 1.0; a0 = 1.0; f = 3.0;
U1 = 0.0; UN = 2.0;

%Quadratic elements and a1(x) = a1, a0(x) = a0, f(x) = f;

div = 4;   %num of elements
order = 2; %order of the elements (linear: order = 1, quadratic: order = 2)

h = (b-a)/div;

%nodes = (a:h/2:b)';
nodes = linspace(a,b,div*order+1)';

numNodes = size(nodes,1);
elem = [(1:2:numNodes-2)', (2:2:numNodes-1)', (3:2:numNodes)'];
numElem = size(elem,1);
h = (b-a)/div;

K = zeros(numNodes);
F = zeros(numNodes,1);
Q = zeros(numNodes,1);

Ke = a1*[7, -8, 1; -8, 16, -8; 1, -8, 7]/h/3.0 + ...
    a0*h*[4, 2, -1; 2, 16, 2; -1, 2, 4]/30.0;

Fe = f*h*[1; 4; 1]/6.0;
   
for e = 1:numElem
    rows = [elem(e,1); elem(e,2); elem(e,3)];
    cols = rows;
    K(rows,cols) = K(rows,cols)+Ke;
    F(rows) = F(rows) + Fe;
end

%Boundary Conditions
fixedNods = [1,numNodes];
freeNods = setdiff(1:numNodes, fixedNods);

% Natural BC.
Q(freeNods) = 0.0;

%Essential BC;
u = zeros(numNodes,1);
u(1) = U1;
u(numNodes) = UN;

%Reduced System
Fm = Q(freeNods) + F(freeNods) - K(freeNods, fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);

%Solve the reduced system
um = Km\Fm;

%Add the solution of the reduced system to the (global) nodal solution
u(freeNods) = um;

%Post Process
Q = F- K*u;

fprintf('%8s%9s%14s%14s\n','NumNod','x','U','Q')
solution = [(1:numNodes)',nodes,u,Q];
fprintf('%8d%14.6e%14.6e%14.6e\n',solution')

meanValue = sum(u)/numNodes;
fprintf('\nMean Value of the approx. soluton: <U> = %.6e\n', meanValue)



