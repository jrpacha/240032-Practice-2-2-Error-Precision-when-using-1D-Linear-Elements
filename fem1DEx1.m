clearvars
close all

a = 1.0; b = 2.0;
a1 = 1.0; a0 = 1.0; f = -1;
U1 = 0.0; UN = 2.0;

% Linear elements and a1(x) = a1*x, a0(x) = 1, f(x) = f*x;

div = 10;

nodes = linspace(a,b,div+1)';
numNodes = size(nodes,1);
elem = [(1:numNodes-1)', (2:numNodes)'];
numElem = size(elem,1);
h = (b-a)/div;

K = zeros(numNodes);
F = zeros(numNodes,1);
Q = zeros(numNodes,1);
   
for e = 1:numElem
    rows = [elem(e,1); elem(e,2)];
    cols = rows;
    x1 = nodes(rows(1),1); x2 = nodes(rows(2),1);
    Ke = a1*(x1 + x2)*[1,-1;-1,1]/h/2.0 + a0*h*[2,1;1,2]/6;
    Fe = f*h*[2*x1 + x2; x1+2*x2]/6.0;
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
fprintf('\nMean Value of the approx. soluton: <U> = %.4e\n', meanValue)



