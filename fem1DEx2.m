clearvars
close all

a = 0.0; b = 1.0;

%Quadratic elements
%order = 2; (linear: order = 1, quadratic: order = 2, cubic: order = 3,...)

%Coefficients
%a1(x)= a1 = ct, a0(x) = a0 = ct, f(x) = f = ct, being,
a1 = 1.0; a0 = 1.0; f = 3.0;

%Essential B.C.
U1 = 0.0; UN = 2.0; %u(a) = U1, u(b) = U2

%Number of elements
div = 4;   

%nodes = linspace(a,b,div*order+1)';
h = (b-a)/div;
nodes = (a:h/2:b)';
numNodes = size(nodes,1);

elem = [(1:2:numNodes-2)', (2:2:numNodes-1)', (3:2:numNodes)'];
numElem = size(elem,1);

K = zeros(numNodes);
F = zeros(numNodes,1);
Q = zeros(numNodes,1);

%local stiffness matrix: the same for all the elements of lenght h
Ke = a1*[7, -8, 1; -8, 16, -8; 1, -8, 7]/h/3.0 + ...
    a0*h*[4, 2, -1; 2, 16, 2; -1, 2, 4]/30.0;  

%local load vector: the same for all the elements of lenght h
Fe = f*h*[1; 4; 1]/6.0;                        
   
for e = 1:numElem
    rows = [elem(e,1); elem(e,2); elem(e,3)];
    cols = rows;

    % Now: 
    % rows(1) holds the number of the global node corresponding to 
    %         the 1st local node of element e
    % rows(2) holds the number of the global node corresponding to 
    %         the 2nd local node of element e 
    % rows(3) holds the number of the global node corresponding to 
    %         the 3rd local node of element e 
    %
    % Recall that elemet e is a **quadratic** element, so each element has
    % **three** nodes

    %The **local** stiffness matrix of element e must be added to the rows
    %(rows(1), rows(2), rows(3)) and cols (cols(1), cols(2), cols(3)) of the
    %**global** stiffness matrix, so
    K(rows,cols) = K(rows,cols)+Ke;

    %The **local** load vector of element e must be added to the rows
    %(rows(1), rows(2)) of the **global** load vector, so
    F(rows) = F(rows) + Fe;
end

%Boundary Conditions
fixedNods = [1,numNodes];
freeNods = setdiff(1:numNodes, fixedNods);

%Natural B.C.
Q(freeNods) = 0.0;

%Essential B.C.
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

fprintf('%34s\n','Exercise 1')
fprintf('%8s%9s%14s%14s\n','NumNod','x','U','Q') %print table's header
solution = [(1:numNodes)',nodes,u,Q];
fprintf('%8d%14.6e%14.6e%14.6e\n',solution')     %print table's entries                                         

meanValue = sum(u)/numNodes;
fprintf('\nMean Value of the approx. soluton: <U> = %.6e\n', meanValue)
