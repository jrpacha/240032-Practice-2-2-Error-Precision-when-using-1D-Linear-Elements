clearvars
close all

a = 0.0; b = 1.0;
a1 = 1.0; a0 = 1.0; f = -1;
U1 = 0.0; UN = 2.0;

% Linear elements and a1(x)= a1, a0(x) = a0, f(x) = f*x 

U = @(x) -x - 6*exp(1)*sinh(x)/(1-exp(2)); %excat solution

numDiv = [5, 50, 500, 5000];

%fprintf('%8s%8s%15s\n','NumElem','h','Error')
fprintf('%8s%9s%14s\n','NumElem','h','Error')
for div = numDiv
    nodes = linspace(a,b,div+1)';
    numNodes = size(nodes,1);
    elem = [(1:numNodes-1)', (2:numNodes)'];
    numElem = size(elem,1);
    h = (b-a)/div;
    Ke = a1*[1,-1; -1, 1]/h + a0*h*[2,1;1,2]/6;

    K = zeros(numNodes);
    F = zeros(numNodes,1);
    Q = zeros(numNodes,1);

    for e = 1:numElem
        rows = [elem(e,1); elem(e,2)];
        cols = rows;
        x1 = nodes(rows(1),1); x2 = nodes(rows(2),1);
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

    err = norm(u-U(nodes(:,1)), inf);
    %fprintf('%7d%14.6e%14.6e\n',div,h,err)
    %fprintf('%7d%14.6e%14.6e\n',div,h,err)
    fprintf('%8d%14.6e%14.6e\n',div,h,err)
end
