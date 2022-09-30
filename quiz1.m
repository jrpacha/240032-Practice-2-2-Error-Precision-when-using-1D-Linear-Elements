% Practice 2.2 Quiz 1.
% Consider the Boundary Value Problem (BVP)
%
% -((x+3)u')' + u = x, 0 < x < L,
%           u'(0) = u'_0,
%            u(L) = u_L,
%
% By the FEM, taking N = 100 linear elements, compute the nodal soluiton,
% u(k) k=1,...,N+1, of the BVP for L=50, u'_0 = 4.94, u_L = 1 and give
% its maximum value, max u(k).
%
% Hint: The value of the FEM approximation for u at node 51 is 2.4895614e+01 
clearvars
close all

numDiv=100; %Number of divisions
a=0; b=50;
du0=4.94; uL=1.0;    % Boundary conditions 
alpha=1.0; beta=3.0; % a1(x) = alpha*x + beta
gamma=1.0;           % a0(x) = gamma

%Geometry
h=(b-a)/numDiv;
nodes=(a:h:b)';
elem=[1:numDiv;2:numDiv+1]';

numNodes=size(nodes,1);
numElem=size(elem,1);

%Assembly
K=zeros(numNodes);
F=zeros(numNodes,1);
Q=zeros(numNodes,1);

Ke0 = beta*[1,-1;-1,1]/h + gamma*h*[2,1;1,2]/6;

for e=1:numElem
  rows=[elem(e,1), elem(e,2)];
  cols=rows;
  x1=rows(1,1);
  x2=rows(1,2);
  Ke=0.5*alpha*(x1+x2)*[1,-1;-1,1]/h + Ke0;
  Fe=h*[2*x1+x2;x1+2*x2]/6;
  K(rows,cols)=K(rows,cols)+Ke;
  F(rows)=F(rows)+Fe;
end

%B.C.
fixedNods=numNodes;
freeNods=setdiff(1:numNodes,fixedNods);
u=zeros(numNodes,1);

%Natural B.C.
Q(1)=-alpha*du0;
Q(2:numNodes-1)=0; %Not necessary: Q is initialised to 0

%Essential B.C.
u(fixedNods)=uL;

%Reduced system
Qm=F(freeNods)+Q(freeNods)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

%Post Process
Q = K*u - F;

%Print results
for i=1:numNodes
  if i==1
    fprintf('\n%35s\n%3s%10s%18s%16s\n','NODAL SOLUTION','I','X','U(I)','Q(I)')
  end
  fprintf('%3d%16.8e%16.8e%16.8e\n',i,nodes(i,1),u(i,1),Q(i,1))
end

fprintf('\nSolution: max u = %.5e\n',max(u))