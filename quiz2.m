% Practice 2.2 Quiz 2
% Compute the value at the point x = 1 for the solution of the differential
% equation 1.3 u'' + 0.4 u = 0, for 0 < x < 2 with boundary conditions u(0) 
% = 1.3, u'(0) = 2. Use 4 quadratic elements to approximate the solution. 

clearvars
close all

a=0.0;
b=2.0;

a1=1.3;
a0=-0.4;
ff=0.0;

u0=1.3; %Boundary conditions
duL=2;

numDiv=4;

%solExacta:
omega=sqrt(-a0/a1);
A=u0;
B=(duL+A*omega*sin(omega*b))/(omega*cos(omega*b));
U=@(x) A*cos(omega*x)+B*sin(omega*x);

%Geometry
h=(b-a)/numDiv;
nodes=(a:0.5*h:b)';
elem = zeros(numDiv,3);
for e=1:numDiv
    k=2*e-1;
    elem(e,:)=[k,k+1,k+2];
end
numNodes=size(nodes,1);
numElem=size(elem,1);

%Assembly of the global system
K=zeros(numNodes);
F=zeros(numNodes,1);
Q=zeros(numNodes,1); 

Ke=a1/(3*h)*[7,-8,1;-8,16,-8;1,-8,7]+ ...
        a0*h/30*[4,2,-1;2,16,2;-1,2,4];

Fe=ff*h/6.0*[1;4;1];

for e=1:numElem   
    rows=[elem(e,1);elem(e,2);elem(e,3)];
    cols=rows;  
    K(rows,cols)=K(rows,cols)+Ke;
    F(rows,1)=F(rows,1)+Fe;
end

%B.C.
fixedNods=1;
freeNods=setdiff(1:numNodes,fixedNods);

u=zeros(numNodes,1);

%Natural B.C.:
Q(2:numNodes-1)=0.0; %Not necessary: Q was initialised to zero
Q(numNodes)=a1*duL;

%Essential B.C.
u(fixedNods)=u0;

%Reduced system
Qm=Q(freeNods)+F(freeNods)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);
um=Km\Qm;
u(freeNods)=um;
    
%Print the nodal solution, the value of the exact solution at 
%the nodes and the error of the nodal soluiton at each node, 
%defined as err(i) = abs(u(i)-U(nodes(i,1))), for i=1,...9
 
for i=1:numNodes
    if i == 1
        fprintf('\n%42s\n%3s%9s%18s%18s%12s\n',...
            'NODAL SOLUTION','I','x','U(I)','u(x)','ERR')
    end
    fprintf('%3d%16.8e%16.8e%16.8e%13.5e\n',...
     i,nodes(i,1),u(i,1),U(nodes(i,1)),abs(u(i,1)-U(nodes(i,1))))
end

%As x=1 is the position of node 5, the Approximate solution for U(1) is
fprintf('\nx = 1 is the position of node 5, so u(1) %s %.8e\n',...
    char(8776),u(5))
