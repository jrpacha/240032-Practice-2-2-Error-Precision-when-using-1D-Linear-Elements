clearvars
close all

%Current solution
f=@(x) -x-6*exp(1)*sinh(x)/(1-exp(1)^2); 

a=0.0;
b=1.0;
L=b-a;

a1=1.0;
a0=1.0;

numElem=10;

%Geometry
h=(b-a)/numElem;
nodes=a:h:b;
nodes=nodes';
numNod=size(nodes,1);    
elem=zeros(numElem,2); %Connectivity matrix
for e=1:numElem
    elem(e,:)=[e,e+1];
end
%Assembly of the global system
K=zeros(numNod);
F=zeros(numNod,1);
u=zeros(numNod,1);
for e=1:numElem
    x1=nodes(elem(e,1)); %1st. node of elem e
    x2=nodes(elem(e,2)); %2nd. node of elem e
    Ke=0.5*a1*(x1+x2)/h*[1,-1;-1,1]+ ...
        a0*h/6.0*[2,1;1,2];
    rows=[elem(e,1),elem(e,2)];
    cols=rows;
    Fe=-h/6.0*[2*x1+x2;x1+2*x2];
    K(rows,cols)=K(rows,cols)+Ke;
    F(rows,1)=F(rows,1)+Fe;
end
fixedNod=[1,numNod];
freeNod=setdiff(1:numNod,fixedNod);
%Natural B.C.
Q=zeros(numNod,1);
%Essential B.C.
u(1)=0.0;
u(numNod)=2.0;
%Reduced system
Qm=Q(freeNod)-K(freeNod,fixedNod)*u(fixedNod);
Fm=F(freeNod)+Qm;
Km=K(freeNod,freeNod);
um=Km\Fm;
u(freeNod)=um;
    
%Post process
Q=K*u-F;

%Mean value
U=sum(u)/numNod;    
fprintf('Sol. <u> = e%14.6e\n',U)
