clearvars
close all

a=0.0;
b=1.0;
L=b-a;

a1=1.0;
a0=1.0;
ff=3.0;

numElem=4;

%Geometry
h=L/numElem;
nodes=(a:0.5*h:b)';
numNod=size(nodes,1);    
elem=zeros(numElem,3); %Connectivity matrix
for e=1:numElem
    k=2*e-1;
    elem(e,:)=[k,k+1,k+2];
end
%Assembly of the global system
K=zeros(numNod);
F=zeros(numNod,1);
u=zeros(numNod,1);
Ke=a1/(3*h)*[7,-8,1;-8,16,-8;1,-8,7]+ ...
        a0*h/30*[4,2,-1;2,16,2;-1,2,4];
Fe=ff*h/6.0*[1;4;1];
for e=1:numElem   
    rows=[elem(e,1),elem(e,2),elem(e,3)];
    cols=rows;  
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

%Mean value (of the approximated FEM solution).
U=sum(u)/numNod;    
fprintf('Sol. Mean val.of the approx.solution,\n')
fprintf('<u> = %12.6e\n',U)