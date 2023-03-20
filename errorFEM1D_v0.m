clearvars
close all

%Current solution
sol=@(x) -x-6*exp(1)*sinh(x)/(1-exp(1)^2); 

a=0.0;
b=1.0;
L=b-a;

a1=1.0;
a0=1.0;

nDiv = 50;

%fout=fopen('ErrorApprox.txt','w');
fout = 1;
fprintf(fout,'%10s%8s%16s\n','#Elem.','h','Error');

h = (b-a)/nDiv;
nodes = (a:h:b)';
elem = zeros(nDiv,2);
for e=1:nDiv
        elem(e,:)=[e,e+1];
end

numNodes = size(nodes,1);
numElem = size(elem,1);
K=zeros(numNodes);
F=zeros(numNodes,1);
Q=zeros(numNodes,1);
Ke=a1/h*[1,-1;-1,1]+a0*h/6.0*[2,1;1,2];

for e=1:numElem
    rows=[elem(e,1),elem(e,2)];
    cols=rows;
    x1=nodes(rows(1,1),:); %1st. node of elem e.
    x2=nodes(rows(1,2),:); %2nd. node of elem e.
    Fe=-h/6.0*[2*x1+x2;x1+2*x2];
    K(rows,cols)=K(rows,cols)+Ke;
    F(rows,1)=F(rows,1)+Fe;
end

fixedNodes=[1,numNodes];
freeNodes=setdiff(1:numNodes,fixedNodes);

%Natural B.C.:
%set Q to zero. Remark: Q(1) and Q(end), .i.e., the
%components of Q corresponding to the fixed nodes are
%computed in the post-process (see below).
Q(freeNodes)=0.0;

%Essential B.C.:
u=zeros(numNodes,1);
u(1)=0.0;
u(numNodes)=2.0;
    
%Reduced system
Qm=Q(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);
Fm=F(freeNodes)+Qm;
Km=K(freeNodes,freeNodes);
um=Km\Fm;
u(freeNodes)=um;
    
%Post process
Q=K*u-F;
    
%Error w.r.t. the exact solution.
U=sol(nodes(:,1));
error=norm(u-U,inf);%sub-inf norm (=max(abs(U-u))).
fprintf(fout,'%9d%14.6e%14.6e\n',numElem,h,error);

%fclose(fout);
%type('ErrorApprox.txt'); %print out the file to the CW.



