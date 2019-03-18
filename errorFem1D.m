clearvars
close all

%Current solution
f=@(x) -x-6*exp(1)*sinh(x)/(1-exp(1)^2); 

a=0.0;
b=1.0;
L=b-a;

a1=1.0;
a0=1.0;

nDiv=[5,50,500,5000];

fout=fopen('ErrorApprox.txt','w');
fprintf(1,'%10s%8s%16s\n','#Elem.','h','Error')

for numElem=nDiv
    %Geometry
    h=L/numElem;
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
    Ke=a1/h*[1,-1;-1,1]+a0*h/6.0*[2,1;1,2];
        %Note: in this case in point, Ke is the same
        %      for all the elements.
    for e=1:numElem
        rows=[elem(e,1),elem(e,2)];
        cols=rows;
        x1=nodes(rows(1,1),1); %1st. node of elem e
        x2=nodes(rows(1,2),1); %2nd. node of elem e
        Fe=-h/6.0*[2*x1+x2;x1+2*x2];
        K(rows,cols)=K(rows,cols)+Ke;
        F(rows,1)=F(rows,1)+Fe;
    end
    
    fixedNod=[1,numNod];
    freeNod=setdiff(1:numNod,fixedNod);
    
    %Natural B.C.:
    %None, set Q to zero (Q(1) and Q(end) are computed
    %in the post-process (see below).
    Q=zeros(numNod,1);
    
    %Essential B.C.:
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
    
    %Error w.r.t. the exact solution.
    U=f(nodes(:,1));
    error=norm(u-U,inf);%sub-inf norm (=max(abs(U-u)))
    fprintf(1,'%9i%14.6e%14.6e\n',numElem,h,error)
end
fclose(fout);
type('ErrorApprox.txt'); %print out the file to the CW
        