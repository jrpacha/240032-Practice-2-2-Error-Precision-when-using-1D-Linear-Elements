% Prob 12 (class problem's list)

clearvars
close all

a=0;
b=1.2;
numElem=3;
h=(b-a)/numElem;
nodes=(a:h:b)';
numNod=size(nodes,1);
elem=zeros(numElem,2);

A0=1.0e-2;
E=9.0e9;

for e=1:numElem
    elem(e,:)=[e,e+1];
end

K=zeros(numNod);
coef=A0*E/h;
u=zeros(4,1);
for e=1:numElem
    x1=nodes(elem(e,1),1);
    x2=nodes(elem(e,2),1);
    Ke= (1.0+0.25*(x1+x2))*[1 -1;-1 1];
    row=[elem(e,1),elem(e,2)];
    col=row;
    K(row,col) = K(row,col)+Ke; 
end
format short e
K=coef*K

F1=11.50;
F2=25.87;
F4=16.53;
P=4.65;

syms x;
f=@(t) 53.9*(1.0+t/2.0);
F3=int(f(x),0,1.2)-F1-F2-F4;
%F3=eval(F3);
clear x;
F=[F1;F2;F3;F4]

fixedNod=4;
freeNod=setdiff(1:numNod,fixedNod);

%B.C.:
Q=[0.0;P;0.0;0.0]
u(fixedNod)=0.0;

%Reduced system:
Qm=F(freeNod)+Q(freeNod)-K(freeNod,fixedNod)*u(fixedNod);
Km=K(freeNod,freeNod);

um=Km\Qm;
u(freeNod)=um;

fprintf('\n%5s%12s\n','Nod.','Displ.(m)');
fprintf('%3d%14.5e\n',[(1:numNod)',u]');

%Shape functions:
%PsiKI: is the shape function of
% element K (K=1,2,3), node I (I=1,2).
Psi11=@(x) (x-nodes(2,1))/(nodes(1,1)-nodes(2,1));
Psi12=@(x) (x-nodes(1,1))/(nodes(2,1)-nodes(1,1));
Psi21=@(x) (x-nodes(3,1))/(nodes(2,1)-nodes(3,1));
Psi22=@(x) (x-nodes(2,1))/(nodes(3,1)-nodes(2,1));
Psi31=@(x) (x-nodes(4,1))/(nodes(3,1)-nodes(4,1));
Psi32=@(x) (x-nodes(3,1))/(nodes(4,1)-nodes(3,1));

%Approximate displacement at x=0.6: we interpolate using 
%the shape functions Psi21 and Psi22, i.e, the ones associated 
%to the 2nd. elemenet, to which the point x=0.6 belongs.
x=0.6;
U=u(2)*Psi21(x)+u(3)*Psi22(x);

fprintf('\nApproximate displacement at x=%4.1f m: ',x)
fprintf('U = %11.5e m.\n',U)