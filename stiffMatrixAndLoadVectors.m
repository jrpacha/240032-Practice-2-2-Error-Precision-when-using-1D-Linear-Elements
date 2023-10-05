clearvars
close all

% Linear elements and a1(x)= a1, a0(x) = a0, f(x) = f*x 
h = (b-a)/div;
Ke = a1*[1,-1; -1, 1]/h + a0*h*[2,1;1,2]/6;            %outside for loop
x1 = nodes(rows(1),1); x2 = nodes(rows(2),1);          %inside for loop
Fe = f*h*[2*x1 + x2; x1+2*x2]/6.0;                     %inside for loop

% Linear elements and a1(x) = a1*x, a0(x) = a0, f(x) = f*x;
x1 = nodes(rows(1),1); x2 = nodes(rows(2),1);          %inside for loop
Ke = a1*(x1 + x2)*[1,-1;-1,1]/h/2.0 + ...
    a0*h*[2,1;1,2]/6;                                  %inside for loop
Fe = -h*[2*x1 + x2; x1+2*x2]/6.0;                      %inside for loop

%Quadratic elements and a1(x) = a1, a0(x) = a0, f(x) = f;
Ke = a1*[7, -8, 1; -8, 16, -8; 1, -8, 7]/h/3.0 + ...   
    a0*h*[4, 2, -1; 2, 16, 2; -1, 2, 4]/30.0;          %outside for loop
Fe = f*h*[1; 4; 1]/6.0;                                %outside for loop
