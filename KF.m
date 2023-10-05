clearvars
close all

syms x x1 x2 h a1 a0 f0 Psi DPsi
Psi = @(x) [(x-x2)/(x1-x2), (x-x1)/(x2-x1)];
DPsi = @(x) diff(Psi(x),x);

K1 = int(a1*x*DPsi(x)'*DPsi(x),x,x1,x2)
K0 = int(a0*Psi(x)'*Psi(x),x,x1,x2)
F = int(f0*x*Psi(x)',x,x1,x2)