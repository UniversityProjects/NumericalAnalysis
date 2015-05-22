clear all
close all

syms x real

% definiamo la ue che vogliamo noi
% deve annullarsi in o e 1

ue = x*(1-x)*(sin(5*x)+log(1+x^2)-x^3);

% definisco c (> 0)
c = 2 + sin(3*x);


% calcolo f
f = -diff(c*(diff(ue,x)),x);


% vettorizza
cv = vectorize(c);
uev = vectorize(ue);
fv = vectorize(f);