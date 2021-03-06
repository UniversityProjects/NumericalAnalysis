% Solve -u'' = f
% Dirichlet Conditions
% u(0) = u(1) = 0
% with uniform and random mesh
% and trapezoid method
%
% Input:
% N -> Number of nodes
% f -> Force field


clear all
close all

% Nodes'number definition
N = 9;

% random mesh
% x = unique(sort(rand(1,N)));
% N = length(x)+1;

% uniform mesh
x = zeros(1,N-1);

for i=1:N-1
    x(i) = (i/N)^2;
end



% Step definition
% h(i) = x(i) - x(i-1)
h = zeros(1,N);
for i=1:N
    if i==1 %h(1) = x(1)-x(0)
        h(1) = x(1)-0;
    elseif i==N %h(N) = 1-x(N-1)
        h(N) = 1-x(N-1);
    else
        h(i) = x(i)-x(i-1);
    end
end



% Create the Kh matrix
Kh = zeros (N-1, N-1);
for i=1:N-1
    if i==1 % First row
      Kh(1,1) = 1/h(1) + 1/h(2);
      Kh(1,2) = -1/h(2);
    elseif i==N-1 % Last row
        Kh(N-1,N-2) = -1/h(N-1);
        Kh(N-1,N-1) = 1/h(N-1) + 1/h(N);
    else % General row
        Kh(i,i-1) = -1/h(i);
        Kh(i,i) = 1/h(i) + 1/h(i+1);
        Kh(i,i+1) = -1/h(i+1);
    end
end

% Field force definition (Trapezoid)
fh = zeros(N-1,1);
for i=1:N-1
    fh(i) = (h(i)+h(i+1)/2) * f4(x(i));
end


% Solve the linear system
uh = Kh\fh;


% Plot the solution
% Adding the Dirichlet boundary conditions 
% x(0) = x(1) = 0;
% uh' is the transposed of the colum vector uh
% alternatively they could be added as coloums:
% [0; uh; 0]
plot([0 x 1], [0 uh' 0],'.-');

hold on;

xc = linspace(0,1,10*N);
plot(xc,sin(5*pi*xc),'r');

