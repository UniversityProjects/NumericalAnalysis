% Solve -u'' = f
%
% Dirichlet Conditions
% u(0) = u(1) = 0
% with uniform mesh
% and trapezis method
%
% Input:
% N -> Number of nodes
% f -> Force field


clear all
close all

% Nodes'number definition
N = 10;

% Step definition
h = 1/N;

% Create a vector x with the internal points
x = zeros (1, N-1); % Vector Initialisation
for i=1:N-1
    x(i) = i*h;
end

% Create the Kh matrix
Kh = zeros (N-1, N-1);
for i=1:N-1
    if i==1 % First row
      Kh(1,1) = 2/h;
      Kh(1,2) = -1/h;
    elseif i==N-1 % Last row
        Kh(N-1,N-2) = -1/h;
        Kh(N-1,N-1) = 2/h;
    else % General row
        Kh(i,i-1) = -1/h;
        Kh(i,i) = 2/h;
        Kh(i,i+1) = -1/h;
    end
end

% Field force definition
fh = zeros(N-1,1);
for i=1:N-1
    fh(i) = h*f4(x(i));
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

% Plot the exact solution
xc = linspace(0,1,10*N);
plot(xc,sin(5*pi*xc),'r');

