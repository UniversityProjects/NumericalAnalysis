% Solve -(cu')' = f
%
% Dirichlet Neumann Conditions
% u(0) = 0
% u'(1) = 0
% with uniform and random mesh
% With the following methods
% - trapezoid method
% 
%
% Input:
% N -> Number of nodes
% f -> Force field

clear all;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X Definition

N = 100;

% Random Mesh
%x = unique(sort(rand(1,N)));
%N = length(x)+1;

% uniform mesh
x = zeros(1,N);
for i=1:N
    x(i) = (i/N);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h Step Array
% h(i) = x(i) - x(i-1)

h = zeros(1,N); % h Array Definition

% For Loop, h Calculation
for i=1:N
    if i==1 % h(1) = x(1)-x(0), x(0) = 0
        h(1) = x(1)-0;        
    elseif i==N % h(N) = x(N)-x(N-1), x(N) = 1
        h(N) = 1-x(N-1);        
    else
        h(i) = x(i)-x(i-1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kh Matrix

Kh = zeros(N,N); % Kh Matrix Definition

% For Loop, Kh(i) Calculation
for i=1:N-1
    if i==1 % First Row
        Kh(1,1) = +1/h(1) + 1/h(2);
        Kh(1,2) = -1/h(2);
    elseif i==N-1 % Last Row
        Kh(N,N-1) = -1/h(N-1);
        Kh(N,N) = +1/h(N-1) +1/h(N);
    else % Generic Row (with 3 non-zero elements)
        Kh(i,i-1) = -1(xms)/h( i );
        Kh(i,i)   = +1(xms)/h( i ) + 1(xmd)/h(i+1);
        Kh(i,i+1) = -1(xmd)/h(i+1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fh Array

fhT = zeros(N,1); % FhT Array Definition, Trapezoid

% For Loop Trapezoid
for i=1:N-1
    fhT(i) = (h(i)+h(i+1))/2 * f(x(i));
end

fhT(N) = h(N)/2 * f(x(N));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear System Solution

uhT = Kh\fhT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution Plot

% Trapezoid Plot
plot([0 x],[0 uhT'],'.-b')
hold on

% Exact Solution Plot
xc = linspace(0,1,10*N); % Exact Solution Sampling
plot(xc,ue(xc),'r')

 Legend Creation
legend(...
    'f trapezoid',...
    'exact solution');


% Errors
errT = max(abs(uhT' - ue(x) ));

hmax = max(h);


