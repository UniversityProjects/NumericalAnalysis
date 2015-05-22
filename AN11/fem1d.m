% Solve ru_t-(cu_x)_x=f
%
% Dirichlet Conditions
% 			 u(0,t)=alpha
%            u(1,t)=beta
% Initial Condition
%    		 u(x,0) = u_0(x)
%
% with uniform and random mesh
%
% With the following methods
% - trapezoid method
% - medium point method
% - simpson method
%
% Input:
% N -> Number of nodes
% f -> Force field

clear all
close all

N = 5;

% Dirichlet Boundary Condition
%
alpha = 0;
beta = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X Definition (Mesh)

% random mesh
%x = unique(sort(rand(1,N)));
%N = length(x)+1;

% uniform mesh
x = zeros(1,N-1);

for i=1:N-1
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

Kh = zeros(N-1,N-1); % Kh Matrix Definition

% For Loop, Kh(i) Calculation
for i=1:N-1
    if i==1 % First Row
        xms=( 0   + x(1))/2;
        xmd=(x(1) + x(2))/2;
        Kh(1,1) = +c(xms)/h(1) + c(xmd)/h(2);
        Kh(1,2) = -c(xmd)/h(2);
    elseif i==N-1 % Last Row
        xms=(x(N-2) + x(N-1))/2;
        xmd=(x(N-1) +   1   )/2;
        Kh(N-1,N-2) = -c(xms)/h(N-1);
        Kh(N-1,N-1) = +c(xms)/h(N-1) +c(xmd)/h(N);
    else % Generic Row (with 3 non-zero elements)
        xms=(x(i-1) +  x(i) )/2;
        xmd=( x(i)  + x(i+1))/2;
        Kh(i,i-1) = -c(xms)/h( i );
        Kh(i,i)   = +c(xms)/h( i ) + c(xmd)/h(i+1);
        Kh(i,i+1) = -c(xmd)/h(i+1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fh Array

fhT = zeros(N-1,1); % FhT Array Definition, Trapezoid
fhM = zeros(N-1,1); % FhM Array Definition, Medium Point
fhS = zeros(N-1,1); % FhS Array Definition, Simpson


% For Loop Trapezoid
for i=1:N-1
    fhT(i) = (h(i)+h(i+1))/2 * f(x(i));
end


% For Loop Medium Point
for i=1:N-1
    if i==1
        xms=( 0   + x(1))/2;
        xmd=(x(1) + x(2))/2;
        fhM(1) = ...
            h(1)/2 * f(xms)+...
            h(2)/2 * f(xmd);
        fhS(1) = ...
            h(1)/6 * (0 + 4*(f(xms)/2)+ f(x(1))) + ...
            h(2)/6 * (f(x(1)) + 4*(f(xmd)/2) + 0);            
    elseif i==N-1
        xms=(x(N-2) + x(N-1))/2;
        xmd=(x(N-1) +   1   )/2;
        fhM(N-1) = ...
            h(N-1)/2 * f(xms)+...
            h(N)/2   * f(xmd);
    else
        xms=(x(i-1) +  x(i) )/2;
        xmd=( x(i)  + x(i+1))/2;
        fhM(i) = ...
            h(i)/2   * f(xms)+...
            h(i+1)/2 * f(xmd);
    end
end    


% For Loop Simpson
for i=1:N-1
    if i==1
        xms=( 0   + x(1))/2;
        xmd=(x(1) + x(2))/2;
        fhS(1) = ...
            h(1)/6 * (0 + 4*(f(xms)/2)+ f(x(1))) + ...
            h(2)/6 * (f(x(1)) + 4*(f(xmd)/2) + 0);            
    elseif i==N-1
        xms=(x(N-2) + x(N-1))/2;
        xmd=(x(N-1) +   1   )/2;
        fhS(N-1) = ...
            h(N-1)/6 * (0 + 4*(f(xms)/2)+ f(x(N-1))) + ...
            h(N)/6 * (f(x(N-1)) + 4*(f(xmd)/2) + 0);
    else
        xms=(x(i-1) +  x(i) )/2;
        xmd=( x(i)  + x(i+1))/2;
        fhS(i) = ...
            h(i)/6 * (0 + 4*(f(xms)/2)+ f(x(i))) + ...
            h(i+1)/6 * (f(x(i)) + 4*(f(xmd)/2) + 0);
    end
end    


% Edit Fh For Dirichlet Non Omogenous Conditions
xms=( 0   + x(1))/2;
xmd=(x(N-1) +   1   )/2;
fhT(1) = fhT(1) - (-alpha/h(1))*c(xms);
fhM(1) = fhM(1) - (-alpha/h(1))*c(xms);
fhS(1) = fhS(1) - (-alpha/h(1))*c(xms);
fhT(N-1) = fhT(N-1) - (-beta/h(N))*c(xmd);
fhM(N-1) = fhM(N-1) - (-beta/h(N))*c(xmd);
fhS(N-1) = fhS(N-1) - (-beta/h(N))*c(xmd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear System Solution

uhT = Kh\fhT;
uhM = Kh\fhM;
uhS = Kh\fhS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution Plot

lw = 2; % LineWidth

% Trapezoid Plot
plot([0 x 1],[alpha uhT' beta],'.-b','LineWidth',lw)
hold on

% Medium Point Plot
plot([0 x 1],[alpha uhM' beta],'.-g','LineWidth',lw)
hold on

% Simpson Plot
plot([0 x 1],[alpha uhS' beta],'.-k','LineWidth',lw) 
% alternative: plot([0 x 1],[0; uh; 0])
hold on

% Legend Creation
legend(...
    'f trapezoid',...
    'f medium point',...
    'f simpson');
    
title('Stationary Solution');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mh Matrix

Mh = zeros(N-1,N-1); % Mh Matrix Definition

% For Loop, Mh(i) Calculation
% Trapezoid method to obtain a diagonal matrix
for i=1:N-1
    Mh(i,i) = (h(i)+h(i+1))/2 * rho(x(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Solution: Esplicit Euler

% Select the trapedoid f and u to work with
fh = fhT;
uh = uhT;

Tmax = 1; % Maxiumum Time
dt = 0.01; % Time Interval
Nk = round(Tmax/dt); % Step number (use round to obtain an integer)

% Matrix that holds all the uhk arrays
% where uhk(:,k) is the solution at step k (at time k*dt)
uhk = zeros(N-1,Nk); 

% uhk(:0) is the initial condition
% uhk(:,1) must be calculated out of the loop

% uh0 = uhk(:,0) definition
uh0 = zeros(N-1,1);

% uh0 = uhk(:,0) calculation
for i=1:N-1
    uh0(i) = u0(x(i));
end

% uh0 = uhk(:,0) ----> uhk(:,1)
% A is always a constant in our case, but in theory
% can depend upon the time
A = (1/dt)*Mh;
b = fh + ((1/dt)*Mh - Kh)*uh0;
uhk(:,1) = A\b;
    
    
% uhk(:,k) ----> uhk(:,k+1)
% A is always a constant in our case, but in theory
% can depend upon the time
for k=1:Nk
    A = (1/dt)*Mh; % Diagonal Matrix With Trapezoid Method
    b = fh + ((1/dt)*Mh - Kh)*uhk(:,k);
    uhk(:,k+1) = A\b;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Solution: Plot

figure(2);

% Initial data plot
plot([0 x 1],[alpha uh0' beta],'.-b','LineWidth',lw)
hold on;

% uhk plot
for k=1:Nk
    plot([0 x 1],[alpha uhk(:,k)' beta],'.-r','LineWidth',lw)
end

% stationary solution plot
plot([0 x 1],[alpha uh' beta],'.-b','LineWidth',lw)


% Bidimensional Plot
figure(3);


% General Case (i Generic)
for k=0:Nk        
    for i=1:N-2
        X = [x(i)      x(i+1)       x(i+1)       x(i)];
        T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
        if (k == 0) % Time = 0
            U = [uh0(i)    uh0(i+1)     uhk(i+1,k+1) uhk(i,k+1)];
        else % General Time Case
            U = [uhk(i,k)  uhk(i+1,k)   uhk(i+1,k+1) uhk(i,k+1)];
        end
        patch(X,T,U,'w')
    end
end

% i = 0 x_0 = 0
i = 0;
for k=0:Nk        
    X = [0         x(i+1)       x(i+1)       0];
    T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
    if (k == 0) % Time = 0
        U = [alpha    uh0(i+1)     uhk(i+1,k+1) alpha];
    else % General Time Case
        U = [alpha    uhk(i+1,k)   uhk(i+1,k+1) alpha];
    end
    patch(X,T,U,'w')
end

% i = N-1 x_N = 1
i = N-1;
for k=0:Nk        
    X = [x(i)      1             1            x(i)];
    T = [k*dt      k*dt         (k+1)*dt     (k+1)*dt];
    if (k == 0) % Time = 0
        U = [uh0(i)    beta     beta    uhk(i,k+1)];
    else % General Time Case
        U = [uhk(i,k)  beta     beta    uhk(i,k+1)];
    end
    patch(X,T,U,'w')
   
end

% Force Tridimensional View
view(3)

% Title with ODE method name
title('Esplicit Euler');
         
            



