function [hmax, errmaxT, errmaxM, errmaxS] = fem1d(N)

% Solve -(cu')' = f
%
% Dirichlet Conditions
% 			 u(0)=alpha
%            u(1)=beta

% with uniform and random mesh

% With the following methods
% - trapezoid method
% - medium point method
% - simpson method
%
% Input:
% N -> Number of nodes
% f -> Force field


alpha = ue(0);
beta = ue(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEGNAZIONE DELLE X


% Random Mesh
%x = unique(sort(rand(1,N)));
%N = length(x)+1;

% creiamo un vettore x con i punti interni
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


% Edit Fh For Non Omogenous Conditions
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

% Trapezoid Plot
plot([0 x 1],[alpha uhT' beta],'.-b')
hold on

% Medium Point Plot
plot([0 x 1],[alpha uhM' beta],'.-g')
hold on

% Simpson Plot
plot([0 x 1],[alpha uhS' beta],'.-k') % alternative: plot([0 x 1],[0; uh; 0])
hold on

% Exact Solution Plot
xc = linspace(0,1,10*N); % Exact Solution Sampling
plot(xc,ue(xc),'r')


% Legend Creation
legend(...
    'f trapezoid',...
    'f medium point',...
    'f simpson', ...
    'exact solution');


% Solution uh calculated in nodes and in medium points
errmaxT = 0;
errmaxM = 0;
errmaxS = 0;
for i=1:N-1
    % error in x(i)
    errT = abs(uhT(i) - ue(x(i)));
    errM = abs(uhM(i) - ue(x(i)));
    errS = abs(uhS(i) - ue(x(i)));
    if (errT > errmaxT)
        errmaxT = errT;
    end
    if (errM > errmaxM)
        errmaxM = errM;
    end
    if (errS > errmaxS)
        errmaxS = errS;
    end
    % error in x(i)-h(i)/2
    if i==1
        errT = abs((uhT(i)+alpha)/2 - ue(x(i)-h(i)/2));
        errM = abs((uhM(i)+alpha)/2 - ue(x(i)-h(i)/2));
        errS = abs((uhS(i)+alpha)/2 - ue(x(i)-h(i)/2));
    else
         errT = abs((uhT(i)+uhT(i-1))/2 - ue(x(i)-h(i)/2));
         errM = abs((uhM(i)+uhM(i-1))/2 - ue(x(i)-h(i)/2));
         errS = abs((uhS(i)+uhS(i-1))/2 - ue(x(i)-h(i)/2));
    end
    if (errT > errmaxT)
        errmaxT = errT;
    end
    if (errM > errmaxM)
        errmaxM = errM;
    end
    if (errS > errmaxS)
        errmaxS = errS;
    end
end


% medium point last intervall
errT = abs(uhT(N-1)+beta)/2-ue(1-h(N)/2);
errM = abs(uhM(N-1)+beta)/2-ue(1-h(N)/2);
errS = abs(uhS(N-1)+beta)/2-ue(1-h(N)/2);
if (errT > errmaxT)
        errmaxT = errT;
 end
 if (errM > errmaxM)
    errmaxM = errM;
 end
 if (errS > errmaxS)
    errmaxS = errS;
 end
    

hmax = max(h);

end
