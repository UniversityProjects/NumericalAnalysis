% risolviamo ru_t-(cu_x)_x=f
%            u(0,t)=alpha
%            u(1,t)=beta
%            u(x,0) = u_0(x)
% su una mesh non uniforme
%
% dati: N, f

clear all
close all

N = 10;

% Condizioni Al Bordo Di Dirichlet
%
alpha = 0;
beta = 0;

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


% creiamo un vettore con dentro gli h
%
% h_i = x_i - x_{i-1}
%
h = zeros(1,N);
%
for i=1:N
    if i==1
        h(1) = x(1)-0;
        % h(1) = x(1)-x(0)
    elseif i==N
        h(N) = 1-x(N-1);
        % h(N) = x(N)-x(N-1)
    else
        h(i) = x(i)-x(i-1);
    end
end

% matrice Kh
%
Kh = zeros(N-1,N-1);
%
for i=1:N-1
    % costruiamo la riga i-esima
    if i==1
        % riga 1
        xms = ( 0   + x(1))/2;
        xmd = (x(1) + x(2))/2;
        Kh(1,1) = +c(xms)/h(1) + c(xmd)/h(2);
        Kh(1,2) = -c(xmd)/h(2);
    elseif i==N-1
        % riga N-1 (ultima)
        xms = (x(N-2) + x(N-1))/2;
        xmd = (x(N-1) +   1   )/2;
        Kh(N-1,N-2) = -c(xms)/h(N-1);
        Kh(N-1,N-1) = +c(xms)/h(N-1) + c(xmd)/h(N);
    else
        % riga generica i (con 3 elementi)
        xms = (x(i-1) +  x(i) )/2;
        xmd = ( x(i)  + x(i+1))/2;
        Kh(i,i-1) = -c(xms)/h( i );
        Kh(i,i)   = +c(xms)/h( i ) + c(xmd)/h(i+1);
        Kh(i,i+1) = -c(xmd)/h(i+1);
    end
end

% termine noto fh
%
fhT = zeros(N-1,1);
fhM = zeros(N-1,1);
fhS = zeros(N-1,1);
%
% usiamo i trapezi!
%
for i=1:N-1
    fhT(i) = (h(i)+h(i+1))/2 * f(x(i));
end

% punto medio
%
for i=1:N-1
    if i==1
        xms = ( 0   + x(1))/2;
        xmd = (x(1) + x(2))/2;
        fhM(1) = ...
            h(1)/2 * f(xms)+...
            h(2)/2 * f(xmd);    
        fhS(1) = ...
            h(1)/6 * ...
            (0       + 4*(f(xms)/2) + f(x(1)))+...
            h(2)/6 * ...
            (f(x(1)) + 4*(f(xmd)/2) + 0);
    elseif i==N-1
        xms = (x(N-2) + x(N-1))/2;
        xmd = (x(N-1) +   1   )/2;
        fhM(N-1) = ...
            h(N-1)/2 * f(xms)+...
            h(N)/2   * f(xmd);
        fhS(N-1) = ...
            h(N-1)/6 * ...
            (0         + 4*(f(xms)/2) + f(x(N-1)))+...
            h(N)/6 * ...
            (f(x(N-1)) + 4*(f(xmd)/2) + 0);
    else
        xms = (x(i-1) +  x(i) )/2;
        xmd = ( x(i)  + x(i+1))/2;
        fhM(i) = ...
            h(i)/2   * f(xms)+...
            h(i+1)/2 * f(xmd);
        fhS(i) = ...
            h(i)/6 * ...
            (0       + 4*(f(xms)/2) + f(x(i)))+...
            h(i+1)/6 * ...
            (f(x(i)) + 4*(f(xmd)/2) + 0);
    end
end    
%
% modifica per condizioni non omogenee
%
xms = ( 0   + x(1))/2;
fhT(1) = fhT(1) - (-alpha/h(1))*c(xms);
fhM(1) = fhM(1) - (-alpha/h(1))*c(xms);
fhS(1) = fhS(1) - (-alpha/h(1))*c(xms);
%
xmd = (x(N-1) +   1   )/2;
fhT(N-1) = fhT(N-1) - (-beta/h(N))*c(xmd);
fhM(N-1) = fhM(N-1) - (-beta/h(N))*c(xmd);
fhS(N-1) = fhS(N-1) - (-beta/h(N))*c(xmd);
%

% risolviamo il sistema lineare
%
uhT = Kh\fhT;
uhM = Kh\fhM;
uhS = Kh\fhS;
% disegniamo la soluzione
%
lw = 2;
%
plot([0 x 1],[alpha uhT' beta],'.-b','LineWidth',lw)
hold on

plot([0 x 1],[alpha uhM' beta],'.-g','LineWidth',lw)

plot([0 x 1],[alpha uhS' beta],'.-k','LineWidth',lw)

% alternativa:
% plot([0 x 1],[0; uh; 0])

hold on
% campionamenti della soluzione esatta
%xc = linspace(0,1,10*N);
%plot(xc,ue(xc),'r','LineWidth',lw)

%legend(...
%    'f con trapezi',...
%    'f con punto medio',...
%    'f con simpson', ...
%    'soluzione esatta');


legend(...
    'f con trapezi',...
    'f con punto medio',...
    'f con simpson');

title('Soluzione Stazionaria');


% Calcolo Della Matrice Mh
% uso i trapezi in modo da avere Mh diagonale
%
Mh = zeros(N-1,N-1);
%
for i=1:N-1
    Mh(i,i) = (h(i)+h(i+1))/2 * rho(x(i));
end
%
% si ha solo la ODE
%

% scegliamo una fh e una uh

fh = fhT;
uh = uhT;

% Eulero Implicito
%
Tmax = 1;
%
dt = 0.1;
Nk = round(Tmax/0.1);
%
uhk = zeros(N-1,Nk); % Matrice che contiene tutti i vettori uhk per ogni k
%
% uhk(:,k) soluzione al passo k
% ovvero la soluzione al tempo k*dt

% uhk(:0) corrisponde al dato iniziale
% uhk(:,1) la costruiamo a parte

% calcolo uh0 = uhk(:,0)
%
uh0 = zeros(N-1,1);
%
for i=1:N-1
    uh0(i) = u0(x(i));
end

% uh0=uhk(:,0) ----> uhk(:,1)
 A = ((1/dt)*Mh)+Kh; % A in realtà è sempre costante nel nostro caso
 b = fh +(1/dt)*(Mh*uh0);
 uhk(:,1) = A\b;
    
for k=1:Nk
    % uhk(:,k) ----> uhk(:,k+1)
    %   
    A = (1/dt)*Mh+Kh; % A in realtà è sempre costante nel nostro caso
    b = fh +(1/dt)*(Mh*uhk(:,k));
    uhk(:,k+1) = A\b;
end

% disegnare le uhk
figure(2);
% disegnamo il dato iniziale
plot([0 x 1],[alpha uh0' beta],'.-b','LineWidth',lw)
hold on;
% disegnamo le uhk
for k=1:Nk
    plot([0 x 1],[alpha uhk(:,k)' beta],'.-r','LineWidth',lw)
end
% disegnamo la soluzione stazionaria
plot([0 x 1],[alpha uh' beta],'.-b','LineWidth',lw)


% disegno bidimensionale
figure(3);

for k=1:Nk
    for i=1:N-2
        X = [x(i+1)      x(i+1)       x(i)      x(i)];
        T = [k*dt       (k+1)*dt     (k+1)*dt   k*dt];
        U = [uhk(i+1,k) uhk(i+1,k+1) uhk(i,k+1) uhk(i,k)];
        patch(X,T,U,'w')
    end
end

         
            



