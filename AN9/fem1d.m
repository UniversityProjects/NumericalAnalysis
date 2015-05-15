function [hmax, errmaxT, errmaxM, errmaxS] = fem1d(N)


% risolviamo -(cu')'=f
%            u(0)=alpha
%            u(1)=betha
% su una mesh non uniforme
%
% dati: N, f


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
        xms=( 0   + x(1))/2;
        xmd=(x(1) + x(2))/2;
        Kh(1,1) = +c(xms)/h(1) + c(xmd)/h(2);
        Kh(1,2) = -c(xmd)/h(2);
    elseif i==N-1
        % riga N-1 (ultima)
        xms=(x(N-2) + x(N-1))/2;
        xmd=(x(N-1) +   1   )/2;
        Kh(N-1,N-2) = -c(xms)/h(N-1);
        Kh(N-1,N-1) = +c(xms)/h(N-1) +c(xmd)/h(N);
    else
        % riga generica i (con 3 elementi)
        xms=(x(i-1) +  x(i) )/2;
        xmd=( x(i)  + x(i+1))/2;
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


% punto medio e simpson
for i=1:N-1
    if i==1
        xms=( 0   + x(1))/2;
        xmd=(x(1) + x(2))/2;
        % Medium Point
        fhM(1) = ...
            h(1)/2 * f(xms)+...
            h(2)/2 * f(xmd);
        % Simpson
        fhS(1) = ...
            h(1)/6 * (0 + 4*(f(xms)/2)+ f(x(1))) + ...
            h(2)/6 * (f(x(1)) + 4*(f(xmd)/2) + 0);            
    elseif i==N-1
        xms=(x(N-2) + x(N-1))/2;
        xmd=(x(N-1) +   1   )/2;
        % Medium Point
        fhM(N-1) = ...
            h(N-1)/2 * f(xms)+...
            h(N)/2   * f(xmd);
        % Simpson
        fhS(N-1) = ...
            h(N-1)/6 * (0 + 4*(f(xms)/2)+ f(x(N-1))) + ...
            h(N)/6 * (f(x(N-1)) + 4*(f(xmd)/2) + 0);
    else
        xms=(x(i-1) +  x(i) )/2;
        xmd=( x(i)  + x(i+1))/2;
        % Medium Point
        fhM(i) = ...
            h(i)/2   * f(xms)+...
            h(i+1)/2 * f(xmd);
        % Simpson
        fhS(i) = ...
            h(i)/6 * (0 + 4*(f(xms)/2)+ f(x(i))) + ...
            h(i+1)/6 * (f(x(i)) + 4*(f(xmd)/2) + 0);
    end
end    


% modifica condizioni non omogenee
xms=( 0   + x(1))/2;
xmd=(x(N-1) +   1   )/2;
fhT(1) = fhT(1) - (-alpha/h(1))*c(xms);
fhM(1) = fhM(1) - (-alpha/h(1))*c(xms);
fhS(1) = fhS(1) - (-alpha/h(1))*c(xms);
fhT(N-1) = fhT(N-1) - (-beta/h(N))*c(xmd);
fhM(N-1) = fhM(N-1) - (-beta/h(N))*c(xmd);
fhS(N-1) = fhS(N-1) - (-beta/h(N))*c(xmd);


% risolviamo il sistema lineare
%
uhT = Kh\fhT;
uhM = Kh\fhM;
uhS = Kh\fhS;

% disegniamo la soluzione
%
plot([0 x 1],[alpha uhT' beta],'.-b')
hold on

plot([0 x 1],[alpha uhM' beta],'.-g')
hold on

plot([0 x 1],[alpha uhS' beta],'.-k')

% alternativa:
% plot([0 x 1],[0; uh; 0])

hold on
% campionamenti della soluzione esatta
xc = linspace(0,1,10*N);
plot(xc,ue(xc),'r')

legend(...
    'f con trapezi',...
    'f con punto medio',...
    'f con simpson', ...
    'soluzione esatta');


% calcolo la soluzione uh nei nodi e nei punti medi
errmaxT = 0;
errmaxM = 0;
errmaxS = 0;
for i=1:N-1
    % errore in x(i)
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
    % errore in x(i)-h(i)/2
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


% punto medio ultimo intervallo
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
