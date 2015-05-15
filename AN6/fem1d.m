function [hmax, errT, errM, errS] = fem1d(N)


% risolviamo -u''=f
%            u(0)=u(1)=0
% su una mesh non uniforme
%
% dati: N, f


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEGNAZIONE DELLE X

% Random Mesh
x = unique(sort(rand(1,N)));
N = length(x)+1;

% creiamo un vettore x con i punti interni
%x = zeros(1,N-1);

%for i=1:N-1
%    x(i) = (i/N);
%end

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
        Kh(1,1) = +1/h(1) + 1/h(2);
        Kh(1,2) = -1/h(2);
    elseif i==N-1
        % riga N-1 (ultima)
        Kh(N-1,N-2) = -1/h(N-1);
        Kh(N-1,N-1) = +1/h(N-1) +1/h(N);
    else
        % riga generica i (con 3 elementi)
        Kh(i,i-1) = -1/h( i );
        Kh(i,i)   = +1/h( i ) + 1/h(i+1);
        Kh(i,i+1) = -1/h(i+1);
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

% risolviamo il sistema lineare
%
uhT = Kh\fhT;
uhM = Kh\fhM;
uhS = Kh\fhS;

% disegniamo la soluzione
%
plot([0 x 1],[0 uhT' 0],'.-b')
hold on

plot([0 x 1],[0 uhM' 0],'.-g')
hold on

plot([0 x 1],[0 uhS' 0],'.-k')

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


% Errori
errT = max(abs(uhT' - ue(x) ));
errM = max(abs(uhM' - ue(x) ));
errS = max(abs(uhS' - ue(x) ));

hmax = max(h);

end



