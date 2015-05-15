% Plot errors of
% - trapezoid method
% - medium point method
% - simpson method
%


clear all
close all

Nv = [10 20 30 40 50 100 200 300];
i = 0;
hmax=zeros(1,length(Nv));
errT=zeros(1,length(Nv));
errM=zeros(1,length(Nv));
errS=zeros(1,length(Nv));

for N=Nv
    i = i + 1;
    [hmax(i), errT(i), errM(i), errS(i)] = fem1d(N);
end

close all

figure(1)
loglog(...
hmax,errT,'b*-',...
hmax,errM,'g*-',...
hmax,errS,'k*-');
   
