clear all
close all

Nv = [10 20 30 40 50 100 200 300];
i = 0;
hmax=zeros(1,length(Nv));
errmaxT=zeros(1,length(Nv));
errmaxM=zeros(1,length(Nv));
errmaxS=zeros(1,length(Nv));

for N=Nv
    i = i + 1;
    [hmax(i), errmaxT(i), errmaxM(i), errmaxS(i)] = fem1d(N);
end

close all

figure(1)
loglog(...
hmax,errmaxT,'b*-',...
hmax,errmaxM,'g*-',...
hmax,errmaxS,'k*-');
   