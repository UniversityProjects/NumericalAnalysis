function [ ue ] = ue( x )

%ue = sin(5*pi*x);
ue = -x.*(x - 1).*(log(x + 1) + sin(5.*x) + x.^3);



end

