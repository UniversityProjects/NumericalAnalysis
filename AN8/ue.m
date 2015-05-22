function [ ue ] = ue( x )

%ue = sin(5*pi*x);
ue = -x.*(x - 1).*(sin(5.*x) + log(x.^2 + 1) - x.^3);



end

