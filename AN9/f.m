function y=f(x)


y = (sin(3.*x) + 2).*(6.*x + 25.*sin(5.*x) - 2./(x.^2 + 1) + (4.*x.^2)./(x.^2 + 1).^2) - 3.*cos(3.*x).*(5.*cos(5.*x) + (2.*x)./(x.^2 + 1) - 3.*x.^2);

%y = 1;

end
