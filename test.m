a = 1.0;
sigma = 2.0;
r = 0.2:0.001:3;
phi = 0:pi/100:2*pi;
[R,PHI] = meshgrid(r,phi);
X = R.*cos(PHI);
Y = R.*sin(PHI);
SIGMAx = sigma*sqrt(a./(2*R)).*cos(PHI/2).*(1-sin(PHI/2).*sin(3/2*PHI));
contourf(X,Y,SIGMAx)