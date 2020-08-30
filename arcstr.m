function [vptr, t] = arcstr(axialLength, radcurve, res)
Lmax=2*radcurve;
ratL=axialLength/Lmax;
anglecurve=asin(mod(axialLength, Lmax)/Lmax)+.5*pi*floor(ratL);
t=linspace(-anglecurve, anglecurve, res);
x=sin(t)*radcurve;
x=x-x(1);
y=cos(t)*radcurve;
y=y-y(1);
vptr=[x; y; zeros(1, res)]';
end
