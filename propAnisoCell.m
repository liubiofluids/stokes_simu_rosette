function [vV, vOmega, mtrsf, b]=propAnisoCell(omega0, mrstc, a, beta, varargin)
%% torqm is the motor torque (along the x axis); it is normalized by 4pimu
%% All forces here normalized by 4pimu
[mrsst] = init(varargin{:});

mtrsf=mrstc; % normalized by 4pimu 
mcpl=zeros(6,6);
b=zeros(6,1);
mcpl(1,1)=-mrsst(1);
mcpl(1,4)=mrsst(3);
mcpl(1,5)=mrsst(1)*a*sin(beta);
b(1)=-mrsst(3)*omega0;
mcpl(2,2)=-mrsst(2);
mcpl(2,4)=-mrsst(2)*a*sin(beta);
mcpl(2,6)=-mrsst(2)*(a*cos(beta)+mrsst(7)*.5);
mcpl(3,3)=-mrsst(2);
mcpl(3,5)=mrsst(2)*(a*cos(beta)+mrsst(7)*.5);
mcpl(4,1)=-mrsst(3);
mcpl(4,2)=-mrsst(2)*a*sin(beta);
mcpl(4,5)=-mrsst(3)*a*sin(beta);
mcpl(4,4)=-mrsst(5)-mrsst(2)*(a*sin(beta))^2;
mcpl(4,6)=mrsst(2)*a*sin(beta)*(a*cos(beta)+0.5*mrsst(7));
b(4)=mrsst(5)*omega0;
mcpl(5,5)=-mrsst(6);
mcpl(5,1)=mrsst(1)*a*sin(beta);
mcpl(5,3)=mrsst(2)*(a*cos(beta)+.5*mrsst(7));
mcpl(5,4)=-mrsst(3)*a*sin(beta);
mcpl(5,5)=-mrsst(1)*(a*sin(beta))^2-mrsst(2)*(a*cos(beta)+.5*mrsst(7))^2;
b(5)=mrsst(3)*omega0*a*sin(beta);
mcpl(6,6)=-mrsst(6);
mcpl(6,2)=-mrsst(2)*(a*cos(beta)+.5*mrsst(7));
mcpl(6,4)=-mrsst(2)*(a*cos(beta)+.5*mrsst(7))*a*sin(beta);
mcpl(6,6)=-mrsst(2)*(a*cos(beta)+.5*mrsst(7))^2;
mtrsf=mtrsf+mcpl;

vx=mtrsf\b;

vV=vx(1:3);
vOmega=vx(4:6);

end


function [mrsst] = init(varargin)
for i=2:2:nargin
    switch varargin{i-1}
        case 'Resistance'
            mrsst=varargin{i};
     end
end
if ~exist('mrsst', 'var'), mrsst=[0.59, 0.84, 0.024, 0.0, 0.0224, 3.0, 6.0]; end; % default resistance matrix of the flagellum (C. crescentus)

end


