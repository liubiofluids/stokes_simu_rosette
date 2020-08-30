function [matrot, mattr] = EulerRotation(vomg)
phi=norm(vomg);
if phi==0
    matrot=eye(3);
    mattr=eye(3);
    return;
end
zn=vomg/phi;
xn=cross([0,1,0], zn);
if norm(xn)==0
matrot=Rotation_Matrix(-sum(vomg.*[0,1,0]), 2);
mattr=[1,0,0; 0,0,1; 0,1,0];
return;
end
xn=xn/norm(xn);
yn=cross(zn, xn);
mattr=[xn; yn; zn]';
matrot=mattr*Rotation_Matrix(-phi, 3)*mattr';

end
