function [faces, verts, p]=surfArrow2(crdstart, crdend, radline, radarrow, lenarrow, varargin)
[vcolor, cord0, vres]=init(0, varargin{:});
len=norm(crdend-crdstart);
xs=[linspace(0, lenarrow,10), lenarrow, lenarrow, lenarrow, linspace(lenarrow, len, 10), len, len,len]-.5*len;
ys=zeros(size(xs));
zs=zeros(size(xs));
rcirc=[linspace(0, radarrow, 10), linspace(radarrow*.75, radline, 3), radline*ones(1,10), .75*radline, .5*radline, 0]';
[x3d, y3d, z3d]=string2volumeEqu(xs, ys, zs, rcirc, vres(1), vres(2));
vctarrow=(crdend-crdstart);
vctarrow=vctarrow/norm(vctarrow);
if norm(cross(vctarrow, [1,0,0]))==0
if vctarrow(1)<0, matrot=eye(3);else, matrot=-eye(3); end;
else
xp=-vctarrow;
zp=cross(-vctarrow, [1,0,0]);
zp=zp/norm(zp);
yp=cross(zp, xp);
matrot=[xp; yp; zp]';
end

x3r=x3d;
y3r=y3d;
z3r=z3d;
cordmat=matrot*([x3d(:), y3d(:), z3d(:)]');
x3r(:)=cordmat(1,:)';
y3r(:)=cordmat(2,:)';
z3r(:)=cordmat(3,:)';

cord0=matrot*[.5*len, 0, 0]';
x3r=x3r-cord0(1)+crdstart(1); y3r=y3r-cord0(2)+crdstart(2); z3r=z3r-cord0(3)+crdstart(3);
[verts, faces] = rendpatch3( x3r, y3r, z3r, 'period', 1 );
p=patch('Faces', faces, 'Vertices', verts, 'FaceColor', vcolor);
set(p, 'EdgeColor', 'none');

if ~numel(findall(gcf,'Type','Light'))
	lightangle(10, 30);
end
lighting phong;

end


function [vcolor, cord0, vres]=init(vphi, varargin)

for i=2:2:nargin
switch varargin{i-1}
case 'cmap'
cmap=varargin{i};
indx=floor(mod(vphi, 2*pi)/(2*pi)*size(cmap, 1))+1;
vcolor=cmap(indx, :);
case 'drift'
cord0=varargin{i};
case 'Resolution'
vres=varargin{i};
end
end

if ~exist('vcolor', 'var'), vcolor=ones(numel(vphi(:)), 1)*[.65, .75, 1]; end;
if ~exist('cord0', 'var'), cord0=[0,0,0]; end;
if ~exist('vres', 'var'), vres=[128, 64]; elseif numel(vres)<2, vres=[vres, vres]; end;

end


