if ~exist('datfolder', 'var'), datfolder='./'; end;
if ~exist('folder', 'var'), folder='./'; end;

load([datfolder, 'resstwallvh.mat']); %BEM result of the resistance matrix

wallmatrix=struct('mrm', mrm, 'vh', vh);

% lubrication theory result
fNxxFU = @(h) 1.5*(-8/15*log(h)+0.96);
fNzzFU = @(h) 1.5*(h).^(-1); 
fNxxLO = @(h) 2*(-.4*log(h)+0.38);
fNzzLO = @(h) 2*ones(size(h));
fVxyLU = @(h) 2*(-.1*log(h)-.19);
fVxyFO = @(h) 1.5*(-2/15*log(h)-.25);


coefrsst=[1.1, .62, 0.0, 0.84, 0.0224, 3.0, 6.0]; %% resistance matrix elements for the flagellum of C. crescentus 

rectROI=[0   60   1160   600];
vorient0=[1,0,0];
vcell0=[-1,0,0];
pangle=pi/3; 
%vcell0=[-cos(pangle), 0, sin(pangle)];
vorient=vorient0;
vcell=vcell0;
%rot=eye(3);
rot=EulerRotation([0,-pangle,0]);
dt=1;
a=3.;
z0=2;
icutoff=1; % cutting-off index to make the plotted range of x-y within 5 times the rosette size. 
[x,y,z] = sphere;
figure
[f, v, c] =surf2patch(x,y,z,z);

arrowStlet=(rot'*[(a+.2)*vorient0; (a+.2)*vorient0+1.5*norm(vorient0)*[cos(pangle),0, sin(pangle)]]')';
sizearrow=[.1, 0.3, 0.5];
radarrow=0.5;

resbrd=[2048, 2048];
chessboard=zeros(resbrd);
nblock=[8, 8];
wblock=resbrd./nblock;

for i=1:nblock(1);
for j=1:nblock(2);
if mod(i+j, 2)==0
chessboard(floor((i-1)*wblock(1)+1):floor(i*wblock(1)), floor((j-1)*wblock(2)+1):floor(j*wblock(2)))=1;
end
end
end

PSF=fspecial('gaussian', 25, 2);
chessboard=imfilter(chessboard, PSF);
chessboard=cat(3, chessboard, chessboard, chessboard);
bdbar=5;
bdscale=nblock*bdbar;



vb1=cross([0,1,0], vcell0);
if norm(vb1)>0
vb1=vb1/norm(vb1);
vb2=cross(vcell0, vb1);
rotsp=[vb1; vb2; vcell0]';
else
rotsp=eye(3);
end

v=(rotsp*v')'*a;
v=(rot'*v')';

vx=[0,0,z0];
nstep=100000;
zmin=1E-4;
torq0=1;
v1=zeros(nstep, 3);
v1(1,:)=vx;
nctr=20;

figure(gcf);
set(gcf, 'Unit', 'centimeter');
set(gcf, 'Position', [0, 5, 42, 23.5]);
gca1=subplot(1,2,1); hold on; 
gca2=subplot(1,2,2); hold on;
set(gca1, 'Position', [0, .2, .5, .5]);
set(gca2, 'Position', [.5, .2, .5, .5]);
%lightangle(30, 20);
axis equal;
%axis([-200, 200, -200, 200, -5, 5]);
[vptrc, t] = arcstr(2*a, 1000, nctr);

vptrc(:,1)=vptrc(:,1)-a;

%[verts, faces] = string2spheroid(vptrc, a, 'Muliplier', 5, 'Order', 5, 'Resolution', 24);


if ~exist('flagmov', 'var'), flagmov=0; end;

if flagmov
aviobj=VideoWriter([folder, sprintf('tracerosea_%d_%02dth_%d_%02dL_%d.avi', floor(a), floor((a-floor(a))*100+.5), floor(pangle), floor((pangle-floor(pangle))*100+.5), floor(coefrsst(7)))]);
%aviobj.FrameRate=12;
open(aviobj)
end

for i=2:nstep

	NxxFU=fNxxFU(v1(i-1, 3)/a)*a; NyyFU=NxxFU;
	NzzFU=fNzzFU(v1(i-1, 3)/a)*a;
	VxyLU=fVxyLU(v1(i-1, 3)/a)*a^2; VyxLU=VxyLU;
	VxyFO=fVxyFO(v1(i-1, 3)/a)*a^2; VyxFO=VxyFO;
	NxxLO=fNxxLO(v1(i-1, 3)/a)*a^3; NyyLO=NxxLO;
	NzzLO=fNzzLO(v1(i-1, 3)/a)*a^3;

if ~isempty(wallmatrix)
	if v1(i-1,3)/a>=min(wallmatrix.vh) & v1(i-1,3)/a<=max(wallmatrix.vh)
		mdat=2*abs(interp2(1:8, wallmatrix.vh, wallmatrix.mrm, 1:8, v1(i-1,3)/a));
		mdat([1,5])=mdat([1,5])*a;
		mdat([2,3,6,7])=mdat([2,3,6,7])*a^2;
		mdat([4,8])=mdat([4,8])*a^3;
		NxxFU=max([NxxFU, mdat(1)]);NyyFU=NxxFU;
		NzzFU=max([NzzFU, mdat(5)]);
		VxyLU=max([VxyLU, mdat(2)]); VyxLU=VxyLU;
		VxyFO=max([VxyFO, mdat(3)]); VyxFO=VxyFO;
		NxxLO=max([NxxLO, mdat(4)]); NyyLO=NxxLO;
		NzzLO=max([NzzLO, mdat(8)]);
	end

end

mrstc=zeros(6,6);
mrstc(1,[1,5])=[-NxxFU, VxyFO];
mrstc(2,[2,4])=[-NyyFU, -VyxFO];
mrstc(3,3)=-NzzFU;
mrstc(4,[2,4])=[-VxyLU, -NxxLO];
mrstc(5,[1,5])=[VyxLU, -NyyLO];
mrstc(6,6)=-NzzLO;

mrstc(1:3,1:3)=rot'*(mrstc(1:3,1:3)*rot);
mrstc(1:3,4:6)=rot'*(mrstc(1:3,4:6)*rot);
mrstc(4:6,1:3)=rot'*(mrstc(4:6,1:3)*rot);
mrstc(4:6,4:6)=rot'*(mrstc(4:6,4:6)*rot);

[vV, vOmega]=propAnisoCell(-pi, mrstc, a, pangle); 

vV=rot*vV;
vOmega=rot*vOmega;
rot=EulerRotation(vOmega'*dt)*rot;

v1(i,:)=vV(:)'*dt+v1(i-1,:);
v1(i,3)=max([v1(i,3), 2*zmin-v1(i, 3)]);%max([zmin, v1(i,3)]);

vertsrot=(rot*v')'+repmat(v1(i,:), size(v, 1), 1);
if mod(i, 50)==0
gca1=subplot(1,2,1); set(gca1, 'Position', [0, 0, .5, 1]); axis equal ;hold on;

if exist('p', 'var'), delete(p); end;
p=patch('faces', f, 'vertices', vertsrot); 
set(p,'FaceColor', 'Interp', 'CData', c,...
'EdgeColor', 'black');
if exist('pfrc', 'var'), delete(pfrc); end;
crdfrc=(rot*arrowStlet')'+repmat(v1(i,:), 2, 1);

delete(findall(gcf, 'Type', 'Light'));
if v1(i,3)+2*a>a+5
axis([v1(i,1)-10, v1(i,1)+10, v1(i,2)-10, v1(i,2)+10, -a, v1(i,3)+2*a]);
else
axis([v1(i,1)-10, v1(i,1)+10, v1(i,2)-10, v1(i,2)+10, -a, a+5]);
end
view(20, 30);
if exist('p0', 'var'), delete(p0); end;
[xg, yg]=meshgrid(linspace(-.5, .5, resbrd(1))*bdscale(1)+floor(v1(i,1)/(2*bdbar))*(2*bdbar)-v1(i,1), linspace(-.5, .5, resbrd(2))*bdscale(2)+floor(v1(i,2)/(2*bdbar))*(2*bdbar)-v1(i,2));
chessboardr=ones(resbrd(1), resbrd(2));
chessboardr(find(sqrt(xg.^2+yg.^2)<.7*a))=0;
PSF2=fspecial('gaussian', 100, 20);
chessboardr=imfilter(chessboardr, PSF2).*chessboard(:,:,1);
chessboardt=cat(3, chessboardr, chessboardr, .5*chessboard(:,:,1)+.5*chessboardr);
chessboardt=chessboardt*.3+chessboard*.7;
p0=surf(repmat([-.5, .5]*bdscale(1), 2,1)+floor(v1(i,1)/(2*bdbar))*(2*bdbar), repmat([-.5, .5]'*bdscale(2), 1, 2)+floor(v1(i,2)/(2*bdbar))*(2*bdbar), repmat(-a, 2, 2), chessboardt , 'FaceColor', 'texturemap');
camproj('Perspective');
set(gca1, 'color', 'none')
set(gca1, 'Xcolor', 'none')
set(gca1, 'Ycolor', 'none')
set(gca1, 'Zcolor', 'none')
set(gca1,'xtick',[])
set(gca1,'xticklabel',[])
set(gca1,'ytick',[])
set(gca1,'yticklabel',[])
set(gca1,'ztick',[])
set(gca1,'zticklabel',[])

gca2=subplot(1,2,2); set(gca2, 'Position', [.5, .05, .45, .9]); axis equal ;hold on;
plot3(v1(1:i,1), v1(1:i,2), v1(1:i,3));
if exist('pf2', 'var'), delete(pf2); end;
pf2=patch('faces', f, 'vertices', vertsrot); 
set(pf2,'FaceColor', 'Interp', 'CData', c,...
'EdgeColor', 'black');
if exist('pfrc2', 'var'), delete(pfrc2); end;
crdfrc=(rot*arrowStlet')'+repmat(v1(i,:)-vcell*a, 2, 1);

while max(v1(icutoff:i,1))-min(v1(icutoff:i,1))>5*a || max(v1(icutoff:i,2))-min(v1(icutoff:i,2))>5*a
	icutoff=icutoff+1;
end
axis(gca2, [min(v1(icutoff:i,1))-10, max(v1(icutoff:i,1))+10, min(v1(icutoff:i,2))-10, max(v1(icutoff:i,2))+10, min(v1(icutoff:i,3))-2*a, max(v1(icutoff:i,3))+2*a]);
view(0, 90);

box(gca2, 'on');
camproj('Perspective');
delete(findall(gcf, 'Type', 'annotation'));
htxt=annotation('textbox', [0.5, 0.8, 0.1, .1], 'EdgeColor', 'none', 'FontSize', 20, 'String', ['$$d = $$', sprintf('%.2f', v1(i,3)), '$$\mu$$m'],'interpreter','latex');

sptR=axes('parent',gcf,'position',[.7, .7, 0.25, 0.25], 'XAxisLocation','top');
set(sptR, 'xtick', [], 'ytick', [], 'xcolor', 'none', 'ycolor', 'none', 'Color', 'none'); hold on;
figure(gcf); 
if (flagmov)
h=getframe(gcf);
    frame0=h.cdata; %imcrop(h.cdata, rectROI);
writeVideo(aviobj, frame0);
end


pause(0.1);
delete(sptR);
end
end
if exist('aviobj', 'var')
	close(aviobj)
end




