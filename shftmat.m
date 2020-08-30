function nmat=shftmat(omat, row, column)
nmat=zeros(size(omat));
ind=0:1:numel(omat)-1;
ncl=numel(omat(1,:));
nrw=numel(omat(:,1));

nind= mod(floor(ind/nrw)-column, ncl)*nrw + mod (mod(ind, nrw)-row, nrw)+1;
nmat(nind)=omat(:);