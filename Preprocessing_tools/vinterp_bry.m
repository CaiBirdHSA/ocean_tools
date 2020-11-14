function vinterp_bry(bryname,grdname,Zbryname,vname,obcndx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Vertical interpolation from a Z-grid to a sigma-grid in the
%  case of boundary (bry) files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% open the grid file  
% 
ng=netcdf(grdname,'r');
L=length(ng('xi_rho'));
M=length(ng('eta_rho'));
if obcndx==1
  h=ng{'h'}(1,:);
elseif obcndx==2
  h=ng{'h'}(:,L);
elseif obcndx==3
  h=ng{'h'}(M,:);
elseif obcndx==4
  h=ng{'h'}(:,1);
end
close(ng);
%
% open the boundary file  
% 
nc=netcdf(bryname,'write');
theta_s = nc{'theta_s'}(:);
theta_b =  nc{'theta_b'}(:);
hc  =  nc{'hc'}(:);
N =  length(nc('s_rho'));
vtransform = nc{'Vtransform'}(:);
if ~exist('vtransform')
    vtransform=1; %Old Vtransform
    disp([' NO VTRANSFORM parameter found'])
    disp([' USE VTRANSFORM default value vtransform = 1'])
end
%
% open the oa file  
% 
    noa=netcdf(Zbryname,'r');
z=-noa{'Z'}(:);
t=noa{'bry_time'}(:);
tlen=length(t);
Nz0=length(z);
%
% Get the sigma depths
%
zcroco=squeeze(zlevs(h,0.*h,theta_s,theta_b,hc,N,'r',vtransform));
zmin=min(min(zcroco));
zmax=max(max(zcroco));
%
% Check if the min z level is below the min sigma level 
%    (if not add a deep layer)
%
addsurf=max(z)<zmax;
addbot=min(z)>zmin;
if addsurf
 z=[100;z];
end
if addbot
 z=[z;-100000];
end
Nz=min(find(z<zmin));
z=z(1:Nz);
%
% loop on time
%
for l=1:tlen
%for l=1:1
%  disp([' Time index: ',num2str(l),' of total: ',num2str(tlen)])
  var=squeeze(noa{vname}(l,:));
%  disp(['SIZE VAR=',num2str(size(var))])
  if addsurf
    var=cat(1,var(1,:),var);
  end
  if addbot
    var=cat(1,var,var(end,:));
  end
  var=var(1:Nz,:);
  nc{vname}(l,:,:)=ztosigma_1d(flipdim(var,1),zcroco,flipud(z));
  nc{vname}(12+l,:,:)=ztosigma_1d(flipdim(var,1),zcroco,flipud(z));
  nc{vname}(24+l,:,:)=ztosigma_1d(flipdim(var,1),zcroco,flipud(z));
  nc{vname}(36+l,:,:)=ztosigma_1d(flipdim(var,1),zcroco,flipud(z));
end
close(nc);
close(noa);


%DEBUG
%
%size(var)
%var
%
%size(flipdim(var,1))
%flipdim(var,1)
%
%size(z)
%z
%
%size(flipud(z))
%flipud(z)
%
%size(zcroco)
%zcroco
%
%size(ztosigma_1d(flipdim(var,1),zcroco,flipud(z)))
%
%
%%
return
