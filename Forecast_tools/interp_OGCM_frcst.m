function interp_OGCM_frcst(OGCM_prefix,OGCM_name,Roa,interp_method,...
                           lonU,latU,lonV,latV,lonT,latT,Z,tin,...
		           nc_clm,nc_bry,lon,lat,angle,h,tout,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read the local OGCM files and perform the interpolations
%
% Ok, I am lazy and I did not do something special for the bry files...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conserv=1; % same barotropic velocities as the OGCM
%
disp(['  Horizontal interpolation: ',OGCM_name.name])
%
%
% CROCO grid angle
%
cosa=cos(angle);
sina=sin(angle);
%
% Open the OGCM file (1:el,2:s,3:ua,4:va,5:t,6:u,7:v)
%
nc1=netcdf([OGCM_prefix,OGCM_name(1).name]);    %nc1:el
nc2=netcdf([OGCM_prefix,OGCM_name(3).name]);    %nc2:ua
nc3=netcdf([OGCM_prefix,OGCM_name(4).name]);    %nc3:va
nc4=netcdf([OGCM_prefix,OGCM_name(6).name]);    %nc4:u3d
nc5=netcdf([OGCM_prefix,OGCM_name(7).name]);    %nc5:v3d
nc6=netcdf([OGCM_prefix,OGCM_name(5).name]);    %nc4:t3d
nc7=netcdf([OGCM_prefix,OGCM_name(2).name]);    %nc5:s3d
%
% Interpole data on the OGCM Z grid and CROCO horizontal grid
%
%
% Read and extrapole the 2D variables
%
zeta=ext_data_OGCM(nc1,lonT,latT,'surf_el',tin,lon,lat,1,Roa,interp_method);
u2d=ext_data_OGCM(nc2,lonU,latU,'u_barotropic_velocity',tin,lon,lat,1,Roa,interp_method);
v2d=ext_data_OGCM(nc3,lonV,latV,'v_barotropic_velocity',tin,lon,lat,1,Roa,interp_method);
ubar=rho2u_2d(u2d.*cosa+v2d.*sina);
vbar=rho2v_2d(v2d.*cosa-u2d.*sina);
%
% Read and extrapole the 3D variables
%
NZ=length(Z);
[M,L]=size(lon);
dz=gradient(Z);
temp=zeros(NZ,M,L);
salt=zeros(NZ,M,L);
u=zeros(NZ,M,L-1);
v=zeros(NZ,M-1,L);
for k=1:NZ
  if rem(k,10)==0
    disp(['  Level ',num2str(k),' of ',num2str(NZ)])
  end
  u2d=ext_data_OGCM(nc4,lonU,latU,'water_u',tin,lon,lat,...
                    k,Roa,interp_method);
  v2d=ext_data_OGCM(nc5,lonV,latV,'water_v',tin,lon,lat,...
                    k,Roa,interp_method);
  u(k,:,:)=rho2u_2d(u2d.*cosa+v2d.*sina);
  v(k,:,:)=rho2v_2d(v2d.*cosa-u2d.*sina);
  temp(k,:,:)=ext_temp_OGCM(nc6,lonT,latT,'water_temp',tin,lon,lat,...
                            k,Roa,interp_method);
  salt(k,:,:)=ext_data_OGCM(nc7,lonT,latT,'salinity',tin,lon,lat,...
                            k,Roa,interp_method);
end
%
% Close the OGCM file
%
close(nc1);close(nc2);close(nc3);close(nc4)
close(nc5);close(nc6);close(nc7)
%
% Get the CROCO vertical grid
%
disp('  Vertical interpolations')
if ~isempty(nc_clm)
  theta_s=nc_clm{'theta_s'}(:);
  theta_b=nc_clm{'theta_b'}(:);
  hc=nc_clm{'hc'}(:);
  N=length(nc_clm('s_rho'));
  vtransform=nc_clm{'Vtransform'}(:);
    if  ~exist('vtransform')
        vtransform=1; %Old Vtransform
        disp([' NO VTRANSFORM parameter found'])
        disp([' USE TRANSFORM default value vtransform = 1'])
    end
end
if ~isempty(nc_bry)
  theta_s=nc_bry{'theta_s'}(:);
  theta_b=nc_bry{'theta_b'}(:);
  hc=nc_bry{'hc'}(:);
  N=length(nc_bry('s_rho'));
  vtransform=nc_bry{'Vtransform'}(:);
    if  ~exist('vtransform')
        vtransform=1; %Old Vtransform
        disp([' NO VTRANSFORM parameter found'])
        disp([' USE TRANSFORM default value vtransform = 1'])
    end
end
zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
zu=rho2u_3d(zr);
zv=rho2v_3d(zr);
zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
dzr=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=rho2u_3d(dzr);
dzv=rho2v_3d(dzr);
%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
Z=[100;Z;-100000];
u=cat(1,u(1,:,:),u);
u=cat(1,u,u(end,:,:));
v=cat(1,v(1,:,:),v);
v=cat(1,v,v(end,:,:));
temp=cat(1,temp(1,:,:),temp);
temp=cat(1,temp,temp(end,:,:));
salt=cat(1,salt,salt(end,:,:));
salt=cat(1,salt(1,:,:),salt);
% 
% Perform the vertical interpolations 
%
temp=ztosigma(flipdim(temp,1),zr,flipud(Z));
salt=ztosigma(flipdim(salt,1),zr,flipud(Z));
u=ztosigma(flipdim(u,1),zu,flipud(Z));
v=ztosigma(flipdim(v,1),zv,flipud(Z));
%
% Correct the horizontal transport 
% i.e. remove the interpolated tranport and add 
%      the OGCM transport
%
if conserv==1
  u=u-tridim(squeeze(sum(u.*dzu)./sum(dzu)),N);
  v=v-tridim(squeeze(sum(v.*dzv)./sum(dzv)),N);
  u=u+tridim(ubar,N);
  v=v+tridim(vbar,N);
end
%
% Barotropic velocities
%
ubar=squeeze(sum(u.*dzu)./sum(dzu));
vbar=squeeze(sum(v.*dzv)./sum(dzv));
%
%  fill the files
%
if ~isempty(nc_clm)
  nc_clm{'zeta'}(tout,:,:)=zeta;
  nc_clm{'SSH'}(tout,:,:)=zeta;
  nc_clm{'temp'}(tout,:,:,:)=temp;
  nc_clm{'salt'}(tout,:,:,:)=salt;
  nc_clm{'u'}(tout,:,:,:)=u;
  nc_clm{'v'}(tout,:,:,:)=v;
  nc_clm{'ubar'}(tout,:,:,:)=ubar;
  nc_clm{'vbar'}(tout,:,:,:)=vbar;
end
if ~isempty(nc_bry)
  for obcndx=1:4
    if obcndx==1
      nc_bry{'zeta_south'}(tout,:)=zeta(1,:);
      nc_bry{'temp_south'}(tout,:,:)=temp(:,1,:);
      nc_bry{'salt_south'}(tout,:,:)=salt(:,1,:);
      nc_bry{'u_south'}(tout,:,:)=u(:,1,:);
      nc_bry{'v_south'}(tout,:,:)=v(:,1,:);
      nc_bry{'ubar_south'}(tout,:,:)=ubar(1,:);
      nc_bry{'vbar_south'}(tout,:,:)=vbar(1,:);
   elseif obcndx==2
      nc_bry{'zeta_east'}(tout,:)=zeta(:,end);
      nc_bry{'temp_east'}(tout,:,:)=temp(:,:,end);
      nc_bry{'salt_east'}(tout,:,:)=salt(:,:,end);
      nc_bry{'u_east'}(tout,:,:)=u(:,:,end);
      nc_bry{'v_east'}(tout,:,:)=v(:,:,end);
      nc_bry{'ubar_east'}(tout,:,:)=ubar(:,end);
      nc_bry{'vbar_east'}(tout,:,:)=vbar(:,end);
    elseif obcndx==3
      nc_bry{'zeta_north'}(tout,:)=zeta(end,:);
      nc_bry{'temp_north'}(tout,:,:)=temp(:,end,:);
      nc_bry{'salt_north'}(tout,:,:)=salt(:,end,:);
      nc_bry{'u_north'}(tout,:,:)=u(:,end,:);
      nc_bry{'v_north'}(tout,:,:)=v(:,end,:);
      nc_bry{'ubar_north'}(tout,:,:)=ubar(end,:);
      nc_bry{'vbar_north'}(tout,:,:)=vbar(end,:);
    elseif obcndx==4
      nc_bry{'zeta_west'}(tout,:)=zeta(:,1);
      nc_bry{'temp_west'}(tout,:,:)=temp(:,:,1);
      nc_bry{'salt_west'}(tout,:,:)=salt(:,:,1);
      nc_bry{'u_west'}(tout,:,:)=u(:,:,1);
      nc_bry{'v_west'}(tout,:,:)=v(:,:,1);
      nc_bry{'ubar_west'}(tout,:,:)=ubar(:,1);
      nc_bry{'vbar_west'}(tout,:,:)=vbar(:,1);
    end
  end
end
