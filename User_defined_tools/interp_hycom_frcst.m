function interp_hycom_frcst(hycom_name,Roa,interp_method,...
                           lonU,latU,lonV,latV,lonT,latT,Z,tin,...
		           nc_clm,nc_bry,lon,lat,angle,h,tout,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conserv=1; % same barotropic velocities as the OGCM
%
disp(['  Horizontal interpolation: ',hycom_name])
%
%
% CROCO grid angle
%
cosa=cos(angle);
sina=sin(angle);
%
% Open the OGCM file (1:el,2:s,3:ua,4:va,5:t,6:u,7:v)
%
nc=netcdf(hycom_name);    %nc1:el
%
% Interpole data on the OGCM Z grid and CROCO horizontal grid
%
%
% Read and extrapole the 2D variables
%
zeta=ext_data_hycom(nc,hycom_name,lonT,latT,'surf_el',tin,lon,lat,1,Roa,interp_method);
%
%bottom t,s,u,v
%
% s_bot = ext_data_hycom(nc,hycom_name,lonT,latT,'salinity_bottom',tin,lon,lat,1,Roa,interp_method);
% t_bot = ext_data_hycom(nc,hycom_name,lonT,latT,'water_temp_bottom',tin,lon,lat,1,Roa,interp_method);
% u_bot = ext_data_hycom(nc,hycom_name,lonT,latT,'water_u_bottom',tin,lon,lat,1,Roa,interp_method);
% v_bot = ext_data_hycom(nc,hycom_name,lonT,latT,'water_v_bottom',tin,lon,lat,1,Roa,interp_method);
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
  u2d=ext_data_hycom(nc,hycom_name,lonU,latU,'water_u',tin,lon,lat,...
                    k,Roa,interp_method);
  v2d=ext_data_hycom(nc,hycom_name,lonV,latV,'water_v',tin,lon,lat,...
                    k,Roa,interp_method);
  u(k,:,:)=rho2u_2d(u2d.*cosa+v2d.*sina);
  v(k,:,:)=rho2v_2d(v2d.*cosa-u2d.*sina);
  temp(k,:,:)=ext_data_hycom(nc,hycom_name,lonT,latT,'water_temp',tin,lon,lat,...
                            k,Roa,interp_method);
  salt(k,:,:)=ext_data_hycom(nc,hycom_name,lonT,latT,'salinity',tin,lon,lat,...
                            k,Roa,interp_method);
end

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
Z=[100;Z;-15000];
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
%existing u_barotropic velocity
u2d_flag=nc{'u_barotropic_velocity'}(:);
v2d_flag=nc{'v_barotropic_velocity'}(:);
if ~isempty(u2d_flag) && ~isempty(v2d_flag)
    u2d=ext_data_hycom(nc,hycom_name,lonU,latU,'u_barotropic_velocity',tin,lon,lat,1,Roa,interp_method);
    v2d=ext_data_hycom(nc,hycom_name,lonV,latV,'v_barotropic_velocity',tin,lon,lat,1,Roa,interp_method);
    ubar=rho2u_2d(u2d.*cosa+v2d.*sina);
    vbar=rho2v_2d(v2d.*cosa-u2d.*sina);
    if conserv==1
        u=u-tridim(squeeze(sum(u.*dzu)./sum(dzu)),N);
        v=v-tridim(squeeze(sum(v.*dzv)./sum(dzv)),N);
        u=u+tridim(ubar,N);
        v=v+tridim(vbar,N);
    end
else
    ubar=squeeze(sum(u.*dzu)./sum(dzu));
    vbar=squeeze(sum(v.*dzv)./sum(dzv));
end
%
% Close the OGCM file
%
close(nc);
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
