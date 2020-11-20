%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and fill frc and bulk files with GFS data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
tic
crocotools_param
%
makeplot = 0;
it=2;
%
frc_prefix=[frc_prefix,'_GFS_'];
blk_prefix=[blk_prefix,'_GFS_'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% time (in matlab time)
%
today=floor(now);
%
% date in 'Yorig' time
%
rundate_str='2020-09-09';
rundate=datenum(today)-datenum(Yorig,1,1);
%
% GFS data name
%
GFS_dir ='E:\croco_tools\DATA\GFS\2020090800\';
fcst_name=dir([GFS_dir,'gfs_',datestr(datenum(rundate_str)-1,'yyyymmddhh'),'*.nc']);
have_name=dir([GFS_dir,'gfs_3hravg_',datestr(datenum(rundate_str)-1,'yyyymmddhh'),'*.nc']);
%
%
if level==0
  nc_suffix='.nc';
else
  nc_suffix=['.nc.',num2str(level)];
  grdname=[grdname,'.',num2str(level)];
end
%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
close(nc)
cosa=cos(angle);
sina=sin(angle);
%
% Get the GFS grid 
% 
nc=netcdf([fcst_name(1).folder,'\', fcst_name(1).name]);
lon1=nc{'lon_0'}(:);
lat1=nc{'lat_0'}(:);
time=0:3:9;
tavg=sort([1.5:6:9,3:6:9-0.5]);
timeref=datestr(datenum(nc{'UGRD_P0_L103_GLL0'}.initial_time(:)),'yyyy-mm-dd HH:MM:SS');
timeunits=nc{'UGRD_P0_L103_GLL0'}.forecast_time_units(:);
% mask=nc{'mask'}(:);
tlen=length(time);
close(nc);
%
% bulk and forcing files
%
blkname=[blk_prefix,datestr(datenum(rundate_str)-1,'yyyymmdd'),nc_suffix];
disp(['Create a new bulk file: ' blkname])
create_forecast_bulk(blkname,grdname,CROCO_title,time,0);
nc_blk=netcdf(blkname,'write');
frcname=[frc_prefix,datestr(datenum(rundate_str)-1,'yyyymmdd'),nc_suffix];
disp(['Create a new forcing file: ' frcname])
create_forecast_forcing(frcname,grdname,CROCO_title,time,tavg)
nc_frc=netcdf(frcname,'write');
%
% Loop on time
%
missval=nan;
default=nan;
for l=1:tlen
  disp(['time index: ',num2str(l),' of total: ',num2str(tlen)])
  nc=netcdf([fcst_name(l).folder,'\', fcst_name(l).name]);
  var=squeeze(nc{'PRES_P0_L1_GLL0'}(:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_frc{'Pair'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    disp('Reading Pair error');
  end
  
  var=squeeze(nc{'TMP_P0_L103_GLL0'}(:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_frc{'Tair'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method)-273.15;
  else
    disp('Reading Tair error');
  end

  % Relative humidity at 2m
%   var=squeeze(nc{'RH_P0_L103_GLL0'}(:,:));
%   if mean(mean(isnan(var)~=1))
%     var=get_missing_val(lon1,lat1,var,missval,Roa,default);
%     nc_frc{'Qair'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
%   else
%     disp('Reading Qair(Relative humidity) error');
%   end
%   
   % Precipitation rate at surface
  var=squeeze(nc{'PRATE_P0_L1_GLL0'}(:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_frc{'rain'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    disp('Reading prate error');
  end
  
  %Zonal wind speed
  var=squeeze(nc{'UGRD_P0_L103_GLL0'}(:,:));
  if mean(mean(isnan(var)~=1))
    uwnd=get_missing_val(lon1,lat1,var,missval,Roa,default);
    uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
  else
    disp('Reading uwnd error');
  end
  
  %Meridian wind speed
  var=squeeze(nc{'VGRD_P0_L103_GLL0'}(:,:));
  if mean(mean(isnan(var)~=1))
    vwnd=get_missing_val(lon1,lat1,var,missval,Roa,default);
    vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
  else
    disp('Reading vwnd error');
  end
  
  nc_frc{'uwnd'}(l,:,:)=rho2u_2d(uwnd.*cosa+vwnd.*sina);
  nc_frc{'vwnd'}(l,:,:)=rho2v_2d(vwnd.*cosa-uwnd.*sina);
  
end
close(nc);  

for l=1:tlen-1
    disp(['time index: ',num2str(l),' of total: ',num2str(tlen)])
    nc=netcdf([have_name(l).folder,'\', have_name(l).name]);
    if l==1 || l==2
        % cloud fraction
        var=squeeze(nc{'TCDC_P8_L10_GLL0_avg'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'cloud'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
        else
            disp('Reading cloud error');
        end  
  
        %solar shortwave radiation flux
        var=squeeze(nc{'DSWRF_P8_L1_GLL0_avg'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'swrad'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
        else
            disp('Reading swrad error');
        end
  
        %Downward longwave flux
        var=squeeze(nc{'DLWRF_P8_L1_GLL0_avg'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'lwrad_down'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
        else
            disp('Reading lwrad_down error');
        end
    elseif mod(l,2)==1
        % cloud fraction
        var=squeeze(nc{'TCDC_P8_L10_GLL0_avg3h'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'cloud'}(l+1,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
        else
            disp('Reading cloud error');
        end  
  
        %solar shortwave radiation flux
        var=squeeze(nc{'DSWRF_P8_L1_GLL0_avg3h'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'swrad'}(l+1,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
        else
            disp('Reading swrad error');
        end
  
        %Downward longwave flux
        var=squeeze(nc{'DLWRF_P8_L1_GLL0_avg3h'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'lwrad_down'}(l+1,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
        else
            disp('Reading lwrad_down error');
        end
    else
        % cloud fraction
        var=squeeze(nc{'TCDC_P8_L10_GLL0_avg6h'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'cloud'}(l,:,:)=0.5d0*(interp2(lon1,lat1,var,lon,lat,interp_method)+nc_frc{'cloud'}(l-1,:,:));
        else
            disp('Reading cloud error');
        end  
  
        %solar shortwave radiation flux
        var=squeeze(nc{'DSWRF_P8_L1_GLL0_avg6h'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'swrad'}(l,:,:)=0.5d0*(interp2(lon1,lat1,var,lon,lat,interp_method)+nc_frc{'swrad'}(l-1,:,:));
        else
            disp('Reading swrad error');
        end
  
        %Downward longwave flux
        var=squeeze(nc{'DLWRF_P8_L1_GLL0_avg6h'}(:,:));
        if mean(mean(isnan(var)~=1))
            var=get_missing_val(lon1,lat1,var,missval,Roa,default);
            nc_frc{'lwrad_down'}(l,:,:)=0.5d0*(interp2(lon1,lat1,var,lon,lat,interp_method)+nc_frc{'lwrad_down'}(l-1,:,:));
        else
            disp('Reading lwrad_down error');
        end        
    end
        
end
% 
close(nc_frc);
close(nc_blk);
close(nc)

toc








