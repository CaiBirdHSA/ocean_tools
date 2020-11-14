%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start
tic
crocotools_param

%
% Get the lateral boundary conditions
%
make_OGCM_frcst
%
% Get the surface forcing
%
make_GFS
%
% Copy the resulting files
%
%rundate=datenum(date)-datenum(Yorig,1,1); nc_suffix='.nc';
eval(['!cp ',bry_prefix,num2str(rundate),nc_suffix,' ',bry_prefix,'0',nc_suffix])
eval(['!cp ',blk_prefix,num2str(rundate),nc_suffix,' ',blk_prefix,'0',nc_suffix])
eval(['!cp ',frc_prefix,num2str(rundate),nc_suffix,' ',frc_prefix,'0',nc_suffix])
eval(['!cp ',clm_prefix,num2str(rundate),nc_suffix,' ',clm_prefix,'0',nc_suffix])
eval(['!cp ',ini_prefix,num2str(rundate),nc_suffix,' ',ini_prefix,'0',nc_suffix])
%
% Add tidal data in forcing file
%
if add_tides_fcst==1
  disp(['Add tidal data ... '])
  frcname=[CROCO_files_dir,'croco_frc_GFS_0.nc'];
  [Y,M,d,h,mi,s] = datevec(date);
  add_tidal_data(tidename,grdname,frcname,Ntides,tidalrank,...
                                  Yorig,Y,M,coastfileplot)
end
%
%  Set the clock right: 
%  - copy croco_ini.nc in FORECAST/croco_ini.nc if not available
%  - update scrum_time in FORECAST/croco_ini.nc using timezone information
%    In this case, initial scrum_time is the UTC time corresponding to 
%    local midnight of day: now-hdays+1 (scrum_time needs to be a UTC time
%    since all forcing fields are referenced to UTC time).
%
disp('Set the clock right in initial file using Time Zone information')
time=(floor(now)-datenum(Yorig,1,1)-(hdays-1)-timezone/24)*86400;
ininame='FORECAST/croco_ini.nc';
nc=netcdf(ininame,'write');
if isempty(nc)
  disp('No restart file available in CROCO_FILES, copy OGCM inifile')
  eval(['!cp ',ini_prefix,num2str(rundate),nc_suffix,' ',ininame])
  nc=netcdf(ininame,'write');
end
nc{'scrum_time'}(:)=time;
close(nc)
toc
return
