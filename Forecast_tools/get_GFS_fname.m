function fname=get_GFS_fname(time,gfs_run_time,gfstype)
%
%  fname=get_GFS_fname(time,gfs_run_time)
%
%  Give the GFS url for a given date.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% set URL
%
url='https://nomads.ncep.noaa.gov:9090';
%
% set file types
%
if gfstype==0
  gfsname ='fnl';
  gfsname1='fnlflx';     % 1/2 GDAS data
else
  gfsname ='gfs';
  gfsname1='gfs_0p25';  % 1/4 deg res GFS data
  %gfsname1='gfs_0p50';   % 1/2 deg res GFS data
end
%
% Get the date
%
[y,m,d,h,mi,s]=datevec(time);
stry=num2str(y);
if m<10
  strm=['0',num2str(m)];
else
  strm=num2str(m);
end
if d<10
  strd=['0',num2str(d)];
else
  strd=num2str(d);
end
if gfs_run_time < 10
  strh='_0';
else
  strh='_';
end
%
% Get the grid
%
if gfstype==0
  gfsdir =[url,'/dods/',gfsname,'/',gfsname,stry,strm,strd,'/'];
  fname=[gfsdir,gfsname1,strh,num2str(gfs_run_time),'z'];
else
  gfsdir =[url,'/dods/',gfsname1,'/',gfsname,stry,strm,strd,'/'];
  fname=[gfsdir,gfsname1,strh,num2str(gfs_run_time),'z'];
end
fname
