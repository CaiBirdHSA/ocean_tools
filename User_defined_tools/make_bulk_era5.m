clc;close all;clear all;
era5name='./../case/hato/era5_hato.nc';
bulkname='./../case/hato/roms_bulk.nc';
xtype= 'NC_FLOAT';

%%read att
timeunit=ncreadatt(era5name,'time','units');
stri    =strfind(timeunit,'since');
time0   =datenum(timeunit(stri+6:end));
dist    =(double(datenum('2000-01-01 00:00:00'))-time0)*24.d0*3600.d0;

%%read var
longitude = double(ncread(era5name,'longitude'));
latitude =  double(ncread(era5name,'latitude'));
time = double(ncread(era5name,'time'));
u10  = ncread(era5name,'u10');  %10 metre U wind component,m/s
v10  = ncread(era5name,'v10');  %10 metre V wind component,m/s
tcc  = ncread(era5name,'tcc');  %toal cloud cover 1
t2m  = ncread(era5name,'t2m');  %2 metre temperature, K
ssr  = ncread(era5name,'msnswrf');  %Mean surface net short-wave radiation flux, W/m2
sp   = ncread(era5name,'sp');   %surface_air_pressure, Pa
mtpr = ncread(era5name,'mtpr'); %Mean total precipitation rate, kg/m2/s
d2m  = ncread(era5name,'d2m');  %2 metre dewpoint temperature, K
lwrad_down = ncread(era5name,'msdwlwrf');  %Mean surface downward long-wave radiation flux, W/m2

nlon = length(longitude);
nlat = length(latitude);
tlen = length(time);
time = double(time)*3600.000d0-dist;  % convert hours to seconds
t2m  = t2m-273.150;  % convert Kelvin to Celsius

% RH=100*(EXP((17.625*TD)/(243.04+TD))/EXP((17.625*T)/(243.04+T))), TD is dewpoint temperature, T is temperature, both are in Celsius
rh   = 100*(exp((17.625*(d2m-273.15))./(243.04+(d2m-273.15)))./exp((17.625*t2m)./(243.04+t2m)));
rh(rh>100.0)=100.d0;
ssr(ssr<0) = 0.0001d0;
lwrad_down = double(lwrad_down);  % convert J/m2 to Watts/m2

%%create bulk nc
ncid = netcdf.create(bulkname,'64BIT_OFFSET');
lonx = netcdf.defDim(ncid,'longitude',nlon);    
laty = netcdf.defDim(ncid,'latitude',nlat);
timez= netcdf.defDim(ncid,'time',tlen);

timeid= netcdf.defVar(ncid,'time',xtype, [timez]);
lonid = netcdf.defVar(ncid,'longitude', xtype, [lonx]);
latid = netcdf.defVar(ncid,'latitude',  xtype, [laty]);
uid   = netcdf.defVar(ncid,'Uwind', xtype, [lonx laty timez]); %surface U-wind component, m/s
vid   = netcdf.defVar(ncid,'Vwind', xtype, [lonx laty timez]); %surface V-wind component, m/s
tid   = netcdf.defVar(ncid,'Tair', xtype, [lonx laty timez]);  %surface air temperature, Celsius
pid   = netcdf.defVar(ncid,'Pair', xtype, [lonx laty timez]);  %surface air pressure, mb
qid   = netcdf.defVar(ncid,'Qair', xtype, [lonx laty timez]);  %surface air relative humidity, %
cid   = netcdf.defVar(ncid,'cloud',xtype, [lonx laty timez]);  %cloud fraction, 0-1
rid   = netcdf.defVar(ncid,'rain', xtype, [lonx laty timez]);  %rain fall rate, kg/m2/s
swid  = netcdf.defVar(ncid,'swrad',xtype, [lonx laty timez]);  %net shortwave radiation flux, Watts/m2
lwid  = netcdf.defVar(ncid,'lwrad_down',xtype, [lonx laty timez]);  %net shortwave radiation flux, Watts/m2
%sstid = netcdf.defVar(ncid,'SST',  xtype, [lonx laty timez]);  %sea surface temperature, Celsius
%tdid  = netcdf.defVar(ncid,'Tdew', xtype, [lonx laty timez]);  %surface dew-point temperature, Kelvin

netcdf.putAtt(ncid,timeid,'long_name','Time');
netcdf.putAtt(ncid,timeid,'units','seconds since 2000-01-01 00:00:00.0');
%netcdf.putAtt(ncid,timeid,'calendar', 'gregorian');

netcdf.putAtt(ncid,lonid, 'long_name','longitude');
netcdf.putAtt(ncid,lonid, 'units', 'degrees_east');

netcdf.putAtt(ncid,latid, 'long_name','latitude');
netcdf.putAtt(ncid,latid, 'units', 'degrees_north');

netcdf.putAtt(ncid,uid, 'long_name','surface u-wind component');
netcdf.putAtt(ncid,uid, 'units', 'meter second-1');
netcdf.putAtt(ncid,uid, 'time',  'time');
netcdf.putAtt(ncid,uid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,vid, 'long_name','surface v-wind component');
netcdf.putAtt(ncid,vid, 'units', 'meter second-1');
netcdf.putAtt(ncid,vid, 'time',  'time');
netcdf.putAtt(ncid,vid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,tid, 'long_name','surface air temperature');
netcdf.putAtt(ncid,tid, 'units', 'Celsius');
netcdf.putAtt(ncid,tid, 'time',  'time');
netcdf.putAtt(ncid,tid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,pid, 'long_name','surface air pressure');
netcdf.putAtt(ncid,pid, 'units', 'Pa');
netcdf.putAtt(ncid,pid, 'time',  'time');
netcdf.putAtt(ncid,pid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,qid, 'long_name','surface air relative humidity');
netcdf.putAtt(ncid,qid, 'units', 'percentage');
netcdf.putAtt(ncid,qid, 'time',  'time');
netcdf.putAtt(ncid,qid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,cid, 'long_name','cloud fraction');
netcdf.putAtt(ncid,cid, 'units', '');
netcdf.putAtt(ncid,cid, 'time',  'time');
netcdf.putAtt(ncid,cid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,rid, 'long_name','rain fall rate');
netcdf.putAtt(ncid,rid, 'units', 'kilogram meter-2 second-1');
netcdf.putAtt(ncid,rid, 'time',  'time');
netcdf.putAtt(ncid,rid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,swid, 'long_name','solar shortwave radiation flux');
netcdf.putAtt(ncid,swid, 'units', 'watt meter-2');
netcdf.putAtt(ncid,swid, 'time',  'time');
netcdf.putAtt(ncid,swid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,lwid, 'long_name','downwelling longwave radiation flux');
netcdf.putAtt(ncid,lwid, 'units', 'watt meter-2');
netcdf.putAtt(ncid,lwid, 'time',  'time');
netcdf.putAtt(ncid,lwid, 'coordinates', 'longitude latitude');

netcdf.endDef(ncid)

if strcmp(xtype,'NC_FLOAT')
  disp('Is NC_FLOAT');
netcdf.putVar(ncid,timeid,single(time));
netcdf.putVar(ncid,lonid, single(longitude));
netcdf.putVar(ncid,latid, flip(single(latitude),1));
netcdf.putVar(ncid,uid,   flip(single(u10),2));
netcdf.putVar(ncid,vid,   flip(single(v10),2));
netcdf.putVar(ncid,tid,   flip(single(t2m),2));
netcdf.putVar(ncid,pid,   flip(single(sp),2));
netcdf.putVar(ncid,qid,   flip(single(rh),2));
netcdf.putVar(ncid,cid,   flip(single(tcc),2));
netcdf.putVar(ncid,rid,   flip(single(mtpr),2));
netcdf.putVar(ncid,swid,  flip(single(ssr),2));
netcdf.putVar(ncid,lwid,  flip(single(lwrad_down),2));
%netcdf.putVar(ncid,sstid, single(sst));
%netcdf.putVar(ncid,tdid,  single(d2m));
elseif strcmp(xtype,'NC_DOUBLE')
  disp('Is NC_DOUBLE');
netcdf.putVar(ncid,timeid,(time));
netcdf.putVar(ncid,lonid, (longitude));
netcdf.putVar(ncid,latid, flip(latitude,1));
netcdf.putVar(ncid,uid,   flip(u10,2));
netcdf.putVar(ncid,vid,   flip(v10,2));
netcdf.putVar(ncid,tid,   flip(t2m,2));
netcdf.putVar(ncid,pid,   flip(sp,2));
netcdf.putVar(ncid,qid,   flip(rh,2));
netcdf.putVar(ncid,cid,   flip(tcc,2));
netcdf.putVar(ncid,rid,   flip(mtpr,2));
netcdf.putVar(ncid,swid,  flip(ssr,2));
netcdf.putVar(ncid,lwid,  flip(lwrad_down,2));  
else
  disp('Error string for xtype');
end

netcdf.close(ncid)

