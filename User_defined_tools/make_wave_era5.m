clc;close all;clear all;
era5name='./../case/hato/era5_wave.nc';
bulkname='./../case/hato/roms_wave.nc';
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
dwave  = ncread(era5name,'mdww');  %Mean direction of wind waves, degrees
pwave  = ncread(era5name,'mpww');  %Mean period of wind waves, s
awave  = ncread(era5name,'shww');  %Significant height of wind waves, m

nlon = length(longitude);
nlat = length(latitude);
tlen = length(time);
time = double(time)*3600.000d0-dist;  % convert hours to seconds

%%create bulk nc
ncid = netcdf.create(bulkname,'64BIT_OFFSET');
lonx = netcdf.defDim(ncid,'longitude',nlon);    
laty = netcdf.defDim(ncid,'latitude',nlat);
timez= netcdf.defDim(ncid,'time',tlen);

timeid= netcdf.defVar(ncid,'time',xtype, [timez]);
lonid = netcdf.defVar(ncid,'longitude', xtype, [lonx]);
latid = netcdf.defVar(ncid,'latitude',  xtype, [laty]);
did   = netcdf.defVar(ncid,'Dwave', xtype, [lonx laty timez]); %surface U-wind component, m/s
pid   = netcdf.defVar(ncid,'Pwave', xtype, [lonx laty timez]); %surface V-wind component, m/s
hid   = netcdf.defVar(ncid,'Hwave', xtype, [lonx laty timez]);  %surface air temperature, Celsius

netcdf.putAtt(ncid,timeid,'long_name','Time');
netcdf.putAtt(ncid,timeid,'units','seconds since 2000-01-01 00:00:00');
%netcdf.putAtt(ncid,timeid,'calendar', 'gregorian');

netcdf.putAtt(ncid,lonid, 'long_name','longitude');
netcdf.putAtt(ncid,lonid, 'units', 'degrees_east');

netcdf.putAtt(ncid,latid, 'long_name','latitude');
netcdf.putAtt(ncid,latid, 'units', 'degrees_north');

netcdf.putAtt(ncid,did, 'long_name','wind induced wave direction');
netcdf.putAtt(ncid,did, 'units', 'degree');
netcdf.putAtt(ncid,did, 'time',  'time');
netcdf.putAtt(ncid,did, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,pid, 'long_name','wind induced wave period');
netcdf.putAtt(ncid,pid, 'units', 'second');
netcdf.putAtt(ncid,pid, 'time',  'time');
netcdf.putAtt(ncid,pid, 'coordinates', 'longitude latitude');

netcdf.putAtt(ncid,hid, 'long_name','wind-induced significant wave height');
netcdf.putAtt(ncid,hid, 'units', 'meter');
netcdf.putAtt(ncid,hid, 'time',  'time');
netcdf.putAtt(ncid,hid, 'coordinates', 'longitude latitude');

netcdf.endDef(ncid)

if strcmp(xtype,'NC_FLOAT')
  disp('Is NC_FLOAT');
netcdf.putVar(ncid,timeid,single(time));
netcdf.putVar(ncid,lonid, single(longitude));
netcdf.putVar(ncid,latid, flip(single(latitude),1));
netcdf.putVar(ncid,did,   flip(single(dwave),2));
netcdf.putVar(ncid,pid,   flip(single(pwave),2));
netcdf.putVar(ncid,hid,   flip(single(awave),2));
elseif strcmp(xtype,'NC_DOUBLE')
  disp('Is NC_DOUBLE');
netcdf.putVar(ncid,timeid,(time));
netcdf.putVar(ncid,lonid, (longitude));
netcdf.putVar(ncid,latid, flip(latitude,1));
netcdf.putVar(ncid,did,   flip(dwave,2));
netcdf.putVar(ncid,pid,   flip(pwave,2));
netcdf.putVar(ncid,hid,   flip(awave,2));
else
  disp('Error string for xtype');
end

netcdf.close(ncid)

