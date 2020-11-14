function  create_forcing(frcname,grdname,title,smst,Ymin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf forcing file
%       frcname: name of the forcing file
%       grdname: name of the grid file
%       title: title in the netcdf file    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc=netcdf(grdname,'r');
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
close(nc);
Lp=L+1;
Mp=M+1;

ncid=netcdf.create(frcname,'64BIT_OFFSET');
netcdf.close(ncid);
nw = netcdf(frcname, 'clobber');
%result = redef(nw);

%
%  Create dimensions
%

nw('xi_u') = L;
nw('eta_u') = Mp;
nw('xi_v') = Lp;
nw('eta_v') = M;
nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('xi_psi') = L;
nw('eta_psi') = M;
nw('ocean_time') = 48;
%
%  Create variables and attributes
%
nw{'ocean_time'} = ncdouble('ocean_time');
nw{'ocean_time'}.long_name = ncchar('surface momentum stress time');
nw{'ocean_time'}.long_name = 'surface momentum stress time';
nw{'ocean_time'}.units = ncchar('seconds since 1980-01-01 00:00:00');
nw{'ocean_time'}.units = 'seconds since 1980-01-01 00:00:00';
nc{'ocean_time'}.calendar = ncchar('gregorian');
nc{'ocean_time'}.calendar = 'gregorian';
nc{'ocean_time'}.field = 'time, scalar, series';

nw{'sustr'} = ncdouble('ocean_time', 'eta_u', 'xi_u');
nw{'sustr'}.long_name = ncchar('surface u-momentum stress');
nw{'sustr'}.long_name = 'surface u-momentum stress';
nw{'sustr'}.units = ncchar('Newton meter-2');
nw{'sustr'}.units = 'Newton meter-2';

nw{'svstr'} = ncdouble('ocean_time', 'eta_v', 'xi_v');
nw{'svstr'}.long_name = ncchar('surface v-momentum stress');
nw{'svstr'}.long_name = 'surface v-momentum stress';
nw{'svstr'}.units = ncchar('Newton meter-2');
nw{'svstr'}.units = 'Newton meter-2';

nw{'shflux'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'shflux'}.long_name = ncchar('surface net heat flux');
nw{'shflux'}.long_name = 'surface net heat flux';
nw{'shflux'}.units = ncchar('Watts meter-2');
nw{'shflux'}.units = 'Watts meter-2';

nw{'swflux'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'swflux'}.long_name = ncchar('surface freshwater flux (E-P)');
nw{'swflux'}.long_name = 'surface freshwater flux (E-P)';
nw{'swflux'}.units = ncchar('centimeter day-1');
nw{'swflux'}.units = 'centimeter day-1';
nw{'swflux'}.positive = ncchar('net evaporation');
nw{'swflux'}.positive = 'net evaporation';
nw{'swflux'}.negative = ncchar('net precipitation');
nw{'swflux'}.negative = 'net precipitation';

nw{'SST'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'SST'}.long_name = ncchar('sea surface temperature');
nw{'SST'}.long_name = 'sea surface temperature';
nw{'SST'}.units = ncchar('Celsius');
nw{'SST'}.units = 'Celsius';

nw{'SSS'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'SSS'}.long_name = ncchar('sea surface salinity');
nw{'SSS'}.long_name = 'sea surface salinity';
nw{'SSS'}.units = ncchar('PSU');
nw{'SSS'}.units = 'PSU';

nw{'dQdSST'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'dQdSST'}.long_name = ncchar('surface net heat flux sensitivity to SST');
nw{'dQdSST'}.long_name = 'surface net heat flux sensitivity to SST';
nw{'dQdSST'}.units = ncchar('Watts meter-2 Celsius-1');
nw{'dQdSST'}.units = 'Watts meter-2 Celsius-1';

nw{'swrad'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'swrad'}.long_name = ncchar('solar shortwave radiation');
nw{'swrad'}.long_name = 'solar shortwave radiation';
nw{'swrad'}.units = ncchar('Watts meter-2');
nw{'swrad'}.units = 'Watts meter-2';
nw{'swrad'}.positive = ncchar('downward flux, heating');
nw{'swrad'}.positive = 'downward flux, heating';
nw{'swrad'}.negative = ncchar('upward flux, cooling');
nw{'swrad'}.negative = 'upward flux, cooling';

nw{'Awave'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'Awave'}.long_name = ncchar('wind induced wave amplitude');
nw{'Awave'}.long_name = 'wind induced wave amplitude';
nw{'Awave'}.units = ncchar('m');
nw{'Awave'}.units = 'm';

nw{'Dwave'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'Dwave'}.long_name = ncchar('wind induced wave direction');
nw{'Dwave'}.long_name = 'wind induced wave direction';
nw{'Dwave'}.units = ncchar('degree');
nw{'Dwave'}.units = 'degree';

nw{'Pwave'} = ncdouble('ocean_time', 'eta_rho', 'xi_rho');
nw{'Pwave'}.long_name = ncchar('wind induced wave period');
nw{'Pwave'}.long_name = 'wind induced wave period';
nw{'Pwave'}.units = ncchar('second');
nw{'Pwave'}.units = 'second';

%result = endef(nw);

%
% Create global attributes
%

nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.grd_file = ncchar(grdname);
nw.grd_file = grdname;
nw.type = ncchar('CROCO forcing file');
nw.type = 'CROCO forcing file';

%
% Write time variables
%

nw{'ocean_time'}(:) = (datenum(Ymin-1,1,1)-datenum(1980,1,1)+[15:30:365,380:30:365*2,745:30:365*3,365*3+15:30:365*4])*24*3600;

close(nw);
