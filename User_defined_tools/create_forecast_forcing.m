function  create_forecast_forcing(frcname,grdname,title,smst,avgt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf forcing file 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc=netcdf(grdname,'r');
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
close(nc);
Lp=L+1;
Mp=M+1;

% nw = netcdf(frcname, 'clobber');
ncid=netcdf.create(frcname,'64BIT_OFFSET');
netcdf.close(ncid);
nw = netcdf(frcname,'write');
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
nw('time') = length(smst);
nw('cloud_time') = length(avgt);
%
%  Create variables and attributes
%
nw{'time'} = ncdouble('time');
nw{'time'}.long_name = ncchar('surface forcing time');
nw{'time'}.long_name = 'surface forcing time';
nw{'time'}.units = ncchar('hours');
nw{'time'}.units = 'hours';
nw{'time'}.field = ncchar('time, scalar, series');
nw{'time'}.field = 'time, scalar, series';

nw{'cloud_time'} = ncdouble('cloud_time');
nw{'cloud_time'}.long_name = ncchar('average surface forcing time');
nw{'cloud_time'}.long_name = 'average surface forcing time';
nw{'cloud_time'}.units = ncchar('hours');
nw{'cloud_time'}.units = 'hours';
nw{'cloud_time'}.field = ncchar('cloud_time, scalar, series');
nw{'cloud_time'}.field = 'cloud_time, scalar, series';

nw{'Pair'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Pair'}.long_name = ncchar('surface air pressure');
nw{'Pair'}.long_name = 'surface air pressure';
nw{'Pair'}.units = ncchar('Pa');
nw{'Pair'}.units = 'Pa';

nw{'Tair'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Tair'}.long_name = ncchar('surface air temperature');
nw{'Tair'}.long_name = 'surface air temperature';
nw{'Tair'}.units = ncchar('Celsius');
nw{'Tair'}.units = 'Celsius';

nw{'Qair'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Qair'}.long_name = ncchar('surface air relative humidity');
nw{'Qair'}.long_name = 'surface air relative humidity';
nw{'Qair'}.units = ncchar('kg/kg');
nw{'Qair'}.units = 'kg/kg';

nw{'rain'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'rain'}.long_name = ncchar('rain fall rate');
nw{'rain'}.long_name = 'rain fall rate';
nw{'rain'}.units = ncchar('kilogram meter-2 second-1');
nw{'rain'}.units = 'kilogram meter-2 second-1';

nw{'Uwnd'} = ncdouble('time', 'eta_u', 'xi_u');
nw{'Uwnd'}.long_name = ncchar('Xi-component of wind');
nw{'Uwnd'}.long_name = 'Xi-component of wind';
nw{'Uwnd'}.units = ncchar('meter second-1');
nw{'Uwnd'}.units = 'meter second-1';

nw{'Vwnd'} = ncdouble('time', 'eta_v', 'xi_v');
nw{'Vwnd'}.long_name = ncchar('Eta-component of wind');
nw{'Vwnd'}.long_name = 'Eta-component of wind';
nw{'Vwnd'}.units = ncchar('meter second-1');
nw{'Vwnd'}.units = 'meter second-1';

nw{'cloud'} = ncdouble('cloud_time', 'eta_rho', 'xi_rho');
nw{'cloud'}.long_name = ncchar('cloud fraction');
nw{'cloud'}.long_name = 'cloud fraction';
nw{'cloud'}.units = ncchar('');
nw{'cloud'}.units = '';
nw{'cloud'}.time = ncchar('cloud_time');
nw{'cloud'}.time = 'cloud_time';

nw{'swrad'} = ncdouble('cloud_time', 'eta_rho', 'xi_rho');
nw{'swrad'}.long_name = ncchar('solar shortwave radiation flux');
nw{'swrad'}.long_name = 'solar shortwave radiation flux';
nw{'swrad'}.units = ncchar('watts meter-2');
nw{'swrad'}.units = 'watts meter-2';
nw{'swrad'}.time = ncchar('cloud_time');
nw{'swrad'}.time = 'cloud_time';

nw{'lwrad_down'} = ncdouble('cloud_time', 'eta_rho', 'xi_rho');
nw{'lwrad_down'}.long_name = ncchar('downwelling longwave radiation flux');
nw{'lwrad_down'}.long_name = 'downwelling longwave radiation flux';
nw{'lwrad_down'}.units = ncchar('watts meter-2');
nw{'lwrad_down'}.units = 'watts meter-2';
nw{'lwrad_down'}.time = ncchar('cloud_time');
nw{'lwrad_down'}.time = 'cloud_time';

nw{'Awave'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Awave'}.long_name = ncchar('wind induced wave amplitude');
nw{'Awave'}.long_name = 'wind induced wave amplitude';
nw{'Awave'}.units = ncchar('m');
nw{'Awave'}.units = 'm';

nw{'Dwave'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Dwave'}.long_name = ncchar('wind induced wave direction');
nw{'Dwave'}.long_name = 'wind induced wave direction';
nw{'Dwave'}.units = ncchar('degree');
nw{'Dwave'}.units = 'degree';

nw{'Pwave'} = ncdouble('time', 'eta_rho', 'xi_rho');
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
nw.type = ncchar('suface forcing file');
nw.type = 'surface forcing file';

%
% Write time variables
%

nw{'time'}(:) = smst;
nw{'cloud_time'}(:)=avgt;

close(nw);
