function create_forecast_bryfile(bryname,grdname,title,obc,...
                        theta_s,theta_b,hc,N,...
                        time,timeunits,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',bryname])
disp(' ')
if nargin < 10
    disp([' NO VTRANSFORM parameter found'])
    disp([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1; 
end
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file and check the topography
%
nc = netcdf(grdname, 'nowrite');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
close(nc);
hmin=min(min(h(maskr==1)));
if vtransform ==1
  if hc > hmin
    error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
  end
end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the boundary file
%
type = 'BOUNDARY file' ; 
history = 'HYCOM' ;
ncid=netcdf.create(bryname,'64BIT_OFFSET');
netcdf.close(ncid);
nc = netcdf(bryname,'write');
%%result = redef(nc);
%
%  Create dimensions
%
nc('xi_u') = L;
nc('xi_v') = Lp;
nc('xi_rho') = Lp;
nc('eta_u') = Mp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('s_rho') = N;
nc('s_w') = Np;
nc('tracer') = 2;
nc('ocean_time') = length(time);
nc('one') = 1;
%
%  Create variables and attributes
%
nc{'spherical'} = ncchar('one') ;
nc{'spherical'}.long_name = ncchar('grid type logical switch');
nc{'spherical'}.long_name = 'grid type logical switch';
nc{'spherical'}.flag_values = ncchar('T, F');
nc{'spherical'}.flag_values = 'T, F';
nc{'spherical'}.flag_meanings = ncchar('spherical Cartesian');
nc{'spherical'}.flag_meanings = 'spherical Cartesian';
%
nc{'Vtransform'} = ncint('one') ;
nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
nc{'Vstretching'} = ncint('one') ;
nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = ncchar('nondimensional');
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'} = ncdouble('one') ;
nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = ncchar('nondimensional');
nc{'theta_b'}.units = 'nondimensional';
%
nc{'Tcline'} = ncdouble('one') ;
nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nc{'Tcline'}.units = ncchar('meter');
nc{'Tcline'}.units = 'meter';
%
nc{'hc'} = ncdouble('one') ;
nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = ncchar('meter');
nc{'hc'}.units = 'meter';
%
nc{'sc_r'} = ncdouble('s_rho') ;
nc{'sc_r'}.long_name = ncchar('S-coordinate at RHO-points');
nc{'sc_r'}.long_name = 'S-coordinate at RHO-points';
nc{'sc_r'}.valid_min = -1.;
nc{'sc_r'}.valid_max = 0.;
nc{'sc_r'}.positive = ncchar('up');
nc{'sc_r'}.positive = 'up';
if (vtransform == 1)
    nc{'sc_r'}.standard_name = ncchar('ocena_s_coordinate_g1');
    nc{'sc_r'}.standard_name = 'ocena_s_coordinate_g1';
elseif (vtransform == 2)
    nc{'sc_r'}.standard_name = ncchar('ocena_s_coordinate_g2');
    nc{'sc_r'}.standard_name = 'ocena_s_coordinate_g2';
end
nc{'sc_r'}.formula_terms = ncchar('s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
nc{'sc_r'}.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
%
nc{'sc_w'} = ncdouble('s_w') ;
nc{'sc_w'}.long_name = ncchar('S-coordinate at W-points');
nc{'sc_w'}.long_name = 'S-coordinate at W-points';
nc{'sc_w'}.valid_min = -1. ;
nc{'sc_w'}.valid_max = 0. ;
nc{'sc_w'}.positive = ncchar('up');
nc{'sc_w'}.positive = 'up';
if (vtransform == 1)
    nc{'sc_w'}.standard_name = ncchar('ocena_s_coordinate_g1');
    nc{'sc_w'}.standard_name = 'ocena_s_coordinate_g1';
elseif (vtransform == 2)
    nc{'sc_w'}.standard_name = ncchar('ocena_s_coordinate_g2');
    nc{'sc_w'}.standard_name = 'ocena_s_coordinate_g2';
end
nc{'sc_w'}.formula_terms = ncchar('s: s_w C: Cs_w eta: zeta depth: h depth_c: hc');
nc{'sc_w'}.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
%
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_r'}.units = ncchar('nondimensional');
nc{'Cs_r'}.units = 'nondimensional';
nc{'Cs_r'}.valid_min = -1;
nc{'Cs_r'}.valid_max = 0;
%
nc{'Cs_w'} = ncdouble('s_w') ;
nc{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc{'Cs_w'}.long_name = 'S-coordinate stretching curves at W-points';
nc{'Cs_w'}.units = ncchar('nondimensional');
nc{'Cs_w'}.units = 'nondimensional';
nc{'Cs_w'}.valid_min = -1;
nc{'Cs_w'}.valid_max = 0;
%
nc{'ocean_time'} = ncdouble('ocean_time') ;
nc{'ocean_time'}.long_name = ncchar('time for boundary');
nc{'ocean_time'}.long_name = 'time for boundary';
nc{'ocean_time'}.units = ncchar(timeunits);
nc{'ocean_time'}.units = timeunits;
nc{'ocean_time'}.calendar = ncchar('standard');
nc{'ocean_time'}.calendar = 'standard';
%
if obc(1)==1
%
%   Southern boundary
%
  nc{'temp_south'} = ncdouble('ocean_time','s_rho','xi_rho') ;
  nc{'temp_south'}.long_name = ncchar('southern boundary potential temperature');
  nc{'temp_south'}.long_name = 'southern boundary potential temperature';
  nc{'temp_south'}.units = ncchar('Celsius');
  nc{'temp_south'}.units = 'Celsius';
  nc{'temp_south'}.coordinates = ncchar('lon_rho s_rho ocean_time');
  nc{'temp_south'}.coordinates = 'lon_rho s_rho ocean_time';
  nc{'temp_south'}.time = ncchar('ocean_time');
  nc{'temp_south'}.time = 'ocean_time';
%
  nc{'salt_south'} = ncdouble('ocean_time','s_rho','xi_rho') ;
  nc{'salt_south'}.long_name = ncchar('southern boundary salinity');
  nc{'salt_south'}.long_name = 'southern boundary salinity';
  nc{'salt_south'}.units = ncchar('PSU');
  nc{'salt_south'}.units = 'PSU';
  nc{'salt_south'}.coordinates = ncchar('lon_rho s_rho ocean_time');
  nc{'salt_south'}.coordinates = 'lon_rho s_rho ocean_time';
  nc{'salt_south'}.time = ncchar('ocean_time');
  nc{'salt_south'}.time = 'ocean_time';
%
  nc{'u_south'} = ncdouble('ocean_time','s_rho','xi_u') ;
  nc{'u_south'}.long_name = ncchar('southern boundary u-momentum component');
  nc{'u_south'}.long_name = 'southern boundary u-momentum component';
  nc{'u_south'}.units = ncchar('meter second-1');
  nc{'u_south'}.units = 'meter second-1';
  nc{'u_south'}.coordinates = ncchar('lon_u s_rho ocean_time');
  nc{'u_south'}.coordinates = 'lon_u s_rho ocean_time';
  nc{'u_south'}.time = ncchar('ocean_time');
  nc{'u_south'}.time = 'ocean_time';
%
  nc{'v_south'} = ncdouble('ocean_time','s_rho','xi_rho') ;
  nc{'v_south'}.long_name = ncchar('southern boundary v-momentum component');
  nc{'v_south'}.long_name = 'southern boundary v-momentum component';
  nc{'v_south'}.units = ncchar('meter second-1');
  nc{'v_south'}.units = 'meter second-1';
  nc{'v_south'}.coordinates = ncchar('lon_v s_rho ocean_time');
  nc{'v_south'}.coordinates = 'lon_v s_rho ocean_time';
  nc{'v_south'}.time = ncchar('ocean_time');
  nc{'v_south'}.time = 'ocean_time';  
%
  nc{'ubar_south'} = ncdouble('ocean_time','xi_u') ;
  nc{'ubar_south'}.long_name = ncchar('southern boundary vertically integrated u-momentum component');
  nc{'ubar_south'}.long_name = 'southern boundary vertically integrated u-momentum component';
  nc{'ubar_south'}.units = ncchar('meter second-1');
  nc{'ubar_south'}.units = 'meter second-1';
  nc{'ubar_south'}.coordinates = ncchar('lon_u ocean_time');
  nc{'ubar_south'}.coordinates = 'lon_u ocean_time';
  nc{'ubar_south'}.time = ncchar('ocean_time');
  nc{'ubar_south'}.time = 'ocean_time';
%
  nc{'vbar_south'} = ncdouble('ocean_time','xi_rho') ;
  nc{'vbar_south'}.long_name = ncchar('southern boundary vertically integrated v-momentum component');
  nc{'vbar_south'}.long_name = 'southern boundary vertically integrated v-momentum component';
  nc{'vbar_south'}.units = ncchar('meter second-1');
  nc{'vbar_south'}.units = 'meter second-1';
  nc{'vbar_south'}.coordinates = ncchar('lon_v ocean_time');
  nc{'vbar_south'}.coordinates = 'lon_v ocean_time';
  nc{'vbar_south'}.time = ncchar('ocean_time');
  nc{'vbar_south'}.time = 'ocean_time';
%
  nc{'zeta_south'} = ncdouble('ocean_time','xi_rho') ;
  nc{'zeta_south'}.long_name = ncchar('southern boundary sea surface height');
  nc{'zeta_south'}.long_name = 'southern boundary sea surface height';
  nc{'zeta_south'}.units = ncchar('meter');
  nc{'zeta_south'}.units = 'meter';
  nc{'zeta_south'}.coordinates = ncchar('lon_rho ocean_time');
  nc{'zeta_south'}.coordinates = 'lon_rho ocean_time';
  nc{'zeta_south'}.time = ncchar('ocean_time');
  nc{'zeta_south'}.time = 'ocean_time';
%
end
%
if obc(2)==1
%
%   Eastern boundary
%
  nc{'temp_east'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'temp_east'}.long_name = ncchar('eastern boundary potential temperature');
  nc{'temp_east'}.long_name = 'eastern boundary potential temperature';
  nc{'temp_east'}.units = ncchar('Celsius');
  nc{'temp_east'}.units = 'Celsius';
  nc{'temp_east'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'temp_east'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'temp_east'}.time = ncchar('ocean_time');
  nc{'temp_east'}.time = 'ocean_time';
%
  nc{'salt_east'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'salt_east'}.long_name = ncchar('eastern boundary salinity');
  nc{'salt_east'}.long_name = 'eastern boundary salinity';
  nc{'salt_east'}.units = ncchar('PSU');
  nc{'salt_east'}.units = 'PSU';
  nc{'salt_east'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'salt_east'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'salt_east'}.time = ncchar('ocean_time');
  nc{'salt_east'}.time = 'ocean_time';
  %
  nc{'mud_east_01'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'mud_east_01'}.long_name = ncchar('eastern boundary mud_01');
  nc{'mud_east_01'}.long_name = 'eastern boundary mud_01';
  nc{'mud_east_01'}.units = ncchar('kilogram meter-3');
  nc{'mud_east_01'}.units = 'kilogram meter-3';
  nc{'mud_east_01'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'mud_east_01'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'mud_east_01'}.time = ncchar('ocean_time');
  nc{'mud_east_01'}.time = 'ocean_time';
  %
  nc{'sand_east_01'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'sand_east_01'}.long_name = ncchar('eastern boundary sand_01');
  nc{'sand_east_01'}.long_name = 'eastern boundary sand_01';
  nc{'sand_east_01'}.units = ncchar('kilogram meter-3');
  nc{'sand_east_01'}.units = 'kilogram meter-3';
  nc{'sand_east_01'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'sand_east_01'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'sand_east_01'}.time = ncchar('ocean_time');
  nc{'sand_east_01'}.time = 'ocean_time';
  %
  nc{'sand_east_02'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'sand_east_02'}.long_name = ncchar('eastern boundary sand_02');
  nc{'sand_east_02'}.long_name = 'eastern boundary sand_02';
  nc{'sand_east_02'}.units = ncchar('kilogram meter-3');
  nc{'sand_east_02'}.units = 'kilogram meter-3';
  nc{'sand_east_02'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'sand_east_02'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'sand_east_02'}.time = ncchar('ocean_time');
  nc{'sand_east_02'}.time = 'ocean_time';
%
  nc{'u_east'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'u_east'}.long_name = ncchar('eastern boundary u-momentum component');
  nc{'u_east'}.long_name = 'eastern boundary u-momentum component';
  nc{'u_east'}.units = ncchar('meter second-1');
  nc{'u_east'}.units = 'meter second-1';
  nc{'u_east'}.coordinates = ncchar('lat_u s_rho ocean_time');
  nc{'u_east'}.coordinates = 'lat_u s_rho ocean_time';
  nc{'u_east'}.time = ncchar('ocean_time');
  nc{'u_east'}.time = 'ocean_time';
%
  nc{'v_east'} = ncdouble('ocean_time','s_rho','eta_v') ;
  nc{'v_east'}.long_name = ncchar('eastern boundary v-momentum component');
  nc{'v_east'}.long_name = 'eastern boundary v-momentum component';
  nc{'v_east'}.units = ncchar('meter second-1');
  nc{'v_east'}.units = 'meter second-1';
  nc{'v_east'}.coordinates = ncchar('lat_v s_rho ocean_time');
  nc{'v_east'}.coordinates = 'lat_v s_rho ocean_time';
  nc{'v_east'}.time = ncchar('ocean_time');
  nc{'v_east'}.time = 'ocean_time';
%
  nc{'ubar_east'} = ncdouble('ocean_time','eta_rho') ;
  nc{'ubar_east'}.long_name = ncchar('eastern boundary vertically integrated u-momentum component');
  nc{'ubar_east'}.long_name = 'eastern boundary vertically integrated u-momentum component';
  nc{'ubar_east'}.units = ncchar('meter second-1');
  nc{'ubar_east'}.units = 'meter second-1';
  nc{'ubar_east'}.coordinates = ncchar('lat_u ocean_time');
  nc{'ubar_east'}.coordinates = 'lat_u ocean_time';
  nc{'ubar_east'}.time = ncchar('ocean_time');
  nc{'ubar_east'}.time = 'ocean_time';
%
  nc{'vbar_east'} = ncdouble('ocean_time','eta_v') ;
  nc{'vbar_east'}.long_name = ncchar('eastern boundary vertically integrated v-momentum component');
  nc{'vbar_east'}.long_name = 'eastern boundary vertically integrated v-momentum component';
  nc{'vbar_east'}.units = ncchar('meter second-1');
  nc{'vbar_east'}.units = 'meter second-1';
  nc{'vbar_east'}.coordinates = ncchar('lat_v ocean_time');
  nc{'vbar_east'}.coordinates = 'lat_v ocean_time';
  nc{'vbar_east'}.time = ncchar('ocean_time');
  nc{'vbar_east'}.time = 'ocean_time';
%
  nc{'zeta_east'} = ncdouble('ocean_time','eta_rho') ;
  nc{'zeta_east'}.long_name = ncchar('eastern boundary sea surface height');
  nc{'zeta_east'}.long_name = 'eastern boundary sea surface height';
  nc{'zeta_east'}.units = ncchar('meter');
  nc{'zeta_east'}.units = 'meter';
  nc{'zeta_east'}.coordinates = ncchar('lat_rho ocean_time');
  nc{'zeta_east'}.coordinates = 'lat_rho ocean_time';
  nc{'zeta_east'}.time = ncchar('ocean_time');
  nc{'zeta_east'}.time = 'ocean_time';
%
end
%
if obc(3)==1
%
%   Northern boundary
%
  nc{'temp_north'} = ncdouble('ocean_time','s_rho','xi_rho') ;
  nc{'temp_north'}.long_name = ncchar('northern boundary potential temperature');
  nc{'temp_north'}.long_name = 'northern boundary potential temperature';
  nc{'temp_north'}.units = ncchar('Celsius');
  nc{'temp_north'}.units = 'Celsius';
  nc{'temp_north'}.coordinates = ncchar('lon_rho s_rho ocean_time');
  nc{'temp_north'}.coordinates = 'lon_rho s_rho ocean_time';
  nc{'temp_north'}.time = ncchar('ocean_time');
  nc{'temp_north'}.time = 'ocean_time';
%
  nc{'salt_north'} = ncdouble('ocean_time','s_rho','xi_rho') ;
  nc{'salt_north'}.long_name = ncchar('northern boundary salinity');
  nc{'salt_north'}.long_name = 'northern boundary salinity';
  nc{'salt_north'}.units = ncchar('PSU');
  nc{'salt_north'}.units = 'PSU';
  nc{'salt_north'}.coordinates = ncchar('lon_rho s_rho ocean_time');
  nc{'salt_north'}.coordinates = 'lon_rho s_rho ocean_time';
  nc{'salt_north'}.time = ncchar('ocean_time');
  nc{'salt_north'}.time = 'ocean_time';
%
  nc{'u_north'} = ncdouble('ocean_time','s_rho','xi_u') ;
  nc{'u_north'}.long_name = ncchar('northern boundary u-momentum component');
  nc{'u_north'}.long_name = 'northern boundary u-momentum component';
  nc{'u_north'}.units = ncchar('meter second-1');
  nc{'u_north'}.units = 'meter second-1';
  nc{'u_north'}.coordinates = ncchar('lon_u s_rho ocean_time');
  nc{'u_north'}.coordinates = 'lon_u s_rho ocean_time';
  nc{'u_north'}.time = ncchar('ocean_time');
  nc{'u_north'}.time = 'ocean_time';
%
  nc{'v_north'} = ncdouble('ocean_time','s_rho','xi_rho') ;
  nc{'v_north'}.long_name = ncchar('northern boundary v-momentum component');
  nc{'v_north'}.long_name = 'northern boundary v-momentum component';
  nc{'v_north'}.units = ncchar('meter second-1');
  nc{'v_north'}.units = 'meter second-1';
  nc{'v_north'}.coordinates = ncchar('lon_v s_rho ocean_time');
  nc{'v_north'}.coordinates = 'lon_v s_rho ocean_time';
  nc{'v_north'}.time = ncchar('ocean_time');
  nc{'v_north'}.time = 'ocean_time';
%
  nc{'ubar_north'} = ncdouble('ocean_time','xi_u') ;
  nc{'ubar_north'}.long_name = ncchar('northern boundary vertically integrated u-momentum component');
  nc{'ubar_north'}.long_name = 'northern boundary vertically integrated u-momentum component';
  nc{'ubar_north'}.units = ncchar('meter second-1');
  nc{'ubar_north'}.units = 'meter second-1';
  nc{'ubar_north'}.coordinates = ncchar('lon_u ocean_time');
  nc{'ubar_north'}.coordinates = 'lon_u ocean_time';
  nc{'ubar_north'}.time = ncchar('ocean_time');
  nc{'ubar_north'}.time = 'ocean_time';
%
  nc{'vbar_north'} = ncdouble('ocean_time','xi_rho') ;
  nc{'vbar_north'}.long_name = ncchar('northern boundary vertically integrated v-momentum component');
  nc{'vbar_north'}.long_name = 'northern boundary vertically integrated v-momentum component';
  nc{'vbar_north'}.units = ncchar('meter second-1');
  nc{'vbar_north'}.units = 'meter second-1';
  nc{'vbar_north'}.coordinates = ncchar('lon_v ocean_time');
  nc{'vbar_north'}.coordinates = 'lon_v ocean_time';
  nc{'vbar_north'}.time = ncchar('ocean_time');
  nc{'vbar_north'}.time = 'ocean_time';

  nc{'zeta_north'} = ncdouble('ocean_time','xi_rho') ;
  nc{'zeta_north'}.long_name = ncchar('northern boundary sea surface height');
  nc{'zeta_north'}.long_name = 'northern boundary sea surface height';
  nc{'zeta_north'}.units = ncchar('meter');
  nc{'zeta_north'}.units = 'meter';
  nc{'zeta_north'}.coordinates = ncchar('lon_rho ocean_time');
  nc{'zeta_north'}.coordinates = 'lon_rho ocean_time';
  nc{'zeta_north'}.time = ncchar('ocean_time');
  nc{'zeta_north'}.time = 'ocean_time';
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc{'temp_west'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'temp_west'}.long_name = ncchar('western boundary potential temperature');
  nc{'temp_west'}.long_name = 'western boundary potential temperature';
  nc{'temp_west'}.units = ncchar('Celsius');
  nc{'temp_west'}.units = 'Celsius';
  nc{'temp_west'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'temp_west'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'temp_west'}.time = ncchar('ocean_time');
  nc{'temp_west'}.time = 'ocean_time';
%
  nc{'salt_west'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'salt_west'}.long_name = ncchar('western boundary salinity');
  nc{'salt_west'}.long_name = 'western boundary salinity';
  nc{'salt_west'}.units = ncchar('PSU');
  nc{'salt_west'}.units = 'PSU';
  nc{'salt_west'}.coordinates = ncchar('lat_rho s_rho ocean_time');
  nc{'salt_west'}.coordinates = 'lat_rho s_rho ocean_time';
  nc{'salt_west'}.time = ncchar('ocean_time');
  nc{'salt_west'}.time = 'ocean_time';
%
  nc{'u_west'} = ncdouble('ocean_time','s_rho','eta_rho') ;
  nc{'u_west'}.long_name = ncchar('western boundary u-momentum component');
  nc{'u_west'}.long_name = 'western boundary u-momentum component';
  nc{'u_west'}.units = ncchar('meter second-1');
  nc{'u_west'}.units = 'meter second-1';
  nc{'u_west'}.coordinates = ncchar('lat_u s_rho ocean_time');
  nc{'u_west'}.coordinates = 'lat_u s_rho ocean_time';
  nc{'u_west'}.time = ncchar('ocean_time');
  nc{'u_west'}.time = 'ocean_time';
%
  nc{'v_west'} = ncdouble('ocean_time','s_rho','eta_v') ;
  nc{'v_west'}.long_name = ncchar('western boundary v-momentum component');
  nc{'v_west'}.long_name = 'western boundary v-momentum component';
  nc{'v_west'}.units = ncchar('meter second-1');
  nc{'v_west'}.units = 'meter second-1';
  nc{'v_west'}.coordinates = ncchar('lat_v s_rho ocean_time');
  nc{'v_west'}.coordinates = 'lat_v s_rho ocean_time';
  nc{'v_west'}.time = ncchar('ocean_time');
  nc{'v_west'}.time = 'ocean_time';
%
  nc{'ubar_west'} = ncdouble('ocean_time','eta_rho') ;
  nc{'ubar_west'}.long_name = ncchar('western boundary vertically integrated u-momentum component');
  nc{'ubar_west'}.long_name = 'western boundary vertically integrated u-momentum component';
  nc{'ubar_west'}.units = ncchar('meter second-1');
  nc{'ubar_west'}.units = 'meter second-1';
  nc{'ubar_west'}.coordinates = ncchar('lat_u ocean_time');
  nc{'ubar_west'}.coordinates = 'lat_u ocean_time';
  nc{'ubar_west'}.time = ncchar('ocean_time');
  nc{'ubar_west'}.time = 'ocean_time';
%
  nc{'vbar_west'} = ncdouble('ocean_time','eta_v') ;
  nc{'vbar_west'}.long_name = ncchar('western boundary vertically integrated v-momentum component');
  nc{'vbar_west'}.long_name = 'western boundary vertically integrated v-momentum component';
  nc{'vbar_west'}.units = ncchar('meter second-1');
  nc{'vbar_west'}.units = 'meter second-1';
  nc{'vbar_west'}.coordinates = ncchar('lat_v ocean_time');
  nc{'vbar_west'}.coordinates = 'lat_v ocean_time';
  nc{'vbar_west'}.time = ncchar('ocean_time');
  nc{'vbar_west'}.time = 'ocean_time';
%
  nc{'zeta_west'} = ncdouble('ocean_time','eta_rho') ;
  nc{'zeta_west'}.long_name = ncchar('western boundary sea surface height');
  nc{'zeta_west'}.long_name = 'western boundary sea surface height';
  nc{'zeta_west'}.units = ncchar('meter');
  nc{'zeta_west'}.units = 'meter';
  nc{'zeta_west'}.coordinates = ncchar('lat_rho ocean_time');
  nc{'zeta_west'}.coordinates = 'lat_rho ocean_time';
  nc{'zeta_west'}.time = ncchar('ocean_time');
  nc{'zeta_west'}.time = 'ocean_time';
%
end
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(bryname);
nc.clim_file = bryname;
nc.grd_file = ncchar(grdname);
nc.grd_file = grdname;
nc.type = ncchar(type);
nc.type = type;
nc.history = ncchar(history);
nc.history = history;
%
% Leave define mode
%
%%result = endef(nc);
%
% Compute S coordinates
%
[sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%disp(['vtransform=',num2str(vtransform)])
%
% Write variables
%
nc{'spherical'}(:)='T';
nc{'Vtransform'}(:)=vtransform;
nc{'Vstretching'}(:)=4;
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'sc_r'}(:) = sc_r;
nc{'sc_w'}(:) = sc_w;
nc{'Cs_r'}(:) = Cs_r ; 
nc{'Cs_w'}(:) = Cs_w;
nc{'ocean_time'}(:) =  time;  
if obc(1)==1
  nc{'u_south'}(:) =  0; 
  nc{'v_south'}(:) =  0; 
  nc{'ubar_south'}(:) =  0; 
  nc{'vbar_south'}(:) =  0; 
  nc{'zeta_south'}(:) =  0; 
  nc{'temp_south'}(:) =  0; 
  nc{'salt_south'}(:) =  0;
end 
if obc(2)==1
  nc{'u_east'}(:) =  0; 
  nc{'v_east'}(:) =  0; 
  nc{'ubar_east'}(:) =  0; 
  nc{'vbar_east'}(:) =  0; 
  nc{'zeta_east'}(:) =  0; 
  nc{'temp_east'}(:) =  0; 
  nc{'salt_east'}(:) =  0;
  nc{'mud_east_01'}(:)  = 0.001;
  nc{'sand_east_01'}(:) = 0.001;
  nc{'sand_east_02'}(:) = 0.0005;
end 
if obc(3)==1
  nc{'u_north'}(:) =  0; 
  nc{'v_north'}(:) =  0; 
  nc{'ubar_north'}(:) =  0; 
  nc{'vbar_north'}(:) =  0; 
  nc{'zeta_north'}(:) =  0; 
  nc{'temp_north'}(:) =  0; 
  nc{'salt_north'}(:) =  0;
end 
if obc(4)==1
  nc{'u_west'}(:) =  0; 
  nc{'v_west'}(:) =  0; 
  nc{'ubar_west'}(:) =  0; 
  nc{'vbar_west'}(:) =  0; 
  nc{'zeta_west'}(:) =  0; 
  nc{'temp_west'}(:) =  0; 
  nc{'salt_west'}(:) =  0;
end 
close(nc)
return


