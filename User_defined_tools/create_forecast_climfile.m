function create_forecast_climfile(clmname,grdname,title,...
                         theta_s,theta_b,hc,N,...
                         time,timeunits,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   cycle        Length (days) for cycling the climatology.(Real)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',clmname])
disp(' ')
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file
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
%  Create the climatology file
%
type = 'CLIMATOLOGY file' ; 
history = 'ROMS Climfile' ;
ncid=netcdf.create(clmname,'64BIT_OFFSET');
netcdf.close(ncid);
nc = netcdf(clmname,'write');
% result = redef(nc);
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
%  Create variables
%
nc{'spherical'} = ncchar('one') ;
nc{'Vtransform'} = ncint('one') ;
nc{'Vstretching'} = ncint('one') ;
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_b'} = ncdouble('one') ;
nc{'Tcline'} = ncdouble('one') ;
nc{'hc'} = ncdouble('one') ;
nc{'sc_r'} = ncdouble('s_rho') ;
nc{'sc_w'} = ncdouble('s_w') ;
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'Cs_w'} = ncdouble('s_w') ;
nc{'ocean_time'} = ncdouble('ocean_time') ;
nc{'temp'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'salt'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'u'} = ncdouble('ocean_time','s_rho','eta_u','xi_u') ;
nc{'v'} = ncdouble('ocean_time','s_rho','eta_v','xi_v') ;
nc{'ubar'} = ncdouble('ocean_time','eta_u','xi_u') ;
nc{'vbar'} = ncdouble('ocean_time','eta_v','xi_v') ;
nc{'SSH'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'zeta'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
%
%  Create attributes
%
nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
nc{'spherical'}.long_name = ncchar('grid type logical switch');
nc{'spherical'}.long_name = 'grid type logical switch';
nc{'spherical'}.flag_values = ncchar('T, F');
nc{'spherical'}.flag_values = 'T, F';
nc{'spherical'}.flag_meanings = ncchar('spherical Cartesian');
nc{'spherical'}.flag_meanings = 'spherical Cartesian';
%
nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = ncchar('nondimensional');
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = ncchar('nondimensional');
nc{'theta_b'}.units = 'nondimensional';
%
nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nc{'Tcline'}.units = ncchar('meter');
nc{'Tcline'}.units = 'meter';
%
nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = ncchar('meter');
nc{'hc'}.units = 'meter';
%
nc{'sc_r'}.long_name = ncchar('S-coordinate at RHO-points');
nc{'sc_r'}.long_name = 'S-coordinate at RHO-points';
nc{'sc_r'}.valid_min = -1.;
nc{'sc_r'}.valid_max = 0.;
nc{'sc_r'}.positive = ncchar('up');
nc{'sc_r'}.positive = 'up';
if (vtransform ==1)
    nc{'sc_r'}.standard_name = ncchar('ocena_s_coordinate_g1');
    nc{'sc_r'}.standard_name = 'ocena_s_coordinate_g1';
elseif (vtransform ==2)
    nc{'sc_r'}.standard_name = ncchar('ocena_s_coordinate_g2');
    nc{'sc_r'}.standard_name = 'ocena_s_coordinate_g2';     
end
nc{'sc_r'}.formula_terms = ncchar('s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
nc{'sc_r'}.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
%
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
nc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_r'}.units = ncchar('nondimensional');
nc{'Cs_r'}.units = 'nondimensional';
nc{'Cs_r'}.valid_min = -1;
nc{'Cs_r'}.valid_max = 0;
%
nc{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc{'Cs_w'}.long_name = 'S-coordinate stretching curves at W-points';
nc{'Cs_w'}.units = ncchar('nondimensional');
nc{'Cs_w'}.units = 'nondimensional';
nc{'Cs_w'}.valid_min = -1;
nc{'Cs_w'}.valid_max = 0;
%
nc{'ocean_time'}.long_name = ncchar('time for climatology');
nc{'ocean_time'}.long_name = 'time for climatology';
nc{'ocean_time'}.units = ncchar(timeunits);
nc{'ocean_time'}.units = timeunits;
nc{'ocean_time'}.calendar = ncchar('gregorian');
nc{'ocean_time'}.calendar = 'gregorian';
%
nc{'temp'}.long_name = ncchar('potential temperature');
nc{'temp'}.long_name = 'potential temperature';
nc{'temp'}.units = ncchar('Celsius');
nc{'temp'}.units = 'Celsius';
nc{'temp'}.time = ncchar('temp_time');
nc{'temp'}.time = 'temp_time';
nc{'temp'}.coordinates = ncchar('lon_rho lat_rho s_rho temp_time');
nc{'temp'}.coordinates = 'lon_rho lat_rho s_rho temp_time';
%
nc{'salt'}.long_name = ncchar('salinity');
nc{'salt'}.long_name = 'salinity';
nc{'salt'}.units = ncchar('PSU');
nc{'salt'}.units = 'PSU';
nc{'salt'}.time = ncchar('salt_time');
nc{'salt'}.time = 'salt_time';
nc{'salt'}.coordinates = ncchar('lon_rho lat_rho s_rho salt_time');
nc{'salt'}.coordinates = 'lon_rho lat_rho s_rho salt_time';
%
nc{'u'}.long_name = ncchar('u-momentum component');
nc{'u'}.long_name = 'u-momentum component';
nc{'u'}.units = ncchar('meter second-1');
nc{'u'}.units = 'meter second-1';
nc{'u'}.time = ncchar('uclm_time');
nc{'u'}.time = 'uclm_time';
nc{'u'}.coordinates = ncchar('lon_u lat_u s_rho uclm_time');
nc{'u'}.coordinates = 'lon_u lat_u s_rho u_time';
%
nc{'v'}.long_name = ncchar('v-momentum component');
nc{'v'}.long_name = 'v-momentum component';
nc{'v'}.units = ncchar('meter second-1');
nc{'v'}.units = 'meter second-1';
nc{'v'}.time = ncchar('vclm_time');
nc{'v'}.time = 'vclm_time';
nc{'v'}.coordinates = ncchar('lon_v lat_v s_rho vclm_time');
nc{'v'}.coordinates = 'lon_v lat_v s_rho vclm_time';
%
nc{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
nc{'ubar'}.long_name = 'vertically integrated u-momentum component';
nc{'ubar'}.units = ncchar('meter second-1');
nc{'ubar'}.units = 'meter second-1';
nc{'ubar'}.time = ncchar('uclm_time');
nc{'ubar'}.time = 'uclm_time';
nc{'ubar'}.coordinates = ncchar('lon_u lat_u uclm_time');
nc{'ubar'}.coordinates = 'lon_u lat_u uclm_time';
%
nc{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
nc{'vbar'}.long_name = 'vertically integrated v-momentum component';
nc{'vbar'}.units = ncchar('meter second-1');
nc{'vbar'}.units = 'meter second-1';
nc{'vbar'}.time = ncchar('vclm_time');
nc{'vbar'}.time = 'vclm_time';
nc{'vbar'}.coordinates = ncchar('lon_v lat_v vclm_time');
nc{'vbar'}.coordinates = 'lon_v lat_v vclm_time';
%
nc{'SSH'}.long_name = ncchar('sea surface height');
nc{'SSH'}.long_name = 'sea surface height';
nc{'SSH'}.units = ncchar('meter');
nc{'SSH'}.units = 'meter';
nc{'SSH'}.time = ncchar('zeta_time');
nc{'SSH'}.time = 'zeta_time';
nc{'SSH'}.coordinates = ncchar('lon_rho lat_rho zeta_time');
nc{'SSH'}.coordinates = 'lon_rho lat_rho zeta_time';
%
nc{'zeta'}.long_name = ncchar('sea surface height');
nc{'zeta'}.long_name = 'sea surface height';
nc{'zeta'}.units = ncchar('meter');
nc{'zeta'}.units = 'meter';
nc{'zeta'}.time = ncchar('zeta_time');
nc{'zeta'}.time = 'zeta_time';
nc{'zeta'}.coordinates = ncchar('lon_rho lat_rho zeta_time');
nc{'zeta'}.coordinates = 'lon_rho lat_rho zeta_time';
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(clmname);
nc.clim_file = clmname;
nc.grd_file = ncchar(grdname);
nc.grd_file = grdname;
nc.type = ncchar(type);
nc.type = type;
nc.history = ncchar(history);
nc.history = history;
%
% Leave define mode
%
%result = endef(nc);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
[sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%disp(['vtransform=',num2str(vtransform)])

% cff1=1./sinh(theta_s);
% cff2=0.5/tanh(0.5*theta_s);
% sc_r=((1:N)-N-0.5)/N;
% Cs_r=(1.-theta_b)*cff1*sinh(theta_s*sc_r)...
%     +theta_b*(cff2*tanh(theta_s*(sc_r+0.5))-0.5);
% sc_w=((0:N)-N)/N;
% Cs_w=(1.-theta_b)*cff1*sinh(theta_s*sc_w)...
%     +theta_b*(cff2*tanh(theta_s*(sc_w+0.5))-0.5);


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
nc{'Cs_r'}(:) =  Cs_r; 
nc{'Cs_w'}(:) = Cs_w;
nc{'ocean_time'}(:) =  time; 
nc{'u'}(:) =  0; 
nc{'v'}(:) =  0; 
nc{'ubar'}(:) =  0; 
nc{'vbar'}(:) =  0; 
nc{'SSH'}(:) =  0; 
nc{'zeta'}(:) =  0; 
nc{'temp'}(:) =  0; 
nc{'salt'}(:) =  0; 
close(nc)
return