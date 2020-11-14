function create_forecast_inifile(inifile,gridfile,title,...
                         theta_s,theta_b,hc,N,time,timeunits,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',inifile])
if nargin < 9
   disp([' NO VTRANSFORM parameter found'])
   disp([' USE TRANSFORM default value vtransform = 1'])
   vtransform = 1; 
end
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file
%
nc=netcdf(gridfile,'r');
h=nc{'h'}(:);  
mask=nc{'mask_rho'}(:);
close(nc);
hmin=min(min(h(mask==1)));
if vtransform ==1
    if hc > hmin
        error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
    end
end
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'HYCOM' ;
ncid=netcdf.create(inifile,'64BIT_OFFSET');
netcdf.close(ncid);
nc = netcdf(inifile,'write');
%result = redef(nc);
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
nc('ocean_time') = 1;
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
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'ocean_time'} = ncdouble('ocean_time') ;
nc{'u'} = ncdouble('ocean_time','s_rho','eta_u','xi_u') ;
nc{'v'} = ncdouble('ocean_time','s_rho','eta_v','xi_v') ;
nc{'ubar'} = ncdouble('ocean_time','eta_u','xi_u') ;
nc{'vbar'} = ncdouble('ocean_time','eta_v','xi_v') ;
nc{'zeta'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'temp'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'salt'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
%
%  Create attributes
%
nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
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
nc{'sc_r'}.units = ncchar('nondimensional');
nc{'sc_r'}.units = 'nondimensional';
nc{'sc_r'}.valid_min = -1;
nc{'sc_r'}.valid_max = 0;
%
nc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_r'}.units = ncchar('nondimensional');
nc{'Cs_r'}.units = 'nondimensional';
nc{'Cs_r'}.valid_min = -1;
nc{'Cs_r'}.valid_max = 0;
%
nc{'ocean_time'}.long_name = ncchar('time since initialization');
nc{'ocean_time'}.long_name = 'time since initialization';
nc{'ocean_time'}.units = ncchar(timeunits);
nc{'ocean_time'}.units = timeunits;
nc{'ocean_time'}.calendar = ncchar('gregorian');
nc{'ocean_time'}.calendar = 'gregorian';
%
nc{'u'}.long_name = ncchar('u-momentum component');
nc{'u'}.long_name = 'u-momentum component';
nc{'u'}.units = ncchar('meter second-1');
nc{'u'}.units = 'meter second-1';
nc{'u'}.time = ncchar('ocean_time');
nc{'u'}.time = 'ocean_time';
%
nc{'v'}.long_name = ncchar('v-momentum component');
nc{'v'}.long_name = 'v-momentum component';
nc{'v'}.units = ncchar('meter second-1');
nc{'v'}.units = 'meter second-1';
nc{'v'}.time = ncchar('ocean_time');
nc{'v'}.time = 'ocean_time';
%
nc{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
nc{'ubar'}.long_name = 'vertically integrated u-momentum component';
nc{'ubar'}.units = ncchar('meter second-1');
nc{'ubar'}.units = 'meter second-1';
nc{'ubar'}.time = ncchar('ocean_time');
nc{'ubar'}.time = 'ocean_time';
%
nc{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
nc{'vbar'}.long_name = 'vertically integrated v-momentum component';
nc{'vbar'}.units = ncchar('meter second-1');
nc{'vbar'}.units = 'meter second-1';
nc{'vbar'}.time = ncchar('ocean_time');
nc{'vbar'}.time = 'ocean_time';
%
nc{'zeta'}.long_name = ncchar('free-surface');
nc{'zeta'}.long_name = 'free-surface';
nc{'zeta'}.units = ncchar('meter');
nc{'zeta'}.units = 'meter';
nc{'zeta'}.time = ncchar('ocean_time');
nc{'zeta'}.time = 'ocean_time';
%
nc{'temp'}.long_name = ncchar('potential temperature');
nc{'temp'}.long_name = 'potential temperature';
nc{'temp'}.units = ncchar('Celsius');
nc{'temp'}.units = 'Celsius';
nc{'temp'}.time = ncchar('ocean_time');
nc{'temp'}.time = 'ocean_time';
%
nc{'salt'}.long_name = ncchar('salinity');
nc{'salt'}.long_name = 'salinity';
nc{'salt'}.units = ncchar('PSU');
nc{'salt'}.units = 'PSU';
nc{'salt'}.time = ncchar('ocean_time');
nc{'salt'}.time = 'ocean_time';
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(inifile);
nc.clim_file = inifile;
nc.grd_file = ncchar(gridfile);
nc.grd_file = gridfile;
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
nc{'sc_r'}(:) =  sc_r; 
nc{'Cs_r'}(:) =  Cs_r; 
nc{'ocean_time'}(:) = time*1.000000000000000d0; 
nc{'u'}(:) =  0; 
nc{'v'}(:) =  0; 
nc{'zeta'}(:) =  0; 
nc{'ubar'}(:) =  0; 
nc{'vbar'}(:) =  0; 
nc{'temp'}(:) =  0; 
nc{'salt'}(:) =  0; 
%
% Synchronize on disk
%
close(nc);
return