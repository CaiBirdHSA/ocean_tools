function create_sediment_inifile(inifile,gridfile,title,...
                         theta_s,theta_b,hc,N,time,clobber,vtransform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',inifile])
if nargin < 10
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
NB=5;     %layer number of bed 
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'ROMS Sediment' ;
ncid=netcdf.create(inifile,'64BIT_OFFSET');
netcdf.close(ncid);
nc = netcdf(inifile,clobber);
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
nc('s_bed') = NB;
nc('one') = 1;
%
%  Create variables
%
nc{'spherical'} = ncchar('one') ;
nc{'Vtransform'} = ncint('one') ;
nc{'Vstretching'} = ncint('one') ;
nc{'tstart'} = ncdouble('one') ;
nc{'tend'} = ncdouble('one') ;
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_b'} = ncdouble('one') ;
nc{'Tcline'} = ncdouble('one') ;
nc{'hc'} = ncdouble('one') ;
nc{'sc_r'} = ncdouble('s_rho') ;
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'ocean_time'} = ncdouble('ocean_time') ;
nc{'scrum_time'} = ncdouble('ocean_time') ;
nc{'u'} = ncdouble('ocean_time','s_rho','eta_u','xi_u') ;
nc{'v'} = ncdouble('ocean_time','s_rho','eta_v','xi_v') ;
nc{'ubar'} = ncdouble('ocean_time','eta_u','xi_u') ;
nc{'vbar'} = ncdouble('ocean_time','eta_v','xi_v') ;
nc{'zeta'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'temp'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'salt'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'mud_01'}  = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'sand_01'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'sand_02'} = ncdouble('ocean_time','s_rho','eta_rho','xi_rho') ;
nc{'mudfrac_01'}  = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'sandfrac_01'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'sandfrac_02'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'mudmass_01'}  = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'sandmass_01'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'sandmass_02'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'bedload_Umud_01'}  = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bedload_Usand_01'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bedload_Usand_02'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bedload_Vmud_01'}  = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bedload_Vsand_01'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bedload_Vsand_02'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bed_thickness'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'bed_age'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'bed_porosity'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'bed_biodiff'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'bed_tau_crit'} = ncdouble('ocean_time','s_bed','eta_rho','xi_rho') ;
nc{'grain_diameter'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'grain_density'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'settling_vel'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'erosion_stress'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'ripple_length'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'ripple_height'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'bed_wave_amp'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'rdrag'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'rdrag2'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'ZoBot'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_def'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_app'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_Nik'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_bio'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_bedform'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_bedload'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'Zo_wbl'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
nc{'active_layer_thickness'} = ncdouble('ocean_time','eta_rho','xi_rho') ;
%
%  Create attributes
%
nc{'Vtransform'}.long_name = ncchar('vertical terrain-following transformation equation');
nc{'Vtransform'}.long_name = 'vertical terrain-following transformation equation';
%
nc{'Vstretching'}.long_name = ncchar('vertical terrain-following stretching function');
nc{'Vstretching'}.long_name = 'vertical terrain-following stretching function';
%
nc{'tstart'}.long_name = ncchar('start processing day');
nc{'tstart'}.long_name = 'start processing day';
nc{'tstart'}.units = ncchar('day');
nc{'tstart'}.units = 'day';
%
nc{'tend'}.long_name = ncchar('end processing day');
nc{'tend'}.long_name = 'end processing day';
nc{'tend'}.units = ncchar('day');
nc{'tend'}.units = 'day';
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
nc{'ocean_time'}.long_name = ncchar('time since 1980-01-01 00:00:00');
nc{'ocean_time'}.long_name = 'time since 1980-01-01 00:00:00';
nc{'ocean_time'}.units = ncchar('seconds since 1980-01-01 00:00:00');
nc{'ocean_time'}.units = 'seconds since 1980-01-01 00:00:00';
%
nc{'u'}.long_name = ncchar('u-momentum component');
nc{'u'}.long_name = 'u-momentum component';
nc{'u'}.units = ncchar('meter second-1');
nc{'u'}.units = 'meter second-1';
%
nc{'v'}.long_name = ncchar('v-momentum component');
nc{'v'}.long_name = 'v-momentum component';
nc{'v'}.units = ncchar('meter second-1');
nc{'v'}.units = 'meter second-1';
%
nc{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
nc{'ubar'}.long_name = 'vertically integrated u-momentum component';
nc{'ubar'}.units = ncchar('meter second-1');
nc{'ubar'}.units = 'meter second-1';
%
nc{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
nc{'vbar'}.long_name = 'vertically integrated v-momentum component';
nc{'vbar'}.units = ncchar('meter second-1');
nc{'vbar'}.units = 'meter second-1';
%
nc{'zeta'}.long_name = ncchar('free-surface');
nc{'zeta'}.long_name = 'free-surface';
nc{'zeta'}.units = ncchar('meter');
nc{'zeta'}.units = 'meter';
%
nc{'temp'}.long_name = ncchar('potential temperature');
nc{'temp'}.long_name = 'potential temperature';
nc{'temp'}.units = ncchar('Celsius');
nc{'temp'}.units = 'Celsius';
%
nc{'salt'}.long_name = ncchar('salinity');
nc{'salt'}.long_name = 'salinity';
nc{'salt'}.units = ncchar('PSU');
nc{'salt'}.units = 'PSU';
%
nc{'mud_01'}.long_name = ncchar('suspended cohesive sediment');
nc{'mud_01'}.long_name = 'suspended cohesive sediment';
nc{'mud_01'}.units = ncchar('kilogram meter-3');
nc{'mud_01'}.units = 'kilogram meter-3';
%
nc{'sand_01'}.long_name = ncchar('suspended noncohesive sediment 1');
nc{'sand_01'}.long_name = 'suspended noncohesive sediment 1';
nc{'sand_01'}.units = ncchar('kilogram meter-3');
nc{'sand_01'}.units = 'kilogram meter-3';
%
nc{'sand_02'}.long_name = ncchar('suspended noncohesive sediment 2');
nc{'sand_02'}.long_name = 'suspended noncohesive sediment 2';
nc{'sand_02'}.units = ncchar('kilogram meter-3');
nc{'sand_02'}.units = 'kilogram meter-3';
%
nc{'mudfrac_01'}.long_name = ncchar('suspended noncohesive sediment 2');
nc{'mudfrac_01'}.long_name = 'suspended noncohesive sediment 2';
nc{'mudfrac_01'}.units = ncchar('kilogram meter-3');
nc{'mudfrac_01'}.units = 'kilogram meter-3';
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
nc{'tstart'}(:) =  time; 
nc{'tend'}(:) =  time; 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'sc_r'}(:) =  sc_r; 
nc{'Cs_r'}(:) =  Cs_r; 
nc{'scrum_time'}(1) =  time*24*3600; 
nc{'ocean_time'}(1) =  time*24*3600; 
nc{'u'}(:) =  0; 
nc{'v'}(:) =  0; 
nc{'zeta'}(:) =  0; 
nc{'ubar'}(:) =  0; 
nc{'vbar'}(:) =  0; 
nc{'temp'}(:) =  0; 
nc{'salt'}(:) =  0; 
nc{'mud_01'}(:) =  0; 
nc{'sand_01'}(:) =  0; 
nc{'sand_02'}(:) =  0; 
nc{'mudfrac_01'}(:) =  0; 
nc{'sandfrac_01'}(:) =  0; 
nc{'sandfrac_02'}(:) =  0; 
nc{'mudmass_01'}(:) =  0; 
nc{'sandmass_01'}(:) =  0; 
nc{'sandmass_02'}(:) =  0; 
nc{'bedload_Umud_01'}(:) =  0; 
nc{'bedload_Usand_01'}(:) =  0; 
nc{'bedload_Usand_02'}(:) =  0; 
nc{'bedload_Vmud_01'}(:) =  0; 
nc{'bedload_Vsand_01'}(:) =  0; 
nc{'bedload_Vsand_02'}(:) =  0; 
nc{'bed_thickness'}(:) =  0.1; 
nc{'bed_age'}(:) =  864000; 
nc{'bed_porosity'}(:) =  0.02; 
nc{'bed_biodiff'}(:) =  0.01; 
nc{'bed_tau_crit'}(:) =  0.04; 
nc{'grain_diameter'}(:) =  0.000032; 
nc{'grain_density'}(:) =  1032; 
nc{'settling_vel'}(:) =  0.0002; 
nc{'erosion_stress'}(:) =  0.002; 
nc{'ripple_length'}(:) =  0.5; 
nc{'ripple_height'}(:) =  0.05; 
nc{'bed_wave_amp'}(:) =  0.01; 
nc{'rdrag'}(:) =  0.001; 
nc{'rdrag2'}(:) = 0.001; 
nc{'ZoBot'}(:) = 0.0001; 
nc{'Zo_def'}(:) = 0.0001; 
nc{'Zo_app'}(:) = 0.0001; 
nc{'Zo_Nik'}(:) = 0.001; 
nc{'Zo_bio'}(:) = 0.001; 
nc{'Zo_bedform'}(:) = 0.001; 
nc{'Zo_bedload'}(:) = 0.001; 
nc{'Zo_wbl'}(:) = 0.001; 
nc{'active_layer_thickness'}(:) = 0.02; 
%
% Synchronize on disk
%
close(nc);
return


