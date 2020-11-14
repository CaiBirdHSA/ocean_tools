function nc_add_tides(tidname,grdname,Ntides,start_tide_mjd,components)
%
nc = netcdf(grdname, 'nowrite');
maskr=nc{'mask_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
close(nc);
%
type = 'TIDE file' ; 
ncid=netcdf.create(tidname,'64BIT_OFFSET');
netcdf.close(ncid);
nc = netcdf(tidname,'clobber');
%
%%redef(nc);		% for Octave compatibility
%
%  Add dimension
%
nc('xi_rho') = Lp;
nc('eta_rho') = Mp;
nc('tide_period')=Ntides;
%
%  Add variables and attributes
%
nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho');
nc{'mask_rho'}.long_name = ncchar('mask on RHO-points');
nc{'mask_rho'}.long_name = 'mask on RHO-points';
nc{'mask_rho'}.option_0 = ncchar('land');
nc{'mask_rho'}.option_0 = 'land';
nc{'mask_rho'}.option_1 = ncchar('water');
nc{'mask_rho'}.option_1 = 'water';

nc{'tide_period'} = ncdouble('tide_period');
nc{'tide_period'}.long_name = ncchar('Tide angular period');
nc{'tide_period'}.long_name = 'Tide angular period';
nc{'tide_period'}.units = ncchar('Hours');
nc{'tide_period'}.units = 'Hours';

nc{'tide_Ephase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Ephase'}.long_name = ncchar('Tidal elevation phase angle');
nc{'tide_Ephase'}.long_name = 'Tidal elevation phase angle';
nc{'tide_Ephase'}.units = ncchar('Degrees');
nc{'tide_Ephase'}.units = 'Degrees';

nc{'tide_Eamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Eamp'}.long_name = ncchar('Tidal elevation amplitude');
nc{'tide_Eamp'}.long_name = 'Tidal elevation amplitude';
nc{'tide_Eamp'}.units = ncchar('Meter');
nc{'tide_Eamp'}.units = 'Meter';

nc{'tide_Cmin'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cmin'}.long_name = ncchar('Tidal current ellipse semi-minor axis');
nc{'tide_Cmin'}.long_name = 'Tidal current ellipse semi-minor axis';
nc{'tide_Cmin'}.units = ncchar('Meter second-1');
nc{'tide_Cmin'}.units = 'Meter second-1';

nc{'tide_Cmax'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cmax'}.long_name = ncchar('Tidal current, ellipse semi-major axis');
nc{'tide_Cmax'}.long_name = 'Tidal current, ellipse semi-major axis';
nc{'tide_Cmax'}.units = ncchar('Meter second-1');
nc{'tide_Cmax'}.units = 'Meter second-1';

nc{'tide_Cangle'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cangle'}.long_name = ncchar('Tidal current inclination angle');
nc{'tide_Cangle'}.long_name = 'Tidal current inclination angle';
nc{'tide_Cangle'}.units = ncchar('Degrees between semi-major axis and East');
nc{'tide_Cangle'}.units = 'Degrees between semi-major axis and East';

nc{'tide_Cphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cphase'}.long_name = ncchar('Tidal current phase angle');
nc{'tide_Cphase'}.long_name = 'Tidal current phase angle';
nc{'tide_Cphase'}.units = ncchar('Degrees');
nc{'tide_Cphase'}.units = 'Degrees';

nc{'tide_Pamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Pamp'}.long_name = ncchar('Tidal potential amplitude');
nc{'tide_Pamp'}.long_name = 'Tidal potential amplitude';
nc{'tide_Pamp'}.units = ncchar('Meter');
nc{'tide_Pamp'}.units = 'Meter';

nc{'tide_Pphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Pphase'}.long_name = ncchar('Tidal potential phase angle');
nc{'tide_Pphase'}.long_name = 'Tidal potential phase angle';
nc{'tide_Pphase'}.units = ncchar('Degrees');
nc{'tide_Pphase'}.units = 'Degrees';

nc{'mask_rho'}(:) =  maskr;
nc.date = ncchar(date);
nc.date = date;
nc.type = ncchar(type);
nc.type = type;
nc.start_tide_date=mjd2greg(start_tide_mjd);
nc.start_tide_mjd=start_tide_mjd;
nc.components = ncchar(components);
nc.components = components;

close(nc)
