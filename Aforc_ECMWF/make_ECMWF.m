%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  make_ECMWF_daily.m
% 
%  Create and fill frc and bulk files with ECMWF data.
%  (ERA INTERIM Reanalysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param
%
  frc_prefix=[frc_prefix,'_ECMWF_'];                           
  blk_prefix=[blk_prefix,'_ECMWF_'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if level==0
  nc_suffix='.nc';
else
  nc_suffix=['.nc.',num2str(level)];
  grdname=[grdname,'.',num2str(level)];
end
%
% Get the model grid
%
disp(' ')
disp([' Read the grid in ',grdname])
nc=netcdf(grdname);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
lonu=nc{'lon_u'}(:);
latu=nc{'lat_u'}(:);
lonv=nc{'lon_v'}(:);
latv=nc{'lat_v'}(:);
angle=nc{'angle'}(:);
close(nc);

if makefrc==1 || makeblk==1
  %
  % Get the ECMWF horizontal grids (it should be the same for every month)
  %
  nc=netcdf([ECMWF_dir,'land_Y',num2str(Ymin),'M',num2str(Mmin),'.nc']);
  disp(['Use this land file :',char([ECMWF_dir,'land_Y',num2str(Ymin),'M',num2str(Mmin),'.nc'])])

  lon1=nc{'lon'}(:);
  lat1=nc{'lat'}(:);
  [lon1,lat1]=meshgrid(lon1,lat1);
  
  
  mask=1-squeeze(nc{'land'}(:));
  mask(mask==0)=NaN;  
  close(nc);
  
  %
  %Loop on the years and the months
  %
  disp(['====================='])
  disp(['INTERPOLATION STEP'])
  disp(['====================='])
  disp(['Loop on the years and the months'])

  %
   for Y=Ymin:Ymax
    if Y==Ymin 
      mo_min=Mmin;
    else
      mo_min=1;
    end
    if Y==Ymax
      mo_max=Mmax;
    else
      mo_max=12;
    end
     for M=mo_min:mo_max
      disp(' ')
      disp(['Processing  year ',num2str(Y),...
            ' - month ',num2str(M)])
      disp(' ')
      %-------------------------------------------------------------------%
      %
      % Process time (here in days), with SST file (minimum time step)
      %
      %-------------------------------------------------------------------%
      nc=netcdf([ECMWF_dir,'SSTK_Y',num2str(Y),'M',num2str(M),'.nc']);
      ECMWF_time=nc{'time'}(:);
      close(nc);
      dt=mean(gradient(ECMWF_time));
      disp(['dt=',num2str(dt)])
      %-----------------------------------------------------------
      %Variable overlapping timesteps : 2 at the beginning and 2 at the end
      %------------------------------------------------------------
      tlen0=length(ECMWF_time);
      disp(['tlen0=',num2str(tlen0)])
      freq=4; % 6-hourly  (freq=1 for daily)
      itolap=freq*itolap_ecmwf;
      tlen=tlen0+2*itolap;
      disp(['tlen=',num2str(tlen)])
      disp(['Overlap is ',num2str(itolap_ecmwf),' days'])
      disp(['Overlap is ',num2str(itolap_ecmwf),' days before and after'])     
      time=0*(1:tlen);
      time(itolap+1:tlen0+itolap)=ECMWF_time;   
      disp(['====================='])
      disp('Compute time for croco file')
      disp(['====================='])
      for aa=1:itolap
        time(aa)=time(itolap+1)-(itolap+1-aa)*dt;
      end
      for aa=1:itolap
        time(tlen0+itolap+aa)=time(tlen0+itolap)+aa*dt;
      end
      %-------------------------------------------------------------------%
      %
      % Create the CROCO forcing files
      %
      % ------------------------------------------------------------------%
      %
      disp(['====================='])
      disp('Create the frc/blk netcdf file')
      disp(['====================='])
      %
      blkname=[blk_prefix,'Y',num2str(Y),...
               'M',num2str(M),nc_suffix];
      frcname=[frc_prefix,'Y',num2str(Y),...
               'M',num2str(M),nc_suffix];
      if makeblk==1
        disp(['Create a new bulk file: ' blkname])
        create_bulk(blkname,grdname,CROCO_title,time,0);
	    disp([' '])
      end
      if makefrc==1
        disp(['Create a new forcing file: ' frcname])
	disp([' '])
        create_forcing(frcname,grdname,CROCO_title,...
                       time,0,0,...
                       0,0,0,...
                       0,0,0,0,0,0)
      end
      %
      % Add the tides (needs to be tested for this version of make_ECMWF)
      %
      if add_tides==1
        add_tidal_data(tidename,grdname,frcname,Ntides,tidalrank,...
                                            Yorig,Y,M,coastfileplot)
      end
      %
      % Open the CROCO forcing files
      if makefrc==1
        nc_frc=netcdf(frcname,'write');
      else
        nc_frc=[];
      end
      if makeblk==1
        nc_blk=netcdf(blkname,'write');
      else
        nc_blk=[];
      end
      %
      % Check if there are ECMWF files for the previous Month
      Mm=M-1;
      Ym=Y;
      if Mm==0
        Mm=12;
        Ym=Y-1;
      end
      %
      fname = [ECMWF_dir,'TP_Y',num2str(Ym),'M',num2str(Mm),'.nc'];
      %nc=netcdf([ECMWF_dir,'TP_Y',num2str(Ym),'M',num2str(Mm),'.nc']);
      %
      disp(' ')
      disp('======================================================')
      disp('Perform interpolations for the previous month')      
      disp('======================================================')
      disp(' ')
      if exist(fname)==0
        disp(['No data for the previous month: using current month'])
        tndx=1;
        Mm=M;
        Ym=Y;
      else
        nc=netcdf(fname);
        tndx=length(nc('time'));
	    if makefrc==1
          for aa=1:itolap
	     nc_frc{'sms_time'}(aa)=nc{'time'}(tndx-(itolap-aa));
          end
        end
        %
        if makeblk==1
          for aa=1:itolap
	      nc_blk{'bulk_time'}(aa)=nc{'time'}(tndx-(itolap-aa));
	  end
        end
        close(nc)
      end
      %
      % Perform interpolations for the previous month or repeat the first one
      %
      for aa=1:itolap
        if exist(fname)==0
	  aa0=aa;
        else
          aa0=tndx-(itolap-aa);
        end
        interp_ECMWF(ECMWF_dir,Ym,Mm,Roa,interp_method,lon1,lat1,...
                     mask,aa0,nc_frc,nc_blk,lon,lat,angle,aa)
      end  
      %######################################################################      
      %   
      disp(' ')
      disp('======================================================')
      disp('Perform interpolations for the current month')
      disp('======================================================')
      disp(' ')

      % Perform interpolations for the current month
      %

      for tndx=1:tlen0
        if mod(tndx,5)==0
          disp(['Step: ',num2str(tndx),' of ',num2str(tlen0)])
        end
	    interp_ECMWF(ECMWF_dir,Y,M,Roa,interp_method,lon1,lat1,...
		    mask,tndx,nc_frc,nc_blk,lon,lat,angle,tndx+itolap)	
      end

      disp(' ')      
      disp('======================================================')
      disp('Perform interpolations for next month')    
      disp('======================================================')
      disp(' ')
      %######################################################################
      % Read ECMWF file for the next month
      %
      Mp=M+1;
      Yp=Y;
      if Mp==13
        Mp=1;
        Yp=Y+1;
      end
      
      fname=[ECMWF_dir,'TP_Y',num2str(Yp),'M',num2str(Mp),'.nc'];
      
      if exist(fname)==0
        disp(['No data for the next month: using current month'])
        tndx=tlen0;
        Mp=M;
        Yp=Y;
      else
	nc=netcdf(fname);
	if makefrc==1
	   disp('sms_time')
	   for tndx=tlen0+itolap+1:tlen;
	      nc_frc{'sms_time'}(tndx)=nc{'time'}(tndx-tlen0-(itolap));
	   end;
	end
	%
        if makeblk==1
	   disp('bulk_time')
           for tndx=tlen0+itolap+1:tlen;
              nc_blk{'bulk_time'}(tndx)=nc{'time'}(tndx-tlen0-(itolap));
           end;
        end
        close(nc)
      end
      %
      % Perform the interpolations for the next month
      %
      disp('Last steps')
      for tndx=tlen0+itolap+1:tlen;
        disp(['tndx= ',num2str(tndx)])
        tout=tndx;
    	disp(['tout=tndx ',num2str(tndx)])
        if Mp==M
          %tin=tlen0;       % persistency if current month is used
          tin=tndx-2*itolap
          disp(['tin=',num2str(tin)])
        else
          tin=tndx-tlen0-itolap;
          disp(['tin=',num2str(tin)])
        end
        interp_ECMWF(ECMWF_dir,Yp,Mp,Roa,interp_method,lon1,lat1,...
		     mask,tin,nc_frc,nc_blk,lon,lat,angle,tout)           
      end;
      %
      % Close the CROCO forcing files
      %
      if ~isempty(nc_frc)
        close(nc_frc);
      end
      if ~isempty(nc_blk)
        close(nc_blk);
      end
    end
   end
 end
%
% Spin-up: (reproduce the first year 'SPIN_Long' times)
% just copy the files for the first year and change the time
%
disp('======================================================')
disp('Add spin up phase')      
if SPIN_Long>0
  M=Mmin-1;
  Y=Ymin-SPIN_Long;
  for month=1:12*SPIN_Long
    M=M+1;
    if M==13
      M=1; 
      Y=Y+1;
    end
    %
    % Forcing files
    %
    if makefrc==1
      %
      % Copy the file
      %
      frcname=[frc_prefix,'Y',num2str(Ymin),'M',num2str(M),nc_suffix];
      frcname2=[frc_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
      disp(['Create ',frcname2]) 
      eval(['!cp ',frcname,' ',frcname2]) 
      %
      % Change the time
      %
      nc=netcdf(frcname2,'write');
      time=nc{'sms_time'}(:)-365.*(Ymin-Y);%+datenum(Yorig,1,1);
      %[y,m,d,h,mi,s]=datevec(time);
      %dy=Ymin-Y;
      %y=y-dy;
      %time=datenum(y,m,d,h,mi,s)-datenum(Yorig,1,1);
      nc{'sms_time'}(:)=time;
      close(nc)
    end
    %
    % Bulk files
    %
    if makeblk==1
      %
      % Copy the file
      %
      blkname=[blk_prefix,'Y',num2str(Ymin),'M',num2str(M),nc_suffix];
      blkname2=[blk_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
      disp(['Create ',blkname2]) 
      eval(['!cp ',blkname,' ',blkname2]) 
      %
      % Change the time
      %
      nc=netcdf(blkname2,'write');
      time=nc{'bulk_time'}(:)-365.*(Ymin-Y);%+datenum(Yorig,1,1);
      %[y,m,d,h,mi,s]=datevec(time);
      %dy=Ymin-Y;
      %y=y-dy;
      %time=datenum(y,m,d,h,mi,s)-datenum(Yorig,1,1);
      nc{'bulk_time'}(:)=time;
      close(nc)
    end
  end
end
%---------------------------------------------------------------
% Create MC file
%---------------------------------------------------------------
%if createmc==1
%    if makefrc
%       make_mc_forcing('ECMWF',21600.,itolap);
%   end
%   if makeblk
%       make_mc_bulk('ECMWF',21600.,itolap);
%   end
%end
%
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
if makeplot==1
  disp(' ')
  disp('======================================================')
  disp(' Make a few plots...')
  slides=[1 12 24 36];
  if makeblk
    test_forcing(blkname,grdname,'tair',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'rhum',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'prate',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'uwnd',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'vwnd',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'wspd',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'radlw',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'radsw',slides,3,coastfileplot)
  end
  if makefrc
    figure
    test_forcing(frcname,grdname,'sustr',slides,3,coastfileplot)
    figure
    test_forcing(frcname,grdname,'svstr',slides,3,coastfileplot)  
    figure
    test_forcing(frcname,grdname,'Awave',slides,3,coastfileplot)
    figure
    test_forcing(frcname,grdname,'Dwave',slides,3,coastfileplot)
    figure
    test_forcing(frcname,grdname,'Pwave',slides,3,coastfileplot)
  end
end

