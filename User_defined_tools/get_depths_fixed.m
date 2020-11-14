function [z]=get_depths_fixed(fname,gname,tindex,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the depths of the sigma levels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc=netcdf(gname);
h=nc{'h'}(:);
close(nc);
%
% open history file
%
nc=netcdf(fname);
zeta=squeeze(nc{'zeta'}(tindex,:,:));
theta_s=5.0;
theta_b=0.4;
Tcline =200;
hc=200;
hmin=25;
s_coord=2;

N=length(nc('s_rho'));
if s_coord==2
 hc=Tcline;
end
close(nc)
%
%
%
if isempty(zeta)
  zeta=0.*h;
end

vtype=type;
if (type=='u')|(type=='v')
  vtype='r';
end
z=zlevs(h,zeta,theta_s,theta_b,hc,N,vtype,s_coord);
if type=='u'
  z=rho2u_3d(z);
end
if type=='v'
  z=rho2v_3d(z);
end
return
