function [z]=get_depths(fname,gname,tindex,type)
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
hmorph=squeeze(nc{'hmorph'}(tindex,:,:));
if ~isempty(hmorph), h=hmorph; end
theta_s=nc.theta_s(:); 
if (isempty(theta_s))
%  disp('Rutgers version')
  theta_s=nc{'theta_s'}(:);
  theta_b=nc{'theta_b'}(:);
  Tcline=nc{'Tcline'}(:);
else 
%  disp('AGRIF/UCLA version');
  theta_b=nc.theta_b(:);
  Tcline=nc.Tcline(:);
  hc=nc.hc(:);
end
if (isempty(Tcline))
%  disp('UCLA version');
  hc=nc.hc(:);
else
  hmin=min(min(h));
  hc=min(hmin,Tcline);
end 
N=length(nc('s_rho'));
s_coord=1;
VertCoordType = nc.VertCoordType(:);
if isempty(VertCoordType)
  vtrans=nc{'Vtransform'}(:);
  if ~isempty(vtrans)
    s_coord=vtrans;
  end
elseif VertCoordType=='NEW', 
 s_coord=2;
end
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
