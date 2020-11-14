function data=ext_data_OGCM(nc,X,Y,vname,tndx,lon,lat,k,Roa,interp_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extrapole one horizontal ECCO (or Data) slice on a CROCO grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% extrapolation parameters
%
default=0;
if strcmp(vname,'SAVE') || strcmp(vname,'salt') || strcmp(vname,'salinity')
  default=32.5;
end
%
% Get the CROCO grid extension + a little margin (~ 2 data grid points)
%
dx=max(abs(gradient(X)));
dy=max(abs(gradient(Y)));
dl=2*max([dx dy]);
%
lonmin=min(min(lon))-dl;
lonmax=max(max(lon))+dl;
latmin=min(min(lat))-dl;
latmax=max(max(lat))+dl;
%
% Extract a data subgrid
%
j=find(Y>=latmin & Y<=latmax);
i1=find(X-360>=lonmin & X-360<=lonmax);
i2=find(X>=lonmin & X<=lonmax);
i3=find(X+360>=lonmin & X+360<=lonmax);
if ~isempty(i2)
  x=X(i2);
else
  x=[];
end
if ~isempty(i1)
  x=cat(2,X(i1)-360,x);
end
if ~isempty(i3)
  x=cat(2,x,X(i3)+360);
end
y=Y(j);
%
%  Get dimensions
%
ndims=length(dim(nc{vname}));
%
% Get data (Horizontal 2D matrix)
%
if ~isempty(i2)
  if ndims==2
    data=squeeze(nc{vname}(j,i2));
  elseif ndims==3
    data=squeeze(nc{vname}(tndx,j,i2));
  elseif ndims==4
    data=squeeze(nc{vname}(tndx,k,j,i2));
  else
    error(['Bad dimension number ',num2str(ndims)])
  end
else
  data=[];
end
if ~isempty(i1)
  if ndims==2
    data=cat(2,squeeze(nc{vname}(j,i1)),data);
  elseif ndims==3
    data=cat(2,squeeze(nc{vname}(tndx,j,i1)),data);
  elseif ndims==4
    data=cat(2,squeeze(nc{vname}(tndx,k,j,i1)),data);
  else
    error(['Bad dimension number ',num2str(ndims)])
  end
end
if ~isempty(i3)
  if ndims==2
    data=cat(2,data,squeeze(nc{vname}(j,i3)));
  elseif ndims==3
    data=cat(2,data,squeeze(nc{vname}(tndx,j,i3)));
  elseif ndims==4
    data=cat(2,data,squeeze(nc{vname}(tndx,k,j,i3)));
  else
    error(['Bad dimension number ',num2str(ndims)])
  end
end
%
% Perform the extrapolation
%
[data,interp_flag]=get_missing_val(x,y,data,NaN,Roa,default);
%
% Interpolation on the CROCO grid
%
if interp_flag==0
  data=interp2(x,y,data,lon,lat,'nearest');
else
  data=interp2(x,y,data,lon,lat,interp_method);
end
%
return
