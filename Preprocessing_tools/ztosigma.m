function vnew = ztosigma(var,z,depth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function  vnew = ztosigma(var,z,depth)
%
% This function transform a variable from z to sigma coordinates
%    warning: the variable must be in the form: var(k,j,i)
%
% On Input:
%
%    var     Variable z (3D matrix).
%    z       Sigma depths (m) of RHO- or W-points (3D matrix).
%    depth   z depth (vector; meters, negative).
% 
% On Output:
%
%    vnew    Variable sigma (3D matrix).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ns,Mp,Lp]=size(z);
[Nz]=length(depth);
%
% Find the grid position of the nearest vertical levels
%
for ks=1:Ns
  sigmalev=squeeze(z(ks,:,:));
  thezlevs=0.*sigmalev;
  for kz=1:Nz
    thezlevs(sigmalev>depth(kz))=thezlevs(sigmalev>depth(kz))+1;
  end
  if max(max(thezlevs))>=Nz | min(min(thezlevs))<=0
    disp(['min sigma level = ',num2str(min(min(min(z)))),...
           ' - min z level = ',num2str(min(depth))])
    disp(['max sigma level = ',num2str(max(max(max(z)))),...
           ' - max z level = ',num2str(max(depth))])
  end
  [imat,jmat]=meshgrid((1:Lp),(1:Mp));
  pos=Nz*Mp*(imat-1)+Nz*(jmat-1)+thezlevs;
  z1=depth(thezlevs);
  z2=depth(thezlevs+1);
  v1=var(pos);
  v2=var(pos+1);
  vnew(ks,:,:)=(((v1-v2).*sigmalev+v2.*z1-v1.*z2)./(z1-z2));
end
return
