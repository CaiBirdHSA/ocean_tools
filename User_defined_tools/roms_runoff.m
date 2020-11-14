%% begin
%huanghe; haihe; liaohe;
lon=[119.333,117.896,121.800];
lat=[37.788, 38.924, 40.838 ];
N=15;
flux=[100, 200, 300 ];          %unit=10.^8 m3/year
mud =[30,   40, 50 ];           %unit=kg/m3
flux=[20,21,22; 10,11,12; 5,6,7];

year=2010;
days=[1:3];
direction=[170,170,85];
nriver=length(lon);
nmonth=12;
vshape=ones(nriver,N)/N;
%transport=ones(nriver,nmonth);
transport=repmat(flux'/12,1,12)*10.^8/365/24/3600;
temp=ones(nriver,N,nmonth)*8;
salt=ones(nriver,N,nmonth)*15;
rivermud =repmat(mud', 1, N, nmonth);

%% read grid file
gridn = './bohai_grd.nc';
lon_rho=ncread(gridn,'lon_rho');
lat_rho=ncread(gridn,'lat_rho');
msk_rho=ncread(gridn,'mask_rho');

for i=1:length(lon)
    xi(i)=find(lon(i)<lon_rho(:,1),1,'first');
    et(i)=find(lat(i)<lat_rho(1,:),1,'first');
end

%% write to nc
xtype= 'NC_DOUBLE';
river_time=12;
river=length(lon);

ncid =netcdf.create(['./bohai_',num2str(year),'runoff.nc'],'64BIT_OFFSET');
nx  =netcdf.defDim(ncid,'river_time',river_time);    ny =    netcdf.defDim(ncid,'s_rho',N);
nz =netcdf.defDim(ncid,'river',river);

riverid=netcdf.defVar(ncid,'river',xtype, [nz]);
timeid= netcdf.defVar(ncid,'river_time',xtype, [nx]);
xpid  = netcdf.defVar(ncid,'river_Xposition',xtype, [nz]);
ypid  = netcdf.defVar(ncid,'river_Eposition',xtype, [nz]);
drid  = netcdf.defVar(ncid,'river_direction',xtype, [nz]);
shid  = netcdf.defVar(ncid,'river_Vshape',xtype, [nz ny]);
trid  = netcdf.defVar(ncid,'river_transport',xtype, [nz nx]);
teid  = netcdf.defVar(ncid,'river_temp'     ,xtype, [nz ny nx]);
said  = netcdf.defVar(ncid,'river_salt'     ,xtype, [nz ny nx]);
muid  = netcdf.defVar(ncid,'river_mud_'     ,xtype, [nz ny nx]);

netcdf.putAtt(ncid, timeid,'long_name', 'river runoff time');
netcdf.putAtt(ncid, timeid,'units',['days since ',num2str(year),'-01-01 00:00:00']);
netcdf.putAtt(ncid, timeid,'calendar', 'gregorian');

netcdf.endDef(ncid)
netcdf.putVar(ncid,riverid,[1:river]);
netcdf.putVar(ncid,timeid,days);
netcdf.putVar(ncid,xpid,xi);
netcdf.putVar(ncid,ypid,et);
netcdf.putVar(ncid,drid,direction);
netcdf.putVar(ncid,shid,vshape);
netcdf.putVar(ncid,trid,transport);
netcdf.putVar(ncid,teid,temp);
netcdf.putVar(ncid,said,salt);
netcdf.putVar(ncid,muid,rivermud);
netcdf.close(ncid);