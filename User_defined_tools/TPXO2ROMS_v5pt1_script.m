% A wrapper for function TPXO2ROMS_v5pt1
%addpath(genpath('/storage0/group/Oceanography/JOHNL/myroms/PrepTides'))

%t0 can be any time of day not necessarily 00:00 
%lengthSim is the approximate anticipated length of model run (days) (f & u are computed at mid-point).
%fnGrid is your ROMS grid file.
%TPXO2ROMS will write to fnOut.
%The available choices for ROMSnames are:  'MM' 'MF' 'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4' '2N2' 'S1'
%For example, use ROMSnames={'M2' 'S2' 'MM'}; to use M2, S2, and MM only
clear

%****************Modify to suit*************** 
% TIDE_START=17967;  %From ROMS *.in file (TIME_REF=-2)
% t0=TIDE_START+datenum(1968,5,23);  %TIDE_START in Matlab datenum
date0='2020-09-01 00:00:00';
t0=datenum(date0);
lengthSim=6;  %approximate length of model run in days
fnGrid='E:\croco_tools\croco_tools-forcast\CROCO_FILES\bohai_grid.nc';
fnOut=['E:\croco_tools\croco_tools-forcast\CROCO_FILES\frc_TPXO9_' datestr(t0,'yyyymmdd') '.nc'];
ROMSnames={'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4','2N2'};
%********************************************

%Do not use existing filename as the number of harmonics may differ, causing untold grief. 
if exist(fnOut,'file')
   error('File fnOut exists - must delete it or modify fnOut filename.')
end
TPXO2ROMS_v5pt1(t0,ROMSnames,fnGrid,fnOut,lengthSim)

%Say goodbye
disp('"May all our differences be finite"')
