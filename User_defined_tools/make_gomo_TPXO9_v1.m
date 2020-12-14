% A wrapper for function TPXO2ROMS_v5pt1

%addpath(genpath('/storage0/group/Oceanography/JOHNL/myroms/PrepTides'))

%t0 can be any time of day not necessarily 00:00 
%lengthSim is the approximate anticipated length of model run (days) (f & u are computed at mid-point).
%fnGrid is your ROMS grid file.
%TPXO2ROMS will write to fnOut.
%The available choices for ROMSnames are:  'MM' 'MF' 'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4' '2N2' 'S1'
%For example, use ROMSnames={'M2' 'S2' 'MM'}; to use M2, S2, and MM only
clc; clear all; close all;
crocotools_param

%****************Modify to suit*************** 
% TIDE_START=17967;  %From ROMS *.in file (TIME_REF=-2)
t0=datenum(2019,8,10);  %TIDE_START in Matlab datenum
lengthSim=8;  %approximate length of model run in days
% fnGrid='E:\croco_tools\case\hato\prd_grid_1km.nc';
% fnOut=['.\roms_tide_' datestr(t0,'ddmmmyyyy') '.nc'];
ROMSnames={'M2','S2','N2','K2','K1','O1','P1','Q1','MM','MF','M4','MN4','MS4','2N2','S1'};
gomo_tidename = 'E:\croco_tools\case\lekima\tide_bc.nc';
gomo_gridname = 'E:\croco_tools\case\lekima\LEKIMA.grid.nc';
%********************************************

%Do not use existing filename as the number of harmonics may differ, causing untold grief. 
if exist(gomo_tidename,'file')
    delete(gomo_tidename)
%    error('File fnOut exists - must delete it or modify fnOut filename.')
end
TPXO2ROMS_GOMO(t0,ROMSnames,gomo_gridname,gomo_tidename,lengthSim,DATADIR)

%Say goodbye
disp('"May all our differences be finite"')
