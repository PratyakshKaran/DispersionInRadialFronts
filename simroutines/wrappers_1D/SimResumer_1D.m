%-----------------------------------------------------------------------------------------------------------------------------------
% Simulation Resumer for
% 1D Transient Reactive Transport Solver for Planar/Cylindrical/Spherical Injection using Finite Difference Method (FDM) &
% 1D Steady State Darcy Flow Solver corresponding for Planar/Cylindrical/Spherical Injection using Finite Difference Method (FDM)
% [Non-Function Program]
% 
%	written by --
%	uddipta ghosh (uddipta.ghosh@iitgn.ac.in)
%	pratyaksh karan (pratyakshkaran@gmail.com)
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	initiating cleaned matlab sessions
clc;
clear;
close all;
fclose all;
addpath('../functions');
warning('off','MATLAB:nearlySingularMatrix');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Input Parameters
%	Preparing Resumption Point
fldrrsm =						'D:\Pratyaksh\LocalOnlyRawData_Level_3\1D\SI_Pe_1d00EP3_eta_3d33EM3';
geomdom.transpinit.fldrrsm =	fldrrsm;
load([geomdom.transpinit.fldrrsm,'/soln.mat']);
geomdom.transpinit.fldrrsm =	fldrrsm;
geomdom.transpinit.rsm =		1;				%	initial condition taked as resumed solution of transport equation
geomdom.transpinit.datfl =		dir([geomdom.transpinit.fldrrsm,'/*.dat']);
geomdom.transpinit.datfl =		size(geomdom.transpinit.datfl,1);
geomdom.transpinit.wrcnt =		geomdom.transpinit.datfl/3-1;
simul.tstepsave = 				unique(simul.tstepsave);
geomdom.transpinit.rsmit =		simul.tstepsave(geomdom.transpinit.wrcnt);
geomdom.transpinit.rsmitprev =	simul.tstepsave(geomdom.transpinit.wrcnt-1);
												%	time step from where to resume
geomdom.transpinit.rsmloc =		geomdom.transpinit.fldrrsm;
												%	location from where the initial conditions are to be taken
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Loading Resumption Point
phiAcurrent =					load([geomdom.transpinit.rsmloc,'/phiA_',num2str(geomdom.transpinit.rsmit),'.dat']);
phiBcurrent =					load([geomdom.transpinit.rsmloc,'/phiB_',num2str(geomdom.transpinit.rsmit),'.dat']);
phiCcurrent =					load([geomdom.transpinit.rsmloc,'/phiC_',num2str(geomdom.transpinit.rsmit),'.dat']);
fltr.init.phiA =				phiAcurrent';
fltr.init.phiB =				phiBcurrent';
fltr.init.phiC =				phiCcurrent';
phiAprev =						load([geomdom.transpinit.rsmloc,'/phiA_',num2str(geomdom.transpinit.rsmitprev),'.dat']);
phiBprev =						load([geomdom.transpinit.rsmloc,'/phiB_',num2str(geomdom.transpinit.rsmitprev),'.dat']);
phiCprev =						load([geomdom.transpinit.rsmloc,'/phiC_',num2str(geomdom.transpinit.rsmitprev),'.dat']);
fltr.init.phiAprev =			phiAprev';
fltr.init.phiBprev =			phiBprev';
fltr.init.phiCprev =			phiCprev';
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Species Concentration Solution
soln.transp.phi =								solver1D_reactive(geomdom,simul,fltr,grd,soln.flow);
disp("Solved species evolution");
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
disp(['Finished Realization ',num2str(irealiz)]);
%-----------------------------------------------------------------------------------------------------------------------------------