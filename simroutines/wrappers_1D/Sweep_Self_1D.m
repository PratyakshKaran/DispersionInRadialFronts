%-----------------------------------------------------------------------------------------------------------------------------------
% Wrapper [Variant - Production Sweep] for
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
warning('off','MATLAB:nearlySingularMatrix');
addpath('../functions');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% definings case combinations
namegeom =						{'CI','SI'};
eta =							[0.0E+0,	...
								1.00E-8,	3.33E-8,	1.00E-7,	3.33E-7,	...
								1.00E-6,	3.33E-6,	...
								1.00E-5,	3.33E-5,	1.00E-4,	3.33E-4,	...
								1.00E-3,	3.33E-3,	1.00E-2,	3.33E-2,	...
								1.00E-1,	3.33E-1,	1.00E+0,	3.33E+0,	...
								1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	...
								1.00E+3,	3.33E+3,	1.00E+4,	3.33E+4,	...
								1.00E+5,	3.33E+5,	1.00E+6		];
etaname =						{'0d0EP0',	...
								'1d00EM8',	'3d33EM8',	'1d00EM7',	'3d33EM7',	...
								'1d00EM6',	'3d33EM6',	...
								'1d00EM5',	'3d33EM5',	'1d00EM4',	'3d33EM4',	...
								'1d00EM3',	'3d33EM3',	'1d00EM2',	'3d33EM2',	...
								'1d00EM1',	'3d33EM1',	'1d00EP0',	'3d33EP0',	...
								'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	...
								'1d00EP3',	'3d33EP3',	'1d00EP4',	'3d33EP4',	...
								'1d00EP5',	'3d33EP5',	'1d00EP6'	};
Pe =							[1.00E+0,	3.33E+0,	1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	1.00E+3		];
Pename =						{'1d00EP0',	'3d33EP0',	'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	'1d00EP3'	};
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% obtaining grid skewness
Pe_b =							[1.00E+0,	1.00E+3];
eta_b =							[1.00E-3,	1.00E+3];
[Pe_mg_b, eta_mg_b] =			meshgrid(Pe_b,eta_b);
[Pe_mg, eta_mg] =				meshgrid(Pe,eta(2:end));
xgridskew_CI_b =				[3.10, 2.50; 2.45, 2.45];
xgridskew_SI_b =				[2.80, 2.45; 3.15, 2.50];	
xgridskew_CI_b_noDisp =			[3.10, 2.50];
xgridskew_SI_b_noDisp =			[2.80, 2.55];
xgridskew_CI =					interp2(Pe_mg_b,eta_mg_b,xgridskew_CI_b,Pe_mg,eta_mg,'spline');
xgridskew_SI =					interp2(Pe_mg_b,eta_mg_b,xgridskew_SI_b,Pe_mg,eta_mg,'spline');
xgridskew_CI_noDisp =			interp1(Pe_b,xgridskew_CI_b_noDisp,Pe,'pchip');
xgridskew_SI_noDisp =			interp1(Pe_b,xgridskew_SI_b_noDisp,Pe,'pchip');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
for igeom = 2
for iPe = 1
for ieta = 14
% for iPe = [1,5]
% for ieta = [14,22] [1,14,18,22,30]
%-----------------------------------------------------------------------------------------------------------------------------------
%	Input Parameters
simul.localonly =					1;					%	save output in a local-only folder
simul.localonlylevel =				1;					%	local-only folder location to save output data
%-----------------------------------------------------------------------------------------------------------------------------------
%	***high-Pe***
%	Simulation Parameters
%	Simulation Solver Parameters
	simul.solveGPU = 				0;					%	whether to solve linear system on GPU (1 = yes, 0 = no)
	simul.nx =						1+10;				%	number of grid points in the x-coordinate
	simul.nt =						1+50000;			%	number of time steps
	simul.nt2write =				324;				%	number of time steps to write to file
	simul.tstepprevs =				1;					%	number of previous time solutions utilized in time derivative approxmn
														%	(is effectively a measure of the order of FD approximation of time)
	simul.confmap.space =			0;					%	whether space is conformally mapped
	simul.confmap.time =			0;					%	whether time is conformally mapped
	simul.confmap.spacegridfromreal=0;					%	whether to map conformed space grid from real time grid
	simul.confmap.timegridfromreal=	0;					%	whether to map conformed time grid from real time grid
	if (eta(ieta)~=0)
		simul.xgridskew =			double(strcmp(namegeom{igeom},'CI'))*xgridskew_CI(ieta-1,iPe) +		...
									double(strcmp(namegeom{igeom},'SI'))*xgridskew_SI(ieta-1,iPe);
	else
		simul.xgridskew =			double(strcmp(namegeom{igeom},'CI'))*xgridskew_CI_noDisp(iPe) +	...
									double(strcmp(namegeom{igeom},'SI'))*xgridskew_SI_noDisp(iPe);
	end
														%	skewing parameter for xgrid [SI -> Pe end eta end]
														%	(higher than unity index - denser near start)
														%	(lower than unity index - denser near end)
	simul.tgridskew =				3.0;				%	(higher than unity index - denser near start)
														%	(lower than unity index - denser near end)
	if (simul.confmap.time == 0)
		simul.t2savemode =			'log';				%	mode of selecting time steps to write to file
	else
		simul.t2savemode =			'uniform';			%	mode of selecting time steps to write to file
	end
														%	uniform - uniformly distributed from start to end
														%	log - logarithmically distributed with more points near start time
														%	revlog - logarithmically distributed with more points near end time
														%	manual - manually specified simul.tstepsave (to be specified below)
	simul.STLtol =					1.00E-05;			%	tolerance of update diminution for Source Term Linearization Iterations
	simul.STLrelax =				1.00;				%	relaxation of solution update
	simul.errSTLmode =				1;					%	mode of computing STL error
														%	1 - standard deviation from last iteration solution
														%	2 -	standard deviation from last iteration solution times number of pts
														%	3 - maximum difference from last iteration
	simul.STexplicit =				0;					%	whether to explicitly solve reaction term (i.e. the source term ST)
	simul.iterateC =				1;					%	whether to iterate for solution of C
	if (simul.confmap.time == 0)
		simul.dt =					2.00E-4;				%	non-dimensionalized time step size
	else
		simul.dt =					((1.00E-2-eps)/(simul.nt-1));
														%	non-dimensionalized time step size
	end
	simul.compactgeom =				1;					%	make the geometry compact to optimize computation time
	simul.compactgeomfactotal =		25*(1+9*(igeom-1));	%	geometry compaction factor for front location
	simul.compactgeomsubfacwidth =	25;					%	geometry compaction factor for front width
	if (simul.confmap.space == 0)
		simul.limitdx =				1;					%	override specified point count to limit dx step size (1 = yes, 2 = no)
		simul.dxLmax =				1.00E-05;			%	maximum allowed dx at left end (enforced when simul.limitdx is 1)
		simul.dxRmax =				1.00E+03;			%	maximum allowed dx at right end (enforced when simul.limitdx is 1)
		simul.dxRmax_comprtoA =		1.00E-03;			%	maximum allowed dx at right end as compared to A
														%	(enforced when simul.limitdx is 1)
	end
	if (simul.confmap.time == 0)
		simul.limitdt =				1;					%	override specified point count to limit dt step size (1 = yes, 2 = no)
		simul.optimizetgrid =		1;					%	whether to optimize t-grid size and skewness
		simul.dtstartmax =			1.00E-14;			%	maximum allowed dt at start time (enforced when simul.limitdt is 1)
		simul.dtendmax =			1.00E-03;			%	maximum allowed dt at end time (enforced when simul.limitdt is 1)
	end
	simul.frontparamsrel =			1;					%	are front parameters presented normalized to Lx (1 = yes, 0 = no)
														%	(when yes, the front parameters are multiplied to Lx before proceeding
														%	to simulation)
	simul.injectparamsrel =			0;					%	are injection parameters presented normalized to Lx (1 = yes, 0 = no)
%	Simulation System Parameters
	simul.frcdispers = 				double(strcmp(namegeom{igeom},'PI'));	
														%	forced dispersion (applicable when the frame moves alongwith front)
	simul.justinject =				double(~strcmp(namegeom{igeom},'PI'));
														%	whether species A has just been injected at start time
														%	(if set to 0, generate a front, else, start injection at 
														%	injection plane/line/point)
	simul.isflowrate =				1;					%	is flowrate given (1) or inlet pressure head (0)
	simul.givenvel =				1;					%	0 - solve Darcy solver to get velocity field or 1 - take user-specified 
														%	velocity profile
	simul.nodiffusion = 			0;					%	whether to consider diffusion (0) or not (1)
	simul.casename =				[namegeom{igeom},'_Pe_',Pename{iPe},'_eta_',etaname{ieta},'_gamm_0d1'];
														%	name of case being solved
	simul.derivwarn =				true;				%	whether to warn before switching derivative direction
	simul.compMassC =				1;					%	whether to compute mass of product
	simul.compfront =				1;					%	whether to compute front properties
%	Geometric & Domain Parameters
	geomdom.isradial =				double(strcmp(namegeom{igeom},'CI'))+2*double(strcmp(namegeom{igeom},'SI'));			
														%	is the domain radial or linear
	if (simul.confmap.space == 0)
		geomdom.size.Lx =			1.00E+06;			%	non-dimensionalized x-coordinate length
	else
		geomdom.size.Lx =			1.0-eps;			%	non-dimensionalized x-coordinate length
	end
%	domain boundary conditions for transport - species concentration and its gradient -
%	1 - 	given constant net hydraulic head
%	2 - 	given flow rate (i.e. net hydraulic head gradient)
	geomdom.flowbc.left =			1;					%	near-end transport boundary conditions
	geomdom.flowbc.right =			1;					%	far-end transport boundary conditions
%	domain boundary conditions for transport - species concentration and its gradient - same meanings as for flow -
	geomdom.transpbc.left =			[1,1,0];			%	near-end transport boundary conditions
	geomdom.transpbc.right =		[1,1,0];			%	far-end transport boundary conditions
	geomdom.transpinit.rsm =		0;					%	initial condition taked as resumed solution of transport equation
%	domain permeability switch (see Flow & Transport Parameters for input variables) -
%	0 - 	uniform permeability
%	1 - 	specified permeability based on function of x and y-coordinates
%	2 -		stochastic permeability generated based on stochasticticity specifications (using a function call to permgen)
	geomdom.perm =					0;
%	Flow & Transport Parameters
	fltr.nondim.gamm =				1.00;				%	concentration of B at far end
	fltr.nondim.Pe =				Pe(iPe);			%	expected Peclet number
	fltr.flowbc.Hleft0 =			0.00;				%	pressure head boundary condition constant - left end
	fltr.flowbc.Hright0 =			0.00;				%	pressure head boundary condition constant - far-end
% 	fltr.flowbc.Q =					double(~strcmp(namegeom{igeom},'PI'));				
% 														%	flow rate
	fltr.transpbc.phileft0 =		[1.00,0.00,0.00];	%	species concentration boundary condition constant - left end
	fltr.transpbc.phiright0 =		fltr.nondim.gamm*	...
									[0.00,1.00,0.00];	%	species concentration boundary condition constant - right end
	fltr.transpbc.frcdispers = 		0.0;				%	forced velocity for disperpsion (i.e. horizontal speed of frame)
 	fltr.transpbc.inject.k1 =		1.00E-06;			%	transport central injection constant - central injection zone size
	fltr.transpbc.inject.singular = 1;					%	whether to implement injection transport b.c. at exactly r=0
 	fltr.init.front.k1 =			0.5;				%	reactive front constant - front zone advancement
	fltr.init.front.k2 =			1.00E-05;			%	reactive front constant - front zone step width
														%	central injection
	fltr.perm.Kperm0 = 				1.00;				%	permeability constant - characteristic permeability
	fltr.perm.Kvarrang =	 		0.00;				%	permeability constant - range of permeability variation
	fltr.eta =						eta(ieta);			%	expected longitudinal dispersivity factor
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Calling Preprocessor
	preprocess1D;
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Flow Solution
%	Obtaining Flow Solution by Darcy Solver
	if (simul.givenvel == 0)
		[soln.flow.H, soln.flow.vel, fltr.perm.Kperm] =		...
										solver1D_darcy(geomdom,simul,fltr,grd);
	else
		fltr.perm.Kperm =				fltr.Kperm;
		soln.flow.H =					zeros(simul.nx,1);
		soln.flow.vel.v =				double(geomdom.isradial~=0)*(1./(grd.x.^double(geomdom.isradial)));
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Saving Output
 	save([simul.outfoldrname,'/soln.mat']);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Transport Solution
 	solnreactive =					solver1D_reactive(geomdom,simul,fltr,grd,soln.flow);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Saving Mass and Front Properties
 	save([simul.outfoldrname,'/solnreactive.mat'],'solnreactive');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Displaying Progress
	disp(simul.casename);
%-----------------------------------------------------------------------------------------------------------------------------------
end
end
end
%-----------------------------------------------------------------------------------------------------------------------------------