%-----------------------------------------------------------------------------------------------------------------------------------
% Wrapper [Variant - Template] for
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
namegeom =		{'Plane','Line','Point'};
eta =			[0, 0.01, 0.1, 1, 10, 1E5];
etaname =		{'No', 'VLow', 'Low', 'Med', 'High', '1E5'};
gamm =			[1, 0.33, 3, 0.1, 10];
gammname =		{'1', '0d33', '3', '0d1', '10'};
Q =				[1, 10, 25, 100, 1000, 0.1];
Qname =			{'1','10','25','100','1000', '1Em1'};
Pe =			[1E0, 1E1, 1E2, 1E3, 0.1];
Pename =		{'1E0', '1E1', '1E2', '1E3', '1Em1'};
Da =			[1E-1, 1E0, 1E1, 1E2];
Daname =		{'1Em1', '1E0', '1E1', '1E2'};
%-----------------------------------------------------------------------------------------------------------------------------------

for iPe = 1:length(Pe)
for iDa = 1:length(Da)
for igeom = 1:1
for ieta = length(eta):length(eta)
for igamm = 1:length(gamm)
for iQ = 1:1+(igeom-1)*(length(Q)-1)
%-----------------------------------------------------------------------------------------------------------------------------------
%	Input Parameters
simul.localonly =					1;					%	save output in a local-only folder
simul.localonlylevel =				1;					%	local-only folder location to save output data
%-----------------------------------------------------------------------------------------------------------------------------------
%	***high-Pe***
%	Simulation Parameters
%	Simulation Solver Parameters
	simul.nx =						50001;				
														%	number of grid points in the x-coordinate
	simul.nt =						10001;		
														%	number of time steps
	simul.solveGPU = 				0;					%	whether to solve linear system on GPU (1 = yes, 0 = no)
	simul.tstepprevs =				1;					%	number of previous time solutions utilized in time derivative approxmn
														%	(is effectively a measure of the order of FD approximation of time)
	simul.nt2write =				101;				%	number of time steps to write to file
	simul.xgridskew =				2.0;				%	skewinf paratial for xgrid 
														%	(higher than unity index - denser near start)
														%	(lower than unity index - denser near end)
	simul.tgridskew =				1.0;				%	skewinf paratial for xgrid
														%	(higher than unity index - denser near start)
														%	(lower than unity index - denser near end)
	simul.t2savemode =				'uniform';			%	mode of selecting time steps to write to file
														%	uniform - uniformly distributed from start to end
														%	log - logarithmically distributed with more points near start time
														%	revlog - logarithmically distributed with more points near end time
														%	manual - manually specified simul.tstepsave (to be specified below)
	simul.STLtol =					1.00E-04;			%	tolerance of update diminution for Source Term Linearization Iterations
	simul.STLrelax =				1.00;				%	relaxation of solution update
	simul.errSTLmode =				1;					%	mode of computing STL error
														%	1 - standard deviation from last iteration solution
														%	2 -	standard deviation from last iteration solution times number of pts
														%	3 - maximum difference from last iteration
	simul.STexplicit =				0;					%	whether to explicitly solve reaction term (i.e. the source term ST)
	simul.dt =						0.01;	
														%	non-dimensionalized time step size
	simul.compactgeom =				0;					%	make the geometry compact to optimize computation time
	simul.limitdx =					1;					%	override specified point count to limit dx step size (1 = yes, 2 = no)
	simul.dxLmax =					0.025;				%	maximum allowed dx at left end (enforced when simul.limitdx is 1)
	simul.dxRmax =					2.500;				%	maximum allowed dx at right end (enforced when simul.limitdx is 1)
	simul.limitdt =					1;					%	override specified point count to limit dt step size (1 = yes, 2 = no)
	simul.dtstartmax =				0.005;				%	maximum allowed dt at start time (enforced when simul.limitdt is 1)
	simul.dtendmax =				0.500;				%	maximum allowed dt at end time (enforced when simul.limitdt is 1)
	simul.frontparamsrel =			0;					%	are front parameters presented normalized to Lx (1 = yes, 0 = no)
														%	(when yes, the front parameters are multiplied to Lx before proceeding
														%	to simulation)
	simul.injectparamsrel =			0;					%	are injection parameters presented normalized to Lx (1 = yes, 0 = no)
%	Simulation System Parameters
	simul.frcdispers = 				double(strcmp(namegeom{igeom},'Plane'));	
														%	forced dispersion (applicable when the frame moves alongwith front)
	simul.justinject =				double(~strcmp(namegeom{igeom},'Plane'));
														%	whether species A has just been injected at start time
														%	(if set to 0, generate a front, else, start injection at 
														%	injection plane/line/point)
	simul.isflowrate =				1;					%	is flowrate given (1) or inlet pressure head (0)
	simul.givenvel =				0;					%	0 - solve Darcy solver to get velocity field or 1 - take user-specified 
														%	velocity profile
	simul.nodiffusion = 			0;					%	whether to consider diffusion (0) or not (1)
	simul.casename =				['DeWit_2020_',namegeom{igeom},'Inj_',etaname{ieta},	...
									'Dispers_',gammname{igamm},'gamm_',Qname{iQ},'Q_',		...
									Pename{iPe},'Pe_',Daname{iDa},'Da'];
														%	name of case being solved
	simul.derivwarn =				true;				%	whether to warn before switching derivative direction
	simul.compMassC =				0;					%	whether to compute mass of product
	simul.compfront =				0;					%	whether to compute front properties
%	Geometric & Domain Parameters
	geomdom.isradial =				igeom-1;			%	is the domain radial or linear
	geomdom.char.Lc =				1.0E-06;			%	characteristic length
	geomdom.char.Hc =				1.0E-08;			%	characteristic pressure head (hydraulic head)
	geomdom.char.vc =				geomdom.char.Hc/geomdom.char.Lc;	
														%	characteristic flow velocity (discharge rate)
	geomdom.char.tc =				geomdom.char.Lc/geomdom.char.vc;	
														%	characteristic time scale
	geomdom.size.Lx =				250;					%	non-dimensionalized x-coordinate length
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
	fltr.phi0 =						1.00E-20;			%	species concentration scale
	fltr.nondim.gamm =				gamm(igamm);		%	concentration of B at far end
	fltr.nondim.Pe_exp =			Pe(iPe);			%	expected Peclet number
	fltr.nondim.Da_exp =			Da(iDa);			%	expected Damkohler Number
	fltr.flowbc.Hleft0 =			0.00;				%	pressure head boundary condition constant - left end
	fltr.flowbc.Hright0 =			0.00;				%	pressure head boundary condition constant - far-end
	fltr.flowbc.Q =					Q(iQ)*double(~strcmp(namegeom{igeom},'Plane'));				
														%	flow rate
	fltr.transpbc.phileft0 =		[1.00,0.00,0.00];	%	species concentration boundary condition constant - left end
	fltr.transpbc.phiright0 =		fltr.nondim.gamm*	...
									[0.00,1.00,0.00];	%	species concentration boundary condition constant - right end
	fltr.transpbc.frcdispers = 		1.0;				%	forced velocity for disperpsion (i.e. horizontal speed of frame)
 	fltr.transpbc.inject.k1 =		1.0E-02;			%	transport central injection constant - central injection zone size
 	fltr.init.front.k1 =			0.5;				%	reactive front constant - front zone advancement
	fltr.init.front.k2 =			1.0E-10;			%	reactive front constant - front zone step width
														%	central injection
	fltr.perm.Kperm0 = 				1.00;				%	permeability constant - characteristic permeability
	fltr.perm.Kvarrang =	 		0.75;				%	permeability constant - range of permeability variation
	fltr.eta_exp =					eta(ieta);			%	expected longitudinal dispersivity factor
	fltr.betaT_exp =				0.10;				%	expected transverse dispersivity factor
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Calling Preprocessor
	preprocess1D;
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Flow Solution
	[soln.flow.H, soln.flow.vel, fltr.perm.Kperm] =		...
									solver1D_darcy(geomdom,simul,fltr,grd);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Saving Output
 	save([simul.outfoldrname,'/soln.mat']);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Transport Solution
 	phi =							solver1D_reactive(geomdom,simul,fltr,grd,soln.flow);
%-----------------------------------------------------------------------------------------------------------------------------------
end
end
end
end
end
end
%-----------------------------------------------------------------------------------------------------------------------------------