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
addpath('../functions');
%-----------------------------------------------------------------------------------------------------------------------------------

namegeom =		{'Plane','Line','Point'};
eta =			[0, 0.01, 0.1, 1, 10];
etaname =		{'No', 'VLow', 'Low', 'Med', 'High'};
gamm =			[0.1, 0.33, 1, 3, 10];
gammname =		{'0d1', '0d33', '1', '3', '10'};
Q =				[1, 10, 25, 100, 1000, 1E10];
Qname =			{'1','10','25','100','1000','Inf'};
Pe =			[1E0, 1E1, 1E2, 1E3];
Pename =		{'1E0', '1E1', '1E2', '1E3'};
Da =			[1E-1, 1E0, 1E1, 1E2];
Daname =		{'1Em1', '1E0', '1E1', '1E2'};

for iPe = 2
for iDa = 3
for igeom = 2
for ieta = [1,4]
for igamm = 4
for iQ = 1
	
%-----------------------------------------------------------------------------------------------------------------------------------
%	Input Parameters
simul.localonly =					1;					%	save output in a local-only folder
simul.localonlylevel =				3;					%	local-only folder location to save output data
%-----------------------------------------------------------------------------------------------------------------------------------
%	***high-Pe***
%	Simulation Parameters
%	Simulation Solver Parameters
	simul.nx =						1+2000*(1+double(strcmp(namegeom{igeom},'Plane')));				
														%	number of grid points in the x-coordinate
	simul.nt =						1+10000*(1+99*double(strcmp(namegeom{igeom},'Point')));		
														%	number of time steps
	simul.solveGPU = 				0;					%	whether to solve linear system on GPU (1 = yes, 0 = no)
	simul.tstepprevs =				1;					%	number of previous time solutions utilized in time derivative approxmn
														%	(is effectively a measure of the order of FD approximation of time)
	simul.nt2write =				1+(simul.nt-1)/(100*(1+99*double(strcmp(namegeom{igeom},'Point'))));
														%	number of time steps to write to file
	simul.xgridskew =				1.0+0.25*double(strcmp(namegeom{igeom},'Line'))+1.50*double(strcmp(namegeom{igeom},'Point'));	
														%	skewinf paratial for xgrid 
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
	simul.maxdt =					1.00;				%	allowed maximum non-dimensionalized time step size (skewed t-grid)
	simul.rescalet =				0;					%	whether to rescale t so that the maximum time step size is simul.maxdt
	simul.dt =						0.01;				%	non-dimensionalized time step size
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
%	Geometric & Domain Parameters
	geomdom.isradial =				igeom-1;			%	is the domain radial or linear
	geomdom.char.Lc =				1.0E-06;			%	characteristic length
	geomdom.char.Hc =				1.0E-08;			%	characteristic pressure head (hydraulic head)
	geomdom.char.vc =				geomdom.char.Hc/geomdom.char.Lc;	
														%	characteristic flow velocity (discharge rate)
	geomdom.char.tc =				geomdom.char.Lc/geomdom.char.vc;	
														%	characteristic time scale
	geomdom.size.Lx =				500+double(strcmp(namegeom{igeom},'Plane'))*(200-500);			
														%	non-dimensionalized x-coordinate length
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
 	fltr.transpbc.inject.k1 =		1.0E-05;			%	transport central injection constant - central injection zone size
 	fltr.init.front.k1 =			geomdom.size.Lx*0.5;%	reactive front constant - front zone advancement
	fltr.init.front.k2 =			1.0E-10;			%	reactive front constant - front zone step width
														%	central injection
	fltr.perm.Kperm0 = 				1.00;				%	permeability constant - characteristic permeability
	fltr.perm.Kvarrang =	 		0.75;				%	permeability constant - range of permeability variation
	fltr.eta_exp =					eta(ieta);			%	expected longitudinal dispersivity factor
	fltr.betaT_exp =				0.10;				%	expected transverse dispersivity factor
%	Computed Flow & Transport Parameters
	fltr.nondim.diffus =			(geomdom.char.Lc^2)/(geomdom.char.tc*fltr.nondim.Pe_exp);
														%	diffusivity
	fltr.nondim.kreact =			(fltr.nondim.Da_exp*fltr.nondim.diffus)/(fltr.phi0*geomdom.char.Lc^2);
														%	reaction rate
	fltr.nondim.Pe =				(geomdom.char.Lc^2)/(geomdom.char.tc*fltr.nondim.diffus);		
														%	Peclet number
	fltr.nondim.Da =				(fltr.nondim.kreact*fltr.phi0*geomdom.char.Lc^2)/(fltr.nondim.diffus);			
														%	Damkohler Number
	fltr.betaL_exp =				fltr.eta_exp*fltr.nondim.Pe;		
														%	expected normalized longitudinal dispersivity
	fltr.eta =						fltr.betaL_exp/fltr.nondim.Pe;
	fltr.alphaL =					(fltr.nondim.diffus*fltr.betaL_exp)/geomdom.char.vc;			
														%	longitudinal dispersivity
	fltr.alphaT =					fltr.alphaL*fltr.betaT_exp;			
														%	transverse dispersivity
	fltr.betaL =					(fltr.alphaL*geomdom.char.vc)/fltr.nondim.diffus;	
														%	normalized longitudinal dispersivity
	fltr.betaT =					fltr.alphaT/fltr.alphaL;
														%	transverse dispersivity factor
%	Computed Geometric & Domain Parameters
	geomdom.size.a =				fltr.transpbc.inject.k1;
	geomdom.size.A =				geomdom.size.Lx;
	if (simul.isflowrate == 1)
		if (geomdom.isradial == 0)
			fltr.flowbc.Hleft0 =	fltr.flowbc.Hright0+fltr.flowbc.Q*geomdom.size.Lx;
		elseif (geomdom.isradial == 1)
			fltr.flowbc.Hleft0 =	fltr.flowbc.Hright0-fltr.flowbc.Q*(log(geomdom.size.a+1E-100)-log(geomdom.size.A));
		else
			fltr.flowbc.Hleft0 =	fltr.flowbc.Hright0+fltr.flowbc.Q*	...
									((geomdom.size.A-(geomdom.size.a+1E-100))/geomdom.size.A*(geomdom.size.a+1E-100));
		end
	end

%	Input Variables with Switches
%	Generating array of time steps to write solution to file
	if (strcmp(simul.t2savemode,'uniform'))
		simul.tstepsave =			floor(linspace(1,simul.nt,simul.nt2write));	
	elseif (strcmp(simul.t2savemode,'log'))
		simul.tstepsave =			floor(exp(linspace(log(1),log(simul.nt),simul.nt2write)));
	elseif (strcmp(simul.t2savemode,'revlog'))
		simul.tstepsave =			floor(exp(linspace(log(1),log(simul.nt),simul.nt2write)));
		simul.tstepsave =			simul.nt-simul.tstepsave;
		simul.tstepsave =			flip(simul.tstepsave);
	elseif (strcmp(simul.t2savemode,'bothlog'))
		simul.tstepsave =			floor(exp(linspace(log(1),log(simul.nt),floor(simul.nt2write/2))));
		simul.tstepsave1 =			floor(exp(linspace(log(1),log(simul.nt),ceil(simul.nt2write/2))));
		simul.tstepsave1 =			simul.nt-simul.tstepsave1;
		simul.tstepsave1 =			flip(simul.tstepsave1);
		simul.tstepsave = 			[simul.tstepsave,simul.tstepsave1];
		simul.tstepsave = 			sort(simul.tstepsave);
	elseif (strcmp(simul.t2savemode,'manual'))
		simul.tstepsave =			linspace(1,100,floor(simul.nt/3));
		simul.tstepsave =			[simul.tstepsave,linspace(101,1000,floor(simul.nt/3))];
		simul.tstepsave =			[simul.tstepsave,linspace(1001,10000,floor(simul.nt/3))];
		simul.nt2write =			length(simul.tstepsave);
	end
	if (simul.tstepsave(1) < 1)
		simul.tstepsave(1) =		1;
	end
	if (simul.tstepsave(end) > simul.nt)
		simul.tstepsave(end) =		simul.nt;
	end	
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
%	Creating Output Folder
	if (simul.localonly == 0)
		simul.outfoldrname =		['../../outdata/1D/',simul.casename];
	else
		simul.outfoldrname =		['../../../../LocalOnlyRawData_Level_',num2str(simul.localonlylevel),'/1D/',simul.casename];
	end
	mkdir('../../outdata/1D');
	mkdir(simul.outfoldrname);
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
%	Grid Generation
	if (geomdom.isradial == 0)
		leftend =					0.0;
		rightend =					geomdom.size.Lx;
	else
		leftend =					geomdom.size.a;
		rightend =					geomdom.size.A;
	end
	grd.x =							leftend+(rightend-leftend)*linspace(0.0,1.0,simul.nx).^simul.xgridskew;
	starttime =						0.0;
	endtime =						simul.dt*(simul.nt-1);
	grd.t =							starttime+(endtime-starttime)*linspace(0.0,1.0,simul.nt).^simul.tgridskew;
	maxdt =							max([grd.t(2)-grd.t(1),grd.t(simul.nt)-grd.t(simul.nt-1)]);
	if (simul.rescalet == 1)
		grd.t =						grd.t*(simul.dt/maxdt);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Permeability Generation
	fltr.Kperm =					fltr.perm.Kperm0*zeros(simul.nx,1);
	if (geomdom.perm == 0)			% uniform permeability
		fltr.Kperm =				fltr.perm.Kperm0*ones(simul.nx,1);
	elseif (geomdom.perm == 1)		% permeability based on a function
		for ix = 1:nx
			fltr.Kperm(ix) =		fltr.perm.Kperm0*(1.0+fltr.perm.Kvarrang*sin(4.0*pi*(grd.x(ix)/geomdom.size.Lx)));
		end
	elseif (geomdom.perm == 2)		% permeability generated using a stochasticity function
		fltr.Kperm =				permgen(fltr.perm,grd,geomdom);
	end
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
%	Initial Front Generation
	fltr.init.phiA =				zeros(simul.nx,1);
	fltr.init.phiB =				zeros(simul.nx,1);
	fltr.init.phiC =				zeros(simul.nx,1);
	if (simul.justinject == 1)
		fltr.init.phiA(1) =			1.0;
		fltr.init.phiB(1:simul.nx)=	fltr.nondim.gamm;
	else
		fltr.init.phiA(1:simul.nx)=	0.5*(1.0-tanh((grd.x(1:simul.nx)-fltr.init.front.k1)/fltr.init.front.k2));
		fltr.init.phiB(1:simul.nx)=	fltr.nondim.gamm*0.5*(1.0+tanh((grd.x(1:simul.nx)-fltr.init.front.k1)/fltr.init.front.k2));
	end
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