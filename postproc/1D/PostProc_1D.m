% initial clean-up and set-up
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');

% what to post-process
postproc.iFlow = 0;						%	plot flow solution
postproc.iTransport = 1; 				%	movie transport solution
postproc.iTransportSS = 1; 				%	snapshots of transport solution
postproc.iTransportAmBEq = 0;			%	movie transport solution (theta)
postproc.iMassC = 1;					%	plot mass of C
postproc.iPDFC = 0;						%	movie PDF of concentration of C (phi_C)
postproc.iFront = 1;					%	plot front properties
postproc.iTermsGDE = 0;					%	plot terms of GDE (at origin and at front)
postproc.iTermsGDEmovie = 0;			%	create movie of terms of GDE
% post-processing option toggles
postproc.vistog = 'on';					%	plot windows visibility toggle
postproc.vistogmov = 'on';				%	plot windows visibility toggle (movies)
postproc.gridtog = 'off';				%	grid toggle
postproc.derivwarn = 'true';			%	derivative computation warnings display toggle
postproc.nTransportSS = 100;			%	number of transport solution of snapshots
postproc.nfourier = 1000;				%	number of terms in Fourier expansion (for analytical comparison)
postproc.iscomparison2D = 0;			%	whether comparing with 2D solution
postproc.iscomparison1D = 0;			%	whether comparing with 1D solution
postproc.isexpected = 0;				%	whether there are expected solutions
postproc.zoomin = 1;					%	zoom into the front location (movies)
postproc.plotmaxfrontR = 0;				%	plot maximum and front values of R using marter
postproc.plotfrontwidth = 1;			%	shade reaction zone width around the reaction front
postproc.timepausevid = 1.0E-5;			%	time for pausing in movie generation
postproc.obtainanalytical = 0;			%	obtain approximate analytical solutions
postproc.fracTransitn = 0.5;			%	factor for defining transition zone
postproc.forcesolnsimil = 1;			%	forcibly obtain similarity solution even if present from previous postprocessing
postproc.forceanalytical = 1;			%	forcibly obtain analytical front solution even if present from previous postprocessing
postproc.forcenumfront = 1;				%	forcibly obtain numerical front solution even if present from previous postprocessing
postproc.forcemassC = 1;				%	forcibly obtain mass of product even if present from previous postprocessing or from 
										%	during-simulation solution
postproc.maxxrange_log = 1E-10;			%	maximum range of x-axis for log scale (transport solution movie)
postproc.thresPedispsimil = 1.0E50;		%	maximum Pe to allow dispersion-dominated similarity solution (CI)
postproc.thresPediff = 5.0E0;			%	maximum Pe to use the Gamma function expression in diffusion-dominated theta-solution
% similarity solution (approximate analytical solution) computation parameters
postproc.simil.nx = 10001;				%	number of grid points on x-axis
postproc.simil.absmaxx = 1E4;			%	maximum admissible value of x
postproc.simil.power = 2.5;				%	x-grid skewness power
postproc.simil.rootfindrelax = 0.75;	%	relaxation in root-finding iterations
postproc.simil.tolrootfind = 1E-5;		%	tolerance for convergence of root-finding iterations
postproc.simil.barQthres = 10;			%	threshold of barQ for considering it to be large (dispersion-free scenario)
postproc.simil.leftBCtuned = 1;			%	whether to tune left-end BC of approximate solution of phi_A (conventional in De Wit)
postproc.simil.rootfindmeth = 'NR';		%	root-finding method - NR or STL
postproc.simil.maxrootfinditer = 1E3;	%	maximum number of admissible iterations for root-finding

% overriding front computation as per requirement
postproc.iFront =						postproc.iFront || postproc.iTermsGDE;
postproc.speciesname =					['A','B','C'];

% simulation output data structure
postproc.namegeom =						{'CI','SI'};
postproc.eta =							[0.0E+0,	...
										1.00E-8,	3.33E-8,	1.00E-7,	3.33E-7,	...
										1.00E-6,	3.33E-6,	1.00E-5,	3.33E-5,	...
										1.00E-4,	3.33E-4,							...
										1.00E-3,	3.33E-3,	1.00E-2,	3.33E-2,	...
										1.00E-1,	3.33E-1,	1.00E+0,	3.33E+0,	...
										1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	...
										1.00E+3,	3.33E+3,	1.00E+4,	3.33E+4,	...
										1.00E+5,	3.33E+5,							...
										1.00E+6];
postproc.etaname =						{'0d0EP0',	...
										'1d00EM8',	'3d33EM8',	'1d00EM7',	'3d33EM7',	...
										'1d00EM6',	'3d33EM6',	'1d00EM5',	'3d33EM5',	...
										'1d00EM4',	'3d33EM4',							...
										'1d00EM3',	'3d33EM3',	'1d00EM2',	'3d33EM2',	...
										'1d00EM1',	'3d33EM1',	'1d00EP0',	'3d33EP0',	...
										'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	...
										'1d00EP3',	'3d33EP3',	'1d00EP4',	'3d33EP4',	...
										'1d00EP5',	'3d33EP5',							...
										'1d00EP6'};
postproc.Pe =							[1.00E+0,	3.33E+0,	1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	1.00E+3		];
postproc.Pename =						{'1d00EP0',	'3d33EP0',	'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	'1d00EP3'	};

postproc.folderloc1Dparent =			'B:\LocalOnlyRawData_Level_1\1D';

for igeom = 2
for iPe = 1
for ieta = 14

% try complete solution
% try

% try to clear previous solution left-overs
try
	clear soln1D;
	close soln_analyze;
	clear front;
catch
end

% name of specific case
postproc.casename =						[postproc.namegeom{igeom},'_Pe_',postproc.Pename{iPe},		...
										'_eta_',postproc.etaname{ieta},'_gamm_0d1'];
% location of specific case data
postproc.folderloc1D =					[postproc.folderloc1Dparent,'/',postproc.casename];

% loading comparison solutions
if (postproc.iscomparison1D == 1)
	postproc.postproc.folderloccompare1Dparent =		...
										'../../../../Validation/1D_Code [Legacy_A]';
	postproc.casenamecompare1D =		['Radial_Injection',DomainSizeName{iDomainSize},'_Rescaled_Flow_Pe',Pename{iPe}];
	postproc.folderloccompare1D =		[postproc.postproc.folderloccompare1Dparent,'/',postproc.casenamecompare1D];
end
if (postproc.iscomparison2D == 1)
	postproc.folderloccompare2Dparent =					...
										'../../../../Validation/2D_Code [Legacy_A]';
	postproc.casenamecompare2D =		['Radial_Injection',DomainSizeName{iDomainSize},'_Pe',Pename{iPe}];
	postproc.folderloccompare2D =		[postproc.folderloccompare2Dparent,'/',postproc.casenamecompare2D];
end

% analyzing data
Analyzer_1D;

% plotting flow solution
if (postproc.iFlow == 1)
	plotFlow;
end

% movie of transport solution
if ((postproc.iTransport == 1) || (postproc.iTransportSS == 1))
	movTransport;
end

% plotting terms of theta equation
if (postproc.iTransportAmBEq == 1)
	movTransportAmB;
end

% plotting mass of product C
if (postproc.iMassC == 1)
	plotMassC;
end

% movie of PDF of concentration of C (phi_C)
if (postproc.iPDFC == 1)
	movPDFC;
end

% plotting front properties
if (postproc.iFront == 1)
	plotFront;
end

% plotting terms of GDE
if (postproc.iTermsGDE == 1)
	plotTermsGDE;
end

% movie of terms of GDE
if (postproc.iTermsGDEmovie == 1)
	movTermsGDE;
end

% catch
% end

end
end
end