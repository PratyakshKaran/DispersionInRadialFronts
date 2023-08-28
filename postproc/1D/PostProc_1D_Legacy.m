% initial clean-up and set-up
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');

% what to post-process
iFlow = 0;							%	plot flow solution
iTransport = 1; 					%	movie transport solution
iTransportSS = 1; 					%	snapshots of transport solution
iTransportAmBEq = 0;				%	movie transport solution (theta)
iMassC = 1;							%	plot mass of C
iPDFC = 0;							%	movie PDF of concentration of C (phi_C)
iFront = 1;							%	plot front properties
iTermsGDE = 0;						%	plot terms of GDE (at origin and at front)
iTermsGDEmovie = 0;					%	create movie of terms of GDE
% post-processing option toggles
vistog = 'on';						%	plot windows visibility toggle
vistogmov = 'on';					%	plot windows visibility toggle (movies)
gridtog = 'on';						%	grid toggle
derivwarn = 'true';					%	derivative computation warnings display toggle
nTransportSS = 25;					%	number of transport solution of snapshots
nfourier = 1000;					%	number of terms in Fourier expansion (for analytical comparison)
iscomparison2D = 0;					%	whether comparing with 2D solution
iscomparison1D = 0;					%	whether comparing with 1D solution
isexpected = 0;						%	TBD
zoomin = 0;							%	zoom into the front location (movies)
plotmaxfrontR = 0;					%	plot maximum and front values of R using marter
plotfrontwidth = 0;					%	shade reaction zone width around the reaction front
timepausevid = 1.0E-5;				%	time for pausing in movie generation
obtainanalytical = 1;				%	obtain approximate analytical solutions
fracTransitn = 0.5;					%	factor for defining transition zone
forcesolnsimil = 0;					%	forcibly obtain similarity solution even if present from previous postprocessing
forceanalytical = 1;				%	forcibly obtain analytical front solution even if present from previous postprocessing
forcenumfront = 0;					%	forcibly obtain numerical front solution even if present from previous postprocessing
forcemassC = 0;						%	forcibly obtain mass of product even if present from previous postprocessing or from 
									%	during-simulation solution
maxxrange_log = 1E-10;				%	maximum range of x-axis for log scale (transport solution movie)
% similarity solution (approximate analytical solution) computation parameters
simil.nx = 1001;					%	number of grid points on x-axis
simil.absmaxx = 1E3;				%	maximum admissible value of x
simil.power = 2.0;					%	x-grid skewness power
simil.rootfindrelax = 0.75;			%	relaxation in root-finding iterations
simil.tolrootfind = 1E-5;			%	tolerance for convergence of root-finding iterations
simil.barQthres = 10;				%	threshold of barQ for considering it to be large (dispersion-free scenario)
simil.leftBCtuned = 1;				%	whether to tune left-end BC of approximate solution of phi_A (conventional in De Wit)
simil.rootfindmeth = 'NR';			%	root-finding method - NR or STL
simil.maxrootfinditer = 200;		%	maximum number of admissible iterations for root-finding

% overriding front computation as per requirement
iFront =							iFront || iTermsGDE;
speciesname =						['A','B','C'];

% simulation output data structure
namegeom =							{'Line','Point'};
eta =								[0.0E+0,	1.00E-1,	1.00E+0,	1.00E+1,	1.00E+2	];
etaname =							{'No',		'Low',		'Med',		'High',		'VHigh'	};
Pe =								[1.0E+0,	1.0E+1,		1.0E+2,		1.0E+3				];
Pename =							{'1E0',		'1E1',		'1E2',		'1E3'				};

folderloc1Dparent =					'B:\Workspace\reactiveporous\outdata\1D\Median';

for igeom = 1
for iPe = 1
for ieta = length(eta)

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
casename =							['DeWit_2020_',namegeom{igeom},'Inj_',etaname{ieta},		...
									'Dispers_1gamm_1Q_',Pename{iPe},'Pe_',Pename{iPe},'Da'];
% location of specific case data
folderloc1D =						[folderloc1Dparent,'/',casename];

% loading comparison solutions
if (iscomparison1D == 1)
	folderloccompare1Dparent =		'../../../../Validation/1D_Code [Legacy_A]';
	casenamecompare1D =				['Radial_Injection',DomainSizeName{iDomainSize},'_Rescaled_Flow_Pe',Pename{iPe}];
	folderloccompare1D =			[folderloccompare1Dparent,'/',casenamecompare1D];
end
if (iscomparison2D == 1)
	folderloccompare2Dparent =		'../../../../Validation/2D_Code [Legacy_A]';
	casenamecompare2D =				['Radial_Injection',DomainSizeName{iDomainSize},'_Pe',Pename{iPe}];
	folderloccompare2D =			[folderloccompare2Dparent,'/',casenamecompare2D];
end

% analyzing data
Analyzer_1D;

% plotting flow solution
if (iFlow == 1)
	plotFlow;
end

% movie of transport solution
if ((iTransport == 1) || (iTransportSS == 1))
	movTransport;
end

% plotting terms of theta equation
if (iTransportAmBEq == 1)
	movTransportAmB;
end

% plotting mass of product C
if (iMassC == 1)
	plotMassC;
end

% movie of PDF of concentration of C (phi_C)
if (iPDFC == 1)
	movPDFC;
end

% plotting front properties
if (iFront == 1)
	plotFront;
end

% plotting terms of GDE
if (iTermsGDEmovie == 1)
	plotTermsGDE;
end

% movie of terms of GDE
if (iTermsGDEmovie == 1)
	movTermsGDE;
end

% saving mass of C computation
% if ((iMassC == 1) || (iPDFC == 1))
% 	save([folderloc1D,'/soln_Mc.mat'],'soln_analyze');
% end

% catch
% end

end
end
end