% initial clean-up and set-up
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');

% whether to additionally post-process front solution
iFront = 1;							%	plot front properties
% post-processing option toggles
vistog = 'on';						%	plot windows visibility toggle
vistogmov = 'on';					%	plot windows visibility toggle (movies)
gridtog = 'on';						%	grid toggle
derivwarn = 'true';					%	derivative computation warnings display toggle
nfourier = 1000;					%	number of terms in Fourier expansion (for analytical comparison)
iscomparison2D = 0;					%	whether comparing with 2D solution
iscomparison1D = 0;					%	whether comparing with 1D solution
zoomin = 0;							%	zoom into the front location (movies)
plotmaxfrontR = 0;					%	plot maximum and front values of R using marter
plotfrontwidth = 0;					%	shade reaction zone width around the reaction front
timepausevid = 1.0E-5;				%	time for pausing in movie generation
fracTransitn = 0.5;					%	factor for defining transition zone
forcenumfront = 1;					%	should the front properties be numerically computed even if post-proc solution is saved
forceMassCnum =	1;					%	should the mass of product be numerically computed even if post-proc solution is saved
obtainanalytical = 0;				%	whether to obtain analytical solution
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
speciesname = ['A','B','C'];

% simulation output data structure
namegeom =						{'CI'};
eta =							[0.0E+0,		1.00E-3,		1.00E-2,		1.00E-1,		...
								1.00E+0,		1.00E+1,		1.00E+2,		1.00E+3			];
etaname =						{'0d0EP0',		'1d00EM3',		'1d00EM2',		'1d00EM1',		...
								'1d00EP0',		'1d00EP1',		'1d00EP2',		'1d00EP3'		};
Pe =							1.00E+1;
Pename =						{'1d00EP2'};
folderloc1Dparent =				'P:\LocalOnlyRawData_Level_0\1D\1D_Production_B';

namenormalize =					{'','_Normalized','_Ratio'};

for iGeom = 1:1

casename =						[namegeom{iGeom},'_Pe_',Pename{1},'_eta_',etaname{1}];
% location of specific case data
folderloc1D =					[folderloc1Dparent,'/',casename];


for iPe = 1
	
for ieta = 1:length(eta)

% try to clear previous solution left-overs
try
	clear soln1D;
	close soln_analyze;
	clear front;
catch
end

% name of specific case
casename =						[namegeom{iGeom},'_Pe_',Pename{iPe},'_eta_',etaname{ieta}];
% location of specific case data
folderloc1D =					[folderloc1Dparent,'/',casename];

% analyzing data
% loading solution and pre-setting differentiation directions
if (isfile([folderloc1D,'/soln_posproc.mat']) ~= 0)
	load([folderloc1D,'/soln_posproc.mat']);
else
	soln1D =								load([folderloc1D,'/soln.mat']);
	soln1D.simul.compMassC =				0;
	soln1D.simul.compfront =				0;
	if (soln1D.simul.compMassC == 1)
		soln1D.solnreactive =				load([folderloc1D,'/solnreactive.mat']);	
		soln1D.solnreactive =				soln1D.solnreactive.solnreactive;
	end
end
% trimming step save array to unique values
soln1D.simul.tstepsave =					unique(soln1D.simul.tstepsave);
if (obtainanalytical == 1)
	if ((isfield(soln1D,'solnsimil') == 0) || (forcesolnsimil == 1))
		soln1D.simul.simil =				simil;
		if (soln1D.geomdom.isradial == 0)
			soln1D.soln.flow.vel.Px =		soln1D.soln.flow.vel.v;
			soln1D.fltr.flowbc.v =			soln1D.soln.flow.vel.Px(1);
			soln1D.solnsimil =				solver_S0(soln1D.simul,soln1D.fltr);
		elseif (soln1D.geomdom.isradial == 1)
			if (soln1D.fltr.eta == 0)
				soln1D.solnsimil =			solver_S1(soln1D.simul,soln1D.fltr);
				soln1D.solnsimil_disp =		soln1D.solnsimil;
			end
			if ((soln1D.fltr.eta ~= 0) && ((soln1D.fltr.nondim.Pe <= 2.0) && (soln1D.fltr.nondim.Pe >= 0.5)))
				soln1D.solnsimil_disp =		solver_S2(soln1D.simul,soln1D.fltr);
				soln1D.solnsimil =			soln1D.solnsimil_disp;
			end
		else
			soln1D.solnsimil_stat =			solver_S3(soln1D.simul,soln1D.fltr);
			soln1D.solnsimil =				solver_S4(soln1D.simul,soln1D.fltr);
		end
	end
end
dirx    =									cell(1,soln1D.simul.nx);
dirx(:) =									{'cd'};
dirx(1) =									{'fd'};
dirx(soln1D.simul.nx) =						{'bd'};
dirt    =									cell(1,soln1D.simul.nt);
dirt(:) =									{'cd'};
dirt(1) =									{'fd'};
dirt(soln1D.simul.nt) =						{'bd'};
% resetting grid to pre-conformed space if solved in conformed space
try
	if (soln1D.simul.confmap.space == 1)
		soln1D.grd.x =						soln1D.grd.x./(1-soln1D.grd.x);
	end
	if (soln1D.simul.confmap.time == 1)
		soln1D.grd.t =						soln1D.grd.t./(1-soln1D.grd.t);
	end
catch
end
% obtaining time grid that is saved to files
front.tsaved =								soln1D.grd.t(soln1D.simul.tstepsave);
% obtaining advective front location
front.barq_adv =							double(soln1D.geomdom.isradial~=0)*	...
											(double(soln1D.geomdom.isradial+1)^(1/double(soln1D.geomdom.isradial+1)));
front.r_adv =								front.barq_adv*soln1D.grd.t.^(1/double(soln1D.geomdom.isradial+1));
% obtaining numerical reaction rate, reaction zone width
if (iFront == 1)
	if ((exist('front','var') == 0) || (forcenumfront == 1))
		neednumfrontsoln =						1;
	else
		if (isfield(soln1D.front,'num') == 0)
			neednumfrontsoln =					1;
		else
			neednumfrontsoln =					0;
		end
	end
	if (neednumfrontsoln == 1)
		if ((soln1D.simul.compfront == 0) || ((forcenumfront == 1) && (soln1D.simul.compfront == 1)))
			ifront.tsaved =						0;
			for it = soln1D.simul.tstepsave
				ifront.tsaved =			ifront.tsaved+1;
				phiCcurrent =			load([folderloc1D,'/phiC_',num2str(it),'.dat']);
				phiAcurrent =			load([folderloc1D,'/phiA_',num2str(it),'.dat']);
				phiBcurrent =			load([folderloc1D,'/phiB_',num2str(it),'.dat']);
				front.num.R =			phiAcurrent.*phiBcurrent;
				front.num.theta =		phiAcurrent-phiBcurrent;
				[front.num.R_f(ifront.tsaved),front.num.inodrf(ifront.tsaved)] =		...
										max(front.num.R);
				[~,front.num.inodrf_orig(ifront.tsaved)] =								...
										min(front.num.theta.^2);
				front.num.r_f(ifront.tsaved)=											...
										soln1D.grd.x(front.num.inodrf(ifront.tsaved));
				front.num.r_f_orig(ifront.tsaved) =			...
										soln1D.grd.x(front.num.inodrf_orig(ifront.tsaved));
				front.num.R_f_orig(ifront.tsaved) =			...
										front.num.R(front.num.inodrf_orig(ifront.tsaved));
				for iwidthdef = 1:2
					IntNr =					trapz(soln1D.grd.x,		...
											(soln1D.grd.x-front.num.r_f_orig(ifront.tsaved)).^2.*front.num.R.*	...
											(soln1D.grd.x.^(double((iwidthdef-1)*soln1D.geomdom.isradial))));
					IntDr =					trapz(soln1D.grd.x,front.num.R.*	...
											(soln1D.grd.x.^(double((iwidthdef-1)*soln1D.geomdom.isradial))));
					front.num.w_f_orig(ifront.tsaved,iwidthdef) =			...
											sqrt(IntNr/IntDr);
				end
				leftswitch =			1;
				rightswitch =			1;
				ixleft =				1;
				ixright =				soln1D.simul.nx;
				for ix = 1:soln1D.simul.nx
					if ((front.num.R(ix) > 0.5*(front.num.R_f(ifront.tsaved))) && (leftswitch == 1))
						ixleft =		ix;
						leftswitch =	0;
					end
					if ((front.num.R(ix) < 0.5*(front.num.R_f(ifront.tsaved))) && (leftswitch == 0) && (rightswitch == 1))
						ixright =		ix-1;
						rightswitch =	0;
					end
				end
				front.num.w_f(ifront.tsaved) =								...
										soln1D.grd.x(ixright)-soln1D.grd.x(ixleft);	
			end
		end
	end
end
% obtaining transition zones
front.trnszon =									ones(2,length(fracTransitn),2,2,2);
if (soln1D.simul.compfront == 1)
	for ifracTransitn = 1:length(fracTransitn)
		[~,front.trnszon(1,ifracTransitn,1,1,1)]=	...
											min((soln1D.solnreactive.xforig-	...
											((fracTransitn(ifracTransitn)*soln1D.fltr.eta)^		...
											(1/double(soln1D.geomdom.isradial)))).^2);
		[~,front.trnszon(1,ifracTransitn,1,2,1)]=	...
											min((soln1D.solnreactive.xforig-	...
											(((1/fracTransitn(ifracTransitn))*soln1D.fltr.eta)^	...
											(1/double(soln1D.geomdom.isradial)))).^2);
		[~,front.trnszon(1,ifracTransitn,2,1,1)]=	...
											min((soln1D.solnreactive.xforig-	...
											((fracTransitn(ifracTransitn)*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^		...
											(1/double(soln1D.geomdom.isradial)))).^2);
		[~,front.trnszon(1,ifracTransitn,2,2,1)]=	...
											min((soln1D.solnreactive.xforig-	...
											(((1/fracTransitn(ifracTransitn))*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^	...
											(1/double(soln1D.geomdom.isradial)))).^2);
	end
	for ilimdef = 1:2
		for iside = 1:2
			[~,front.trnszon(1,ifracTransitn,ilimdef,iside,2)] =	...
											min((front.tsaved-soln1D.grd.t(front.trnszon(1,ifracTransitn,ilimdef,iside,1))).^2);
		end
	end
end
if (iFront == 1)
	if ((soln1D.simul.compfront == 0) || ((forcenumfront == 1) && (soln1D.simul.compfront == 1)))
		for ifracTransitn = 1:length(fracTransitn)
			[~,front.trnszon(2,ifracTransitn,1,1,2)] =			...
											min((front.num.r_f-	...
											((fracTransitn(ifracTransitn)*soln1D.fltr.eta)^								...
											(1/double(soln1D.geomdom.isradial)))).^2);
			[~,front.trnszon(2,ifracTransitn,1,2,2)] =		...
											min((front.num.r_f-	...
											(((1/fracTransitn(ifracTransitn))*soln1D.fltr.eta)^							...
											(1/double(soln1D.geomdom.isradial)))).^2);
			[~,front.trnszon(2,ifracTransitn,2,1,2)] =		...
											min((front.num.r_f-	...
											((fracTransitn(ifracTransitn)*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^		...
											(1/double(soln1D.geomdom.isradial)))).^2);
			[~,front.trnszon(2,ifracTransitn,2,2,2)] =		...
											min((front.num.r_f-	...
											(((1/fracTransitn(ifracTransitn))*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^	...
											(1/double(soln1D.geomdom.isradial)))).^2);
			for ilimdef = 1:2
				for iside = 1:2
					front.trnszon(2,ifracTransitn,ilimdef,iside,1) =	...
											soln1D.simul.tstepsave(front.trnszon(2,ifracTransitn,ilimdef,iside,2));
				end
			end
		end
	end
end

% obtaining mass of product C
if ((exist('soln_analyze','var') == 0) || (forceMassCnum == 1))
	needmassCsoln =							1;
else
	needmassCsoln =							0;
end
if (needmassCsoln == 1)
	soln_analyze.MassC =					zeros(length(soln1D.simul.tstepsave),1);
	soln_analyze.MassC_4mR =				zeros(length(soln1D.simul.tstepsave),1);
	soln_analyze.MassC_4mR_twk=				zeros(length(soln1D.simul.tstepsave),1);
	soln_analyze.barR =						zeros(length(soln1D.simul.tstepsave),1);
	ifront.tsaved = 0;
	for it = soln1D.simul.tstepsave
		ifront.tsaved =						ifront.tsaved+1;
		phiCcurrent =						load([folderloc1D,'/phiC_',num2str(it),'.dat']);
		phiBcurrent =						load([folderloc1D,'/phiB_',num2str(it),'.dat']);
		phiAcurrent =						load([folderloc1D,'/phiA_',num2str(it),'.dat']);
		soln_analyze.MassC(ifront.tsaved) =	((double(soln1D.geomdom.isradial)*2*pi)^double(soln1D.geomdom.isradial~=0))*		...
											trapz(soln1D.grd.x,phiCcurrent.*(soln1D.grd.x.^soln1D.geomdom.isradial));
		soln_analyze.barR(ifront.tsaved) =	((double(soln1D.geomdom.isradial)*2*pi)^double(soln1D.geomdom.isradial~=0))*		...
											trapz(soln1D.grd.x,phiAcurrent.*phiBcurrent.*(soln1D.grd.x.^soln1D.geomdom.isradial));
		if (ifront.tsaved > 1)
			soln_analyze.MassC_4mR(ifront.tsaved) =		...
											0.5*(soln_analyze.barR(ifront.tsaved)+soln_analyze.barR(ifront.tsaved-1))*	...
											(soln1D.grd.t(ifront.tsaved)-soln1D.grd.t(ifront.tsaved));
			soln_analyze.MassC_4mR_twk(ifront.tsaved) =		...
											soln_analyze.MassC_4mR(ifront.tsaved) +			...
											2*pi*phiCcurrent(1)*(soln1D.grd.t(ifront.tsaved)-soln1D.grd.t(ifront.tsaved));
		end
	end
end

% generating comparison data
if (ieta == 1)
	soln1D_base =		soln1D;
	soln_analyze_base =	soln_analyze;
	front_base =		front;
else
	if (sum((front.tsaved-front_base.tsaved).^2) ~= 0)
		soln_analyze.MassC =			pchip(front.tsaved,soln_analyze.MassC,front_base.tsaved);
		soln_analyze.MassC_4mR =		pchip(front.tsaved,soln_analyze.MassC_4mR,front_base.tsaved);
		soln_analyze.MassC_4mR_twk =	pchip(front.tsaved,soln_analyze.MassC_4mR_twk,front_base.tsaved);
		soln_analyze.barR =				pchip(front.tsaved,soln_analyze.barR ,front_base.tsaved);
		front.tsaved =					front_base.tsaved;
	end
	if ((soln1D.simul.compMassC == 1) && (soln1D_base.simul.compMassC == 1))
		if (sum((soln1D.grd.t-soln1D_base.grd.t).^2) ~= 0)
			soln1D.solnreactive.MassC =	pchip(soln1D.grd.t,soln1D.solnreactive.MassC,soln1D_base.grd.t);
			soln1D.solnreactive.MassC_4mR =			...
										pchip(soln1D.grd.t,soln1D.solnreactive.MassC_4mR,soln1D_base.grd.t);
			soln1D.solnreactive.MassC_4mR_twk =		...
										pchip(soln1D.grd.t,soln1D.solnreactive.MassC_4mR_twk,soln1D_base.grd.t);
			soln1D.solnreactive.barR =	pchip(soln1D.grd.t,soln1D.solnreactive.barR,soln1D_base.grd.t);
			soln1D.solnreactive.dMassCdt =			...
										pchip(soln1D.grd.t,soln1D.solnreactive.dMassCdt,soln1D_base.grd.t);
			soln1D.grd.t =				soln1D_base.grd.t;
		end
	end
end
for inormalize = 0:2
	if (inormalize == 0)
		soln_analyze.MassC =				soln_analyze.MassC-soln_analyze_base.MassC;
		soln_analyze.MassC_4mR =			soln_analyze.MassC_4mR-soln_analyze_base.MassC_4mR;
		soln_analyze.MassC_4mR =			soln_analyze.MassC_4mR_twk-soln_analyze_base.MassC_4mR_twk;
		soln_analyze.barR =					soln_analyze.barR-soln_analyze_base.barR;
		if ((soln1D.simul.compMassC == 1) && (soln1D_base.simul.compMassC == 1))
			soln1D.solnreactive.MassC =		soln1D.solnreactive.MassC-soln1D_base.solnreactive.MassC;
			soln1D.solnreactive.MassC_4mR =	soln1D.solnreactive.MassC_4mR-soln1D_base.solnreactive.MassC_4mR;
			soln1D.solnreactive.MassC_4mR_twk =			...
											soln1D.solnreactive.MassC_4mR_twk-soln1D_base.solnreactive.MassC_4mR_twk;
			soln1D.solnreactive.barR =		soln1D.solnreactive.barR-soln1D_base.solnreactive.barR;
			soln1D.solnreactive.dMassCdt =	soln1D.solnreactive.dMassCdt-soln1D_base.solnreactive.dMassCdt;
		end
	elseif (inormalize == 1)
		soln_analyze.MassC =				soln_analyze.MassC./soln_analyze_base.MassC-1;
		soln_analyze.MassC_4mR =			soln_analyze.MassC_4mR./soln_analyze_base.MassC_4mR-1;
		soln_analyze.MassC_4mR =			soln_analyze.MassC_4mR_twk./soln_analyze_base.MassC_4mR_twk-1;
		soln_analyze.barR =					soln_analyze.barR./soln_analyze_base.barR-1;
		if ((soln1D.simul.compMassC == 1) && (soln1D_base.simul.compMassC == 1))
			soln1D.solnreactive.MassC =		soln1D.solnreactive.MassC./soln1D_base.solnreactive.MassC-1;
			soln1D.solnreactive.MassC_4mR =	soln1D.solnreactive.MassC_4mR./soln1D_base.solnreactive.MassC_4mR-1;
			soln1D.solnreactive.MassC_4mR_twk =			...
											soln1D.solnreactive.MassC_4mR_twk./soln1D_base.solnreactive.MassC_4mR_twk-1;
			soln1D.solnreactive.barR =		soln1D.solnreactive.barR./soln1D_base.solnreactive.barR-1;
			soln1D.solnreactive.dMassCdt =	soln1D.solnreactive.dMassCdt./soln1D_base.solnreactive.dMassCdt-1;
		end
	else
		soln_analyze.MassC =				soln_analyze.MassC./soln_analyze_base.MassC;
		soln_analyze.MassC_4mR =			soln_analyze.MassC_4mR./soln_analyze_base.MassC_4mR;
		soln_analyze.MassC_4mR =			soln_analyze.MassC_4mR_twk./soln_analyze_base.MassC_4mR_twk;
		soln_analyze.barR =					soln_analyze.barR./soln_analyze_base.barR-1;
		if ((soln1D.simul.compMassC == 1) && (soln1D_base.simul.compMassC == 1))
			soln1D.solnreactive.MassC =		soln1D.solnreactive.MassC./soln1D_base.solnreactive.MassC;
			soln1D.solnreactive.MassC_4mR =	soln1D.solnreactive.MassC_4mR./soln1D_base.solnreactive.MassC_4mR;
			soln1D.solnreactive.MassC_4mR_twk =			...
											soln1D.solnreactive.MassC_4mR_twk./soln1D_base.solnreactive.MassC_4mR_twk;
			soln1D.solnreactive.barR =		soln1D.solnreactive.barR./soln1D_base.solnreactive.barR;
			soln1D.solnreactive.dMassCdt =	soln1D.solnreactive.dMassCdt./soln1D_base.solnreactive.dMassCdt;
	end
end

% obtaining comparison plots
differenceMassC;

end

disp(casename);

end
end
end