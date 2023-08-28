%-----------------------------------------------------------------------------------------------------------------------------------
% Wrapper [Variant - Production] for
% 1D Stationary/Steady Reactive Transport Solver for Spherical Injection using Finite Difference Method (FDM)
% [Non-Function Program]
%
%	written by --
%	uddipta ghosh (uddipta.ghosh@iitgn.ac.in)
%	pratyaksh karan (pratyakshkaran@gmail.com)
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% initiating cleaned matlab sessions
clc;
clear;
close all;
fclose all;
warning('off','MATLAB:nearlySingularMatrix');
addpath('../functions');
addpath('../../postproc/anasolvers');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% solver switches
solvetransport =				1;
plotstatfront =					0;
plotwidthfront =				0;
plotbarR =						0;
solvesimil =					1;
takeguessfromsimil =			0;
takeguesszero =					0;
takeguessone =					0;
takeguessprev =					0;
nplotnormalized =				2;
vistog =						'on';
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% definings case combinations
namegeom =						{'SI'};
igeom =							1;
eta =							1.00E+3;
etaname =						{'1d00EP3'};
Pe =							1.00E+0;
Pename =						{'1d00EP0'};
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% obtaining grid skewness
Pe_b =							[1.0E+0,	1.00E+3];
eta_b =							[1.00E-3,	1.00E+3];
[Pe_mg_b, eta_mg_b] =			meshgrid(Pe_b,eta_b);
[Pe_mg, eta_mg] =				meshgrid(Pe,eta(2:end));
xgridskew_CI_b =				[3.30, 2.60; 2.62, 2.59];
xgridskew_SI_b =				[3.72, 4.20; 2.59, 3.05];	
xgridskew_CI_b_noDisp =			[3.27, 2.62];
xgridskew_SI_b_noDisp =			[3.70, 3.11];
xgridskew_CI =					interp2(Pe_mg_b,eta_mg_b,xgridskew_CI_b,Pe_mg,eta_mg,'spline');
xgridskew_SI =					interp2(Pe_mg_b,eta_mg_b,xgridskew_SI_b,Pe_mg,eta_mg,'spline');
xgridskew_CI_noDisp =			interp1(Pe_b,xgridskew_CI_b_noDisp,Pe,'pchip');
xgridskew_SI_noDisp =			interp1(Pe_b,xgridskew_SI_b_noDisp,Pe,'pchip');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% Setting Output Folder
simul.outfoldrname =			'../../outdata/sphstat';
mkdir(simul.outfoldrname);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% Input Parameters
simul.localonly =				1;					%	save output in a local-only folder
simul.localonlylevel =			1;					%	local-only folder location to save output data
%-----------------------------------------------------------------------------------------------------------------------------------
% Simulation Parameters
% Simulation Solver Parameters
simul.solveGPU = 				0;					%	whether to solve linear system on GPU (1 = yes, 0 = no)
simul.nx =						50001;				%	number of grid points in the x-coordinate
simul.xgridskew =				2.0;				%	skewing parameter for xgrid [SI -> Pe end eta end]
													%	(higher than unity index - denser near start)
													%	(lower than unity index - denser near end)
simul.STLtol =					1.00E-09;			%	tolerance of update diminution for Source Term Linearization Iterations
simul.STLrelax =				0.90;				%	relaxation of solution update
simul.errSTLmode =				1;					%	mode of computing STL error
													%	1 - standard deviation from last iteration solution
													%	2 -	standard deviation from last iteration solution times number of pts
													%	3 - maximum difference from last iteration
simul.STexplicit =				0;					%	whether to explicitly solve reaction term (i.e. the source term ST)
simul.compactgeom =				0;					%	make the geometry compact to optimize computation time
simul.compactgeomfactotal =		25;					%	geometry compaction factor for front location
simul.compactgeomsubfacwidth =	25;					%	geometry compaction factor for front width
simul.limitdx =					1;					%	override specified point count to limit dx step size (1 = yes, 2 = no)
simul.dxLmax =					1.00E-02;			%	maximum allowed dx at left end (enforced when simul.limitdx is 1)
simul.dxRmax =					1.00E+02;			%	maximum allowed dx at right end (enforced when simul.limitdx is 1)
simul.dxRmax_comprtoA =			1.00E-02;			%	maximum allowed dx at right end as compared to A
													%	(enforced when simul.limitdx is 1)
simul.frontparamsrel =			1;					%	are front parameters presented normalized to Lx (1 = yes, 0 = no)
													%	(when yes, the front parameters are multiplied to Lx before proceeding
													%	to simulation)
simul.injectparamsrel =			0;					%	are injection parameters presented normalized to Lx (1 = yes, 0 = no)
% Simulation System Parameters
simul.justinject =				1;					%	whether species A has just been injected at start time
													%	(if set to 0, generate a front, else, start injection at 
													%	injection plane/line/point)
simul.isflowrate =				1;					%	is flowrate given (1) or inlet pressure head (0)
simul.givenvel =				1;					%	0 - solve Darcy solver to get velocity field or 1 - take user-specified 
													%	velocity profile
simul.nodiffusion = 			0;					%	whether to consider diffusion (0) or not (1)
simul.derivwarn =				true;				%	whether to warn before switching derivative direction
% Geometric & Domain Parameters
geomdom.isradial =				2;					%	is the domain radial or linear
geomdom.size.Lx =				1.00E+07;			%	non-dimensionalized x-coordinate length
% domain boundary conditions for transport - species concentration and its gradient -
% 1 - 	given constant net hydraulic head
% 2 - 	given flow rate (i.e. net hydraulic head gradient)
geomdom.flowbc.left =			1;					%	near-end transport boundary conditions
geomdom.flowbc.right =			1;					%	far-end transport boundary conditions
% domain boundary conditions for transport - species concentration and its gradient - same meanings as for flow -
geomdom.transpbc.left =			[1,1,0];			%	near-end transport boundary conditions
geomdom.transpbc.right =		[1,1,0];			%	far-end transport boundary conditions
% domain permeability switch (see Flow & Transport Parameters for input variables) -
% 0 - 	uniform permeability
% 1 - 	specified permeability based on function of x and y-coordinates
% 2 -		stochastic permeability generated based on stochasticticity specifications (using a function call to permgen)
% Flow & Transport Parameters
fltr.nondim.gamm =				1.00;				%	concentration of B at far end
fltr.flowbc.Hleft0 =			0.00;				%	pressure head boundary condition constant - left end
fltr.flowbc.Hright0 =			0.00;				%	pressure head boundary condition constant - far-end
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
%-----------------------------------------------------------------------------------------------------------------------------------
% Similarity Solution Parameters
simul.simil.nx =				2501;				%	no of grid points to take for similarity solution
simul.simil.absmaxx =			1.0E2;				%	size of z-grid for similarity solution
simul.simil.power =				15.0;				%	power for skewing z-grid for similarity solution
simul.simil.tolrootfind =		1.0E-5;				%	root-finding tolerance for similarity solution
simul.simil.rootfindrelax =		0.9;				%	iteration relaxation in root-finding for similarity solution
simul.simil.leftBCtuned =		1;					%	whether to take -Kz as the left end limit for phi_A for simil solution
simul.simil.rootfindmeth =		'NR';				%	root-finding method for similarity solution
simul.simil.maxrootfinditer =	5E4;				%	maximum number of interactions for similarity solution
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% running preprocessor
% computed geometric & domain parameters
geomdom.size.a =				fltr.transpbc.inject.k1;
geomdom.size.A =				geomdom.size.Lx;
% override space step count to satisfy space step size constraints	
if (simul.limitdx == 1)
	simul.dxRmax =				min([simul.dxRmax,simul.dxRmax_comprtoA*geomdom.size.A]);
	simul.nx =					max([simul.nx,	...
								ceil(1+max([	...
								((geomdom.size.A-geomdom.size.a)/simul.dxLmax)^(1/simul.xgridskew),	...
								1/(1.0-((1.0-(simul.dxRmax/(geomdom.size.A-geomdom.size.a)))^(1/simul.xgridskew)))	...
								]))]);
end
% obtaining absolute front and injection parameters when provided values are relative
if (simul.frontparamsrel == 1)
		fltr.init.front.k1 =	fltr.init.front.k1*geomdom.size.Lx;
		fltr.init.front.k2 =	fltr.init.front.k2*geomdom.size.Lx;
end
if (simul.injectparamsrel == 1)
		fltr.transpbc.inject.k1=fltr.transpbc.inject.k1*geomdom.size.Lx;
end
% space grid generation
leftend =						geomdom.size.a;
rightend =						geomdom.size.A;
grd.x =							leftend+(rightend-leftend)*linspace(0.0,1.0,simul.nx).^simul.xgridskew;
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Pre-Allocating Matrices
soln.ana.rfS =								zeros(length(Pe),length(eta));
soln.ana.thetaS =							zeros(simul.nx,length(Pe),length(eta));
soln.ana.rfS_s =							zeros(length(Pe),length(eta));
soln.ana.rfS_w =							zeros(length(Pe),length(eta));
soln.ana.thetaS_s =							zeros(simul.nx,length(Pe),length(eta));
soln.ana.thetaS_w =							zeros(simul.nx,length(Pe),length(eta));
soln.simil.x =								zeros(simul.simil.nx,length(Pe),length(eta));
soln.simil.rfS =							zeros(length(Pe),length(eta));
soln.simil.wfS =							zeros(length(Pe),length(eta),3);
soln.simil.barRS =							zeros(length(Pe),length(eta));
soln.simil.thetaS =							zeros(simul.simil.nx,length(Pe),length(eta));
soln.simil.phiS =							zeros(simul.simil.nx,length(Pe),length(eta),2);
soln.simil.RS =								zeros(simul.simil.nx,length(Pe),length(eta));
soln.num.rfS =								zeros(length(Pe),length(eta));
soln.num.wfS =								zeros(length(Pe),length(eta),2,2);
soln.num.barRS =							zeros(length(Pe),length(eta));
soln.num.thetaS =							zeros(simul.nx,length(Pe),length(eta));
soln.num.phiS =								zeros(simul.nx,length(Pe),length(eta),2);
soln.num.RS =								zeros(simul.nx,length(Pe),length(eta));
for iPe = 1
%	Pre-Allocating Iteration-Specific Matrices
	if ((solvesimil == 1) && (takeguessfromsimil == 1))
		soln.simil.phiS_interp =			zeros(simul.nx,length(eta),2);
	end
	if (solvetransport == 1)
		for ieta = 1
%			setting paramter values for current loop	
			fltr.nondim.Pe =				Pe(iPe);			%	expected Peclet number
			fltr.eta =						eta(ieta);			%	expected longitudinal dispersivity factor
			simul.casename =				etaname{ieta};		%	case name
%			analytical front location
			soln.ana.rfS(iPe,ieta) =		sqrt(fltr.eta*fltr.nondim.Pe)*tan(sqrt(fltr.eta/(fltr.nondim.Pe))*		...
											log(1-(1/(1+fltr.nondim.gamm))*(1-exp((pi/2)*sqrt(fltr.nondim.Pe/fltr.eta)))));
			soln.ana.rfS_w(iPe,ieta) =		fltr.nondim.Pe/log(1+fltr.nondim.gamm);
			soln.ana.rfS_S(iPe,ieta) =		sqrt(fltr.eta*fltr.nondim.Pe);
%			analytical stationary solution for theta
			soln.ana.thetaS(1:simul.nx,iPe,ieta) =	...
											1-(1+fltr.nondim.gamm)*((1-exp(sqrt(fltr.nondim.Pe/fltr.eta)*	...
											atan(grd.x/sqrt(fltr.eta*fltr.nondim.Pe))))/					...
											(1-exp((pi/2)*sqrt(fltr.nondim.Pe/fltr.eta))));
			soln.ana.thetaS_w(1:simul.nx,iPe,ieta)=	...
											1-(1+fltr.nondim.gamm)*exp(-(fltr.nondim.Pe./grd.x));
			soln.ana.thetaS_s(1:simul.nx,iPe,ieta)=	...
											1-(1+fltr.nondim.gamm)*((2*grd.x)/(pi*sqrt(fltr.nondim.Pe)))*(1/sqrt(fltr.eta));
			if (fltr.eta == 0)
				soln.ana.thetaS(1:simul.nx,iPe,ieta) =	...
											soln.ana.thetaS_w(1:simul.nx,iPe,ieta);
				soln.ana.rfS(iPe,ieta) =	soln.ana.rfS_w(iPe,ieta);
			end
%			analytical stationary solution approximation for phiA and phiB
			if (solvesimil == 1)
				solnS3 =					solver_S3(simul,fltr);
				soln.simil.x(1:simul.simil.nx,iPe,ieta) =	...
											solnS3.r(1:simul.simil.nx);
				soln.simil.phiS(1:simul.simil.nx,iPe,ieta,1) =	...
											solnS3.phiAS(1:simul.simil.nx);
				soln.simil.phiS(1:simul.simil.nx,iPe,ieta,2) =	...
											solnS3.phiBS(1:simul.simil.nx);
				soln.simil.thetaS(1:simul.simil.nx,iPe,ieta) =	...
											solnS3.thetaS(1:simul.simil.nx);
				soln.simil.RS(1:simul.simil.nx,iPe,ieta) =		...
											solnS3.RS(1:simul.simil.nx);
				soln.simil.barRS(iPe,ieta)=	4*pi*trapz(solnS3.r,solnS3.r.^2.*solnS3.RS);
				soln.simil.rfS(iPe,ieta) =	solnS3.r_fS;
				soln.simil.wfS(iPe,ieta,1)=	solnS3.w_fS;
				soln.simil.wfS(iPe,ieta,2)=	solnS3.w_fS_alt;
				soln.simil.wfS(ieta,iPe,3)=	solnS3.w_fS_approx;
			end
			if ((solvesimil == 1) && (takeguessfromsimil == 1))
				soln.simil.phiS_interp(1:simul.nx,ieta,1) =		...
											pchip(solnS3.r,solnS3.phiAS,grd.x);
				soln.simil.phiS_interp(1:simul.nx,ieta,2) =		...
											pchip(solnS3.r,solnS3.phiBS,grd.x);
				fltr.init.phiA =			soln.simil.phiS_interp(1:simul.nx,ieta,1);
				fltr.init.phiB =			soln.simil.phiS_interp(1:simul.nx,ieta,2);
			elseif (takeguesszero == 1)
				fltr.init.phiA =			zeros(1,simul.nx);
				fltr.init.phiB =			zeros(1,simul.nx);
			elseif (takeguessone == 1)
				fltr.init.phiA =			ones(1,simul.nx);
				fltr.init.phiB =			ones(1,simul.nx);
			elseif (takeguessprev == 1)
				try
					fltr.init.phiA =		soln.num.phiS(1:simul.nx,iPe-1,ieta-1,1);
					fltr.init.phiB =		soln.num.phiS(1:simul.nx,iPe-1,ieta-1,2);
				catch
					try
						fltr.init.phiA =	soln.num.phiS(1:simul.nx,iPe-1,ieta,1);
						fltr.init.phiB =	soln.num.phiS(1:simul.nx,iPe-1,ieta,2);
					catch
						try
							fltr.init.phiA=	soln.num.phiS(1:simul.nx,iPe,ieta-1,1);
							fltr.init.phiB=	soln.num.phiS(1:simul.nx,iPe,ieta-1,2);
						catch
							[~,ind_fS] =	min((grd.x-soln.ana.rfS(ieta)).^2);
							fltr.init.phiA=	0.5*(1.0-tanh((grd.x-grd.x(ind_fS))/fltr.init.front.k2));
							fltr.init.phiB=	fltr.nondim.gamm*0.5*(1.0+tanh((grd.x-grd.x(ind_fS))/fltr.init.front.k2));
						end
					end
				end
			else
				[~,ind_fS] =				min((grd.x-soln.ana.rfS(ieta)).^2);
				fltr.init.phiA =			0.5*(1.0-tanh((grd.x-grd.x(ind_fS))/fltr.init.front.k2));
				fltr.init.phiB =			fltr.nondim.gamm*0.5*(1.0+tanh((grd.x-grd.x(ind_fS))/fltr.init.front.k2));
			end
%			numerical stationary solution
%			obtaining flow solution
			if (simul.givenvel == 0)
				[soln.flow.H, soln.flow.vel, fltr.perm.Kperm] =		...
											solver1D_darcy(geomdom,simul,fltr,grd);
			else
				fltr.perm.Kperm =			fltr.perm.Kperm0*ones(simul.nx,1);
				soln.flow.H =				zeros(simul.nx,1);
				soln.flow.vel.v =			double(geomdom.isradial~=0)*(1./(grd.x.^double(geomdom.isradial)));
			end
%			obtaining transport solution
			phitemp =						stationary1D_reactive(geomdom,simul,fltr,grd,soln.flow);
			for ispecies = 1:2
				soln.num.phiS(1:simul.nx,iPe,ieta,ispecies) =	...
											phitemp(1:simul.nx,ispecies);
			end
			soln.num.thetaS(1:simul.nx,iPe,ieta) =				...
											soln.num.phiS(1:simul.nx,iPe,ieta,1)-soln.num.phiS(1:simul.nx,iPe,ieta,2);
			soln.num.RS(1:simul.nx,iPe,ieta) =					...
											soln.num.phiS(1:simul.nx,iPe,ieta,1).*soln.num.phiS(1:simul.nx,iPe,ieta,2);
			soln.num.barRS(iPe,ieta) =		4*pi*trapz(grd.x,grd.x.^2.*phitemp(1:simul.nx,1)'.*phitemp(1:simul.nx,2)');
			[~,ind_fS] =					min(soln.num.thetaS(1:simul.nx,iPe,ieta).^2);
			soln.num.rfS(iPe,ieta) =		grd.x(ind_fS);
			IntNr =							trapz(grd.x,grd.x.^2.*(grd.x-soln.num.rfS(ieta)).^2.* 		...
											permute(soln.num.RS(1:simul.nx,iPe,ieta),[2,1]));
			IntDr =							trapz(grd.x,grd.x.^2.*permute(soln.num.RS(1:simul.nx,iPe,ieta),[2,1]));
			IntNr_alt =						trapz(grd.x,(grd.x-soln.num.rfS(ieta)).^2.*					...
											permute(soln.num.RS(1:simul.nx,iPe,ieta),[2,1]));
			IntDr_alt =						trapz(grd.x,permute(soln.num.RS(1:simul.nx,iPe,ieta),[2,1]));
			IntNr_theta =					trapz(grd.x,grd.x.^2.*(grd.x-soln.num.rfS(iPe,ieta)).^2.*	...
											((double(grd.x<soln.num.rfS(iPe,ieta))*				...
											(1-fltr.nondim.gamm.^2)+fltr.nondim.gamm.^2)-				...
											(permute(soln.num.thetaS(1:simul.nx,iPe,ieta),[2,1])).^2));
			IntDr_theta =					trapz(grd.x,grd.x.^2.*(1-permute(soln.num.thetaS(1:simul.nx,iPe,ieta),[2,1]).^2));
			IntNr_theta_alt =				trapz(grd.x,(grd.x-soln.num.rfS(ieta)).^2.*					...
											((double(grd.x<soln.num.rfS(iPe,ieta))*				...
											(1-fltr.nondim.gamm.^2)+fltr.nondim.gamm.^2)-				...
											(permute(soln.num.thetaS(1:simul.nx,iPe,ieta),[2,1])).^2));
			IntDr_theta_alt =				trapz(grd.x,(1-permute(soln.num.thetaS(1:simul.nx,iPe,ieta),[2,1]).^2));
			soln.num.wfS(iPe,ieta,1,1) =	sqrt(IntNr/IntDr);
			soln.num.wfS(iPe,ieta,2,1) =	sqrt(IntNr_alt/IntDr_alt);
			soln.num.wfS(iPe,ieta,1,2) =	sqrt(IntNr_theta/IntDr_theta);
			soln.num.wfS(iPe,ieta,2,2) =	sqrt(IntNr_theta_alt/IntDr_theta_alt);
%			plotting species solutions
			variantname =	{'','bysqrteta','bysqrtetaPe'};
			for ivariant = 1:1+nplotnormalized
				xnormfac =	sqrt((fltr.eta*(fltr.nondim.Pe^double(ivariant==3)))^double(ivariant~=1));
				figl  =	figure('position',[100,100,1250,500],'visible',vistog); axl = axes(figl,'position',[0.125,0.2,0.80,0.7]);
				set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex');
				set(axl,'xscale','log');
				set(gca,'xlim',[geomdom.size.a,geomdom.size.A]);
				hold on;
				plot(axl,grd.x/xnormfac,soln.ana.thetaS(1:simul.nx,iPe,ieta),'k-','linewidth',2);
				hold on;
				plot(axl,grd.x/xnormfac,soln.num.thetaS(1:simul.nx,iPe,ieta),'k--','linewidth',2);
				plot(axl,grd.x/xnormfac,soln.num.phiS(1:simul.nx,iPe,ieta,1),'g--','linewidth',2);
				plot(axl,grd.x/xnormfac,soln.num.phiS(1:simul.nx,iPe,ieta,2),'b--','linewidth',2);
				plot(axl,grd.x/xnormfac,soln.num.RS(1:simul.nx,iPe,ieta)/max(soln.num.RS(1:simul.nx,iPe,ieta)),'m--','linewidth',2);
				ylimin =	get(axl,'ylim');
				xlimin =	get(axl,'xlim'); xliml = xlimin(1); xlimr = xlimin(1);
				plot(axl,grd.x/xnormfac,soln.ana.thetaS_s(:,iPe,ieta),':','color',[0.1,0.4,0.5],'linewidth',1);
				plot(axl,grd.x/xnormfac,soln.ana.thetaS_w(:,iPe,ieta),'--','color',[0.1,0.4,0.5],'linewidth',1);
				if(ylimin(1)<0)
					ylimin(1) =		ylimin(1)*1.25;
				else
					ylimin(1) =		ylimin(1)*0.75;
				end
				if(ylimin(2)<0)
					ylimin(2) =		ylimin(2)*0.75;
				else
					ylimin(2) =		ylimin(2)*1.25;
				end
				if (solvesimil == 1)
					plot(axl,solnS3.r/xnormfac,solnS3.thetaS,'k:','linewidth',2);
					plot(axl,solnS3.r/xnormfac,solnS3.phiAS,'g:','linewidth',2);
					plot(axl,solnS3.r/xnormfac,solnS3.phiBS,'b:','linewidth',2);
					plot(axl,solnS3.r/xnormfac,solnS3.RS/max(soln.num.RS(1:simul.nx,ieta)),'m:','linewidth',2);
				end
				if ((solvesimil == 1) && (takeguessfromsimil == 1))
					plot(axl,grd.x/xnormfac,soln.simil.phiS_interp(1:simul.nx,ieta,1),'g:','linewidth',1.5);
					plot(axl,grd.x/xnormfac,soln.simil.phiS_interp(1:simul.nx,ieta,2),'b:','linewidth',1.5);
				end
				plot(axl,(soln.ana.rfS(ieta)*ones(1,simul.nx))/xnormfac,linspace(ylimin(1),ylimin(2),simul.nx),			...
				'-','linewidth',0.5,'color',[0.3,0.3,0.3]);
				plot(axl,(soln.num.rfS(ieta)*ones(1,simul.nx))/xnormfac,linspace(ylimin(1),ylimin(2),simul.nx),			...
				'--','linewidth',0.5,'color',[0.3,0.3,0.3]);
				if (solvesimil == 1)
					plot(axl,(solnS3.r_fS*ones(1,simul.nx))/xnormfac,linspace(ylimin(1),ylimin(2),simul.nx),			...
					'--','linewidth',0.5,'color',[0.35,0.3,0]);
				end
				plot(axl,(sqrt(fltr.eta)*ones(1,simul.nx))/xnormfac,linspace(ylimin(1),ylimin(2),simul.nx),					...
				':','linewidth',0.5,'color',[0.5,0.5,0]);
				plot(axl,(sqrt(fltr.eta*fltr.nondim.Pe)*ones(1,simul.nx))/xnormfac,linspace(ylimin(1),ylimin(2),simul.nx),	...
				':','linewidth',0.5,'color',[0,0.5,0.5]);
				xlabel('$r$','interpreter','latex');
				ylabel('$\displaystyle \theta_{(S)},\phi_{i(S)},R_{(S)}$','interpreter','latex');
				title(['$\displaystyle Pe = $',num2str(fltr.nondim.Pe)],'interpreter','latex');
				set(gca,'ylim',ylimin);
				saveas(figl,[simul.outfoldrname,'/species_conc_vs_r',	...
				variantname{ivariant},'_Pe',Pename{iPe},'_eta',etaname{ieta}],'fig');
				close(figl);
			end
		end
	end
%	rfS vs eta
	if (plotstatfront == 1)
		soln.eta_plot =			exp(linspace(log(min(eta)),log(max(eta)),simul.nx));
		soln.ana.rfS_plot =		sqrt(soln.eta_plot*Pe(iPe)).*tan(sqrt(soln.eta_plot/Pe(iPe)).*	...
								log(1-((1/(1+fltr.nondim.gamm))*(1-exp((pi/2)*sqrt(Pe(iPe)./soln.eta_plot))))));
		soln.ana.rfS_w_plot=	Pe(iPe)/log(1+fltr.nondim.gamm)*(soln.eta_plot./soln.eta_plot);
		soln.ana.rfS_s_plot =	sqrt(soln.eta_plot*Pe(iPe));
		figl  =	figure('position',[100,100,1250,500],'visible',vistog); axl = axes(figl,'position',[0.125,0.2,0.80,0.7]);
		set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
		set(axl,'xscale','log'); set(axl,'yscale','log');
		set(gca,'xlim',[min(soln.eta_plot),max(soln.eta_plot)]);
		hold on;
		plot(axl,soln.eta_plot,soln.ana.rfS_plot,'k-','linewidth',2.5);
		hold on;
		ylimin =	get(axl,'ylim');
		if(ylimin(1)<0)
			ylimin(1) =		min([ylimin(1)*1.25,-0.25*(ylimin(2)-ylimin(1))]);
		else
			ylimin(1) =		min([ylimin(1)*0.75,-0.25*(ylimin(2)-ylimin(1))]);
		end
		if(ylimin(2)<0)
			ylimin(2) =		max([ylimin(2)*0.75,0.25*(ylimin(2)-ylimin(1))]);
		else
			ylimin(2) =		max([ylimin(2)*1.25,0.25*(ylimin(2)-ylimin(1))]);
		end
		[~,switchieta] =	min((Pe(iPe)-soln.eta_plot).^2);
		plot(axl,soln.eta_plot(switchieta:end),soln.ana.rfS_s_plot(switchieta:end),'g-','linewidth',1.25);
		plot(axl,soln.eta_plot(1:switchieta),soln.ana.rfS_w_plot(1:switchieta),'b-','linewidth',1.25);
		plot(axl,eta,soln.ana.rfS(iPe,:),'k--','linewidth',1.25);
		plot(axl,eta,soln.ana.rfS_w(iPe,:),'g--','linewidth',0.625);
		plot(axl,eta,soln.ana.rfS_s(iPe,:),'b--','linewidth',0.625);
		set(gca,'ylim',ylimin);
		xlabel('$\displaystyle \eta$','interpreter','latex');
		ylabel('$\displaystyle r_{f(S)}$','interpreter','latex');
		title(['$\displaystyle PeQ = $',num2str(Pe(iPe))],'interpreter','latex');
		legend({'full expression','strong dispersion limit','weak dispersion limit'},	...
		'Interpreter','latex','fontsize',12,'location','eastoutside');
		saveas(figl,[simul.outfoldrname,'/rfS_vs_eta_Pe',Pename{iPe}],'fig');
		close(figl);
	end
%	wfS vs eta
	if (plotwidthfront == 1)
		figl  =	figure('position',[100,100,1250,500],'visible',vistog); axl = axes(figl,'position',[0.125,0.2,0.80,0.7]);
		set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
		set(axl,'xscale','log'); set(axl,'yscale','log');
		set(gca,'xlim',[min(eta),max(eta)]);
		hold on;
		plot(axl,eta,soln.num.wfS(iPe,1:length(eta),1,1),'b-','linewidth',2.5);
		plot(axl,eta,soln.num.wfS(iPe,1:length(eta),1,2),'b--','linewidth',2.5);
		plot(axl,eta,soln.num.wfS(iPe,1:length(eta),2,1),'b-.','linewidth',2.5);
		plot(axl,eta,soln.num.wfS(iPe,1:length(eta),2,2),'b:','linewidth',2.5);
		if (solvesimil == 1)
			plot(axl,eta,soln.simil.wfS(iPe,1:length(eta),1),'r--','linewidth',1.25);
			plot(axl,eta,soln.simil.wfS(iPe,1:length(eta),2),'r-.','linewidth',1.25);
			plot(axl,eta,soln.simil.wfS(iPe,1:length(eta),3),'r:','linewidth',1.25);
		end
		ylimin =	get(axl,'ylim');
		if(ylimin(1)<0)
			ylimin(1) =		min([ylimin(1)*1.25,-0.25*(ylimin(2)-ylimin(1))]);
		else
			ylimin(1) =		min([ylimin(1)*0.75,-0.25*(ylimin(2)-ylimin(1))]);
		end
		if(ylimin(2)<0)
			ylimin(2) =		max([ylimin(2)*0.75,0.25*(ylimin(2)-ylimin(1))]);
		else
			ylimin(2) =		max([ylimin(2)*1.25,0.25*(ylimin(2)-ylimin(1))]);
		end
		set(gca,'ylim',ylimin);
		xlabel('$\displaystyle \eta$','interpreter','latex');
		ylabel('$\displaystyle w_{f(S)}$','interpreter','latex');
		title(['$\displaystyle PeQ = $',num2str(Pe(iPe))],'interpreter','latex');
		legend({'$r^2R$','$R$','$r^2\theta$','$\theta$'},	...
		'Interpreter','latex','fontsize',12,'location','eastoutside');
		saveas(figl,[simul.outfoldrname,'/wfS_vs_eta_Pe',Pename{iPe}],'fig');
		close(figl);
	end
%	barR vs eta
	if (plotbarR == 1)
		figl  =	figure('position',[100,100,1250,500],'visible',vistog); axl = axes(figl,'position',[0.125,0.2,0.80,0.7]);
		set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
		set(axl,'xscale','log'); set(axl,'yscale','log');
		set(gca,'xlim',[min(eta),max(eta)]);
		hold on;
		plot(axl,eta,soln.num.barRS(iPe,1:length(eta)),'b-','linewidth',2.5);
		if (solvesimil == 1)
			plot(axl,eta,soln.simil.barRS(iPe,1:length(eta)),'b--','linewidth',1.25);
		end
		ylimin =	get(axl,'ylim');
		if(ylimin(1)<0)
			ylimin(1) =		min([ylimin(1)*1.25,-0.25*(ylimin(2)-ylimin(1))]);
		else
			ylimin(1) =		min([ylimin(1)*0.75,-0.25*(ylimin(2)-ylimin(1))]);
		end
		if(ylimin(2)<0)
			ylimin(2) =		max([ylimin(2)*0.75,0.25*(ylimin(2)-ylimin(1))]);
		else
			ylimin(2) =		max([ylimin(2)*1.25,0.25*(ylimin(2)-ylimin(1))]);
		end
		set(gca,'ylim',ylimin);
		xlabel('$\displaystyle \eta$','interpreter','latex');
		ylabel('$\displaystyle \bar{R}_{(S)}$','interpreter','latex');
		title(['$\displaystyle PeQ = $',num2str(Pe(iPe))],'interpreter','latex');
		saveas(figl,[simul.outfoldrname,'/barRS_vs_eta_Pe',Pename{iPe}],'fig');
		close(figl);
	end
end
save([simul.outfoldrname,'/Pe',Pename{iPe},'.mat']);
%-----------------------------------------------------------------------------------------------------------------------------------