%-----------------------------------------------------------------------------------------------------------------------------------
% initiating cleaned matlab sessions
clc;
clear;
close all;
fclose all;
warning('off','MATLAB:nearlySingularMatrix');
addpath('../simroutines/functions');
addpath('./anasolvers');
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% specifying parameter values and names
n_eta_in =						41;
n_Pe_in =						21;
etalabel =						'$\log{\eta}$';
Pelabel =						'$\log{Pe}$';
contrastname =					{'ratio', 'diff', 'normdiff','baseval'};
propname =						{'rfS', 'wfS', 'barRS','rfS_num'};
eta =							[0.0,10.^(linspace(-3,6,n_eta_in-1))];
Pe =							flip(10.^(linspace(0,3,n_Pe_in)));
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
simul.outfoldrname =			'../outdata/SIStatHMScaling_gamma_1e1';
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
simul.xgridskew =				3.0;				%	skewing parameter for xgrid [SI -> Pe end eta end]
													%	(higher than unity index - denser near start)
													%	(lower than unity index - denser near end)
simul.STLtol =					1.00E-05;			%	tolerance of update diminution for Source Term Linearization Iterations
simul.STLrelax =				0.80;				%	relaxation of solution update
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
simul.dxRmax =					1.00E+04;			%	maximum allowed dx at right end (enforced when simul.limitdx is 1)
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
geomdom.size.Lx =				1.00E+08;			%	non-dimensionalized x-coordinate length
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
% 2 -	stochastic permeability generated based on stochasticticity specifications (using a function call to permgen)
% Flow & Transport Parameters
fltr.nondim.gamm =				10.0;				%	concentration of B at far end
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
soln.ana.rfS =					zeros(length(Pe),length(eta));
soln.ana.thetaS =				zeros(simul.nx,length(Pe),length(eta));
soln.ana.rfS_s =				zeros(length(Pe),length(eta));
soln.ana.rfS_w =				zeros(length(Pe),length(eta));
soln.ana.thetaS_s =				zeros(simul.nx,length(Pe),length(eta));
soln.ana.thetaS_w =				zeros(simul.nx,length(Pe),length(eta));
%-----------------------------------------------------------------------------------------------------------------------------------
% initializing contrast matrix
contrast_FrProp =				zeros(length(Pe),length(eta),4,3);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------\
% obtaining flow solution
if (simul.givenvel == 0)
	[soln.flow.H, soln.flow.vel, fltr.perm.Kperm] =		...
								solver1D_darcy(geomdom,simul,fltr,grd);
else
	fltr.perm.Kperm =			fltr.perm.Kperm0*ones(simul.nx,1);
	soln.flow.H =				zeros(simul.nx,1);
	soln.flow.vel.v =			double(geomdom.isradial~=0)*(1./(grd.x.^double(geomdom.isradial)));
end
% running simulations for transport
for iPe = 1:length(Pe)
%	Pre-Allocating Iteration-Specific Matrices
	for ieta = 1:length(eta)
%		setting paramter values for current loop	
		fltr.nondim.Pe =				Pe(iPe);			%	expected Peclet number
		fltr.eta =						eta(ieta);			%	expected longitudinal dispersivity factor
%		analytical front location
		soln.ana.rfS(iPe,ieta) =		sqrt(fltr.eta*fltr.nondim.Pe)*tan(sqrt(fltr.eta/(fltr.nondim.Pe))*				...
										log(1-(1/(1+fltr.nondim.gamm))*(1-exp((pi/2)*sqrt(fltr.nondim.Pe/fltr.eta)))));
		soln.ana.rfS_w(iPe,ieta) =		fltr.nondim.Pe/log(1+fltr.nondim.gamm);
		soln.ana.rfS_S(iPe,ieta) =		sqrt(fltr.eta*fltr.nondim.Pe);
%		analytical stationary solution for theta
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
		if (ieta == 1)
			fltr.init.phiA =			0.5*(1.0-tanh((grd.x-soln.ana.rfS(iPe,ieta))/fltr.init.front.k2));
			fltr.init.phiB =			fltr.nondim.gamm*0.5*(1.0+tanh((grd.x-soln.ana.rfS(iPe,ieta))/fltr.init.front.k2));
		else
			fltr.init.phiA =			phitemp(1:simul.nx,1);
			fltr.init.phiB =			phitemp(1:simul.nx,2);
		end
%		numerical stationary solution
		phitemp =						stationary1D_reactive(geomdom,simul,fltr,grd,soln.flow);
		barRS =							4*pi*trapz(grd.x,grd.x.^2.*phitemp(1:simul.nx,1)'.*phitemp(1:simul.nx,2)');
		[~,ind_fS] =					min((phitemp(:,1)-phitemp(:,2)).^2);
		rfS =							grd.x(ind_fS);
		IntNr =							trapz(grd.x,grd.x.^2.*(grd.x-rfS).^2.*phitemp(1:simul.nx,1)'.*phitemp(1:simul.nx,2)');
		wfS =							sqrt(IntNr/barRS);
%		writing contrast values
		if (ieta == 1)
			rfS_base = rfS; barRS_base = barRS; wfS_base = wfS; rfS_ana_base = soln.ana.rfS(iPe,ieta);
		end
		contrast_FrProp(iPe,ieta,1,1) =	soln.ana.rfS(iPe,ieta)/rfS_ana_base;
		contrast_FrProp(iPe,ieta,1,2) =	wfS/wfS_base;
		contrast_FrProp(iPe,ieta,1,3) =	barRS/barRS_base;
		contrast_FrProp(iPe,ieta,1,4) =	rfS/rfS_base;
		contrast_FrProp(iPe,ieta,2,1) =	soln.ana.rfS(iPe,ieta)-rfS_ana_base;
		contrast_FrProp(iPe,ieta,2,2) =	wfS-wfS_base;
		contrast_FrProp(iPe,ieta,2,3) =	barRS-barRS_base;
		contrast_FrProp(iPe,ieta,2,4) =	rfS-rfS_base;
		contrast_FrProp(iPe,ieta,3,1) =	soln.ana.rfS(iPe,ieta)/rfS_ana_base-1;
		contrast_FrProp(iPe,ieta,3,2) =	wfS/wfS_base-1;
		contrast_FrProp(iPe,ieta,3,3) =	barRS/barRS_base-1;
		contrast_FrProp(iPe,ieta,3,4) =	rfS/rfS_base-1;
		contrast_FrProp(iPe,ieta,4,1) =	soln.ana.rfS(iPe,ieta);
		contrast_FrProp(iPe,ieta,4,2) =	wfS;
		contrast_FrProp(iPe,ieta,4,3) =	barRS;
		contrast_FrProp(iPe,ieta,4,4) =	rfS;
		disp(['Done with Pe = ',num2str(Pe(iPe)),', eta = ',num2str(eta(ieta))]);			
	end
end
% generating heatmaps
[logetaeta, logPePe] = meshgrid(log10(eta(2:end)),log10(Pe));
for icontrast = 1:4
	for iprop = 1:4
		pcolor(logetaeta,logPePe,real(log10(contrast_FrProp(:,2:end,icontrast,iprop)))); shading interp;
		saveas(gcf,[simul.outfoldrname,'/SI_HM_',contrastname{icontrast},'_',propname{iprop},'_Log_Real'],'fig');
		close(gcf);
		pcolor(logetaeta,logPePe,imag(log10(contrast_FrProp(:,2:end,icontrast,iprop)))); shading interp;
		saveas(gcf,[simul.outfoldrname,'/SI_HM_',contrastname{icontrast},'_',propname{iprop},'_Log_Imag'],'fig');
		close(gcf);
		pcolor(logetaeta,logPePe,contrast_FrProp(:,2:end,icontrast,iprop)); shading interp;
		saveas(gcf,[simul.outfoldrname,'/SI_HM_',contrastname{icontrast},'_',propname{iprop}],'fig');
		close(gcf);
		figl  =	figure('position',[100,35,1250,700],'visible','on'); axl = axes(figl,'position',[0.125,0.15,0.85,0.800]);
		set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex');
		set(axl,'xscale','log'); set(axl,'yscale','log');
		set(gca,'linewidth',3); hold on;
		xlabel(etalabel,'interpreter','latex');
		for iPe = 1:length(Pe)
			plot(eta(2:end),contrast_FrProp(iPe,2:end,icontrast,iprop),'-','color',(iPe/(2*length(Pe)))*[1,1,1],'linewidth',2);
		end
		saveas(figl,[simul.outfoldrname,'/SI_etaplots_',contrastname{icontrast},'_',propname{iprop}],'fig');
		close(figl);
		figl  =	figure('position',[100,35,1250,700],'visible','on'); axl = axes(figl,'position',[0.125,0.15,0.85,0.800]);
		set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex');
		set(axl,'xscale','log'); set(axl,'yscale','log');
		set(gca,'linewidth',3); hold on;
		xlabel('$\eta/Pe$','interpreter','latex');
		if ((iprop == 1) || (iprop == 4))
			ylabel('$r_{f(S)}/Pe$','interpreter','latex');
			for iPe = 1:length(Pe)
				plot(eta(2:end)/Pe(iPe),contrast_FrProp(iPe,2:end,icontrast,iprop)/Pe(iPe),		...
				'-','color',(iPe/(2*length(Pe)))*[1,1,1],'linewidth',2);
			end
		elseif (iprop == 3)
			ylabel('$\bar{R}_{(S)}$','interpreter','latex');
			for iPe = 1:length(Pe)
				plot(eta(2:end)/Pe(iPe),contrast_FrProp(iPe,2:end,icontrast,iprop),		...
				'-','color',(iPe/(2*length(Pe)))*[1,1,1],'linewidth',2);
			end
		end
		saveas(figl,[simul.outfoldrname,'/SI_etaplots_collapse',contrastname{icontrast},'_',propname{iprop}],'fig');
		close(figl);
		figl  =	figure('position',[100,35,1250,700],'visible','on'); axl = axes(figl,'position',[0.125,0.15,0.85,0.800]);
		set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
		set(axl,'xscale','log'); set(axl,'yscale','log');
		set(gca,'linewidth',3); hold on;
		xlabel(Pelabel,'interpreter','latex');
		for ieta = 1:length(eta)
			plot(Pe,contrast_FrProp(:,ieta,icontrast,iprop),'color',(ieta/(2*length(eta)))*[1,1,1],'linewidth',2);
		end
		saveas(figl,[simul.outfoldrname,'/SI_Peplots_',contrastname{icontrast},'_',propname{iprop}],'fig');
		close(figl);
	end
end
% saving output
save([simul.outfoldrname,'/SI_HM_data.mat']);