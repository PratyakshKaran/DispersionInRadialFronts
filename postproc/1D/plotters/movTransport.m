% generating transport movie

if (postproc.iTransportSS == 1)
	mkdir([postproc.folderloc1D,'/phiSS']);
end
phimovvid = VideoWriter([postproc.folderloc1D,'/Transport_Species.avi']);
phimovvidA = VideoWriter([postproc.folderloc1D,'/Transport_Auxiliary.avi']);
open(phimovvid);
open(phimovvidA);
figl  =	figure('position',[100,35,1250,700],'visible',postproc.vistogmov);
if (soln1D.geomdom.isradial ~= 0)
	axl = axes(figl,'position',[0.125,0.125,0.80,0.325]); hold on;
	axllin = axes(figl,'position',[0.125,0.600,0.80,0.35]); hold on;
else
	axl = axes(figl,'position',[0.125,0.15,0.85,0.800]); hold on;
end
if (soln1D.geomdom.isradial ~= 0)
	set(axllin,'box','on'); grid(postproc.gridtog); set(axllin,'fontsize',25);set(axllin,'ticklabelinterpreter','latex');
end
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25);set(axl,'ticklabelinterpreter','latex');
if (soln1D.geomdom.isradial == 0)
	xliml =			0.0;
	xlimr =			soln1D.geomdom.size.Lx;
else
	xliml =			max([soln1D.geomdom.size.a,postproc.maxxrange_log*soln1D.geomdom.size.A]);
	xlimr =			soln1D.geomdom.size.A;
end
if (soln1D.geomdom.isradial ~= 0)
	set(axllin,'xlim',[soln1D.geomdom.size.a,soln1D.geomdom.size.A]);
	set(axllin,'ylim',[-1.25,1.25]);
end
ylims = [-1.25,1.25];
set(axl,'xlim',[xliml,xlimr]);
set(axl,'ylim',ylims);
if (soln1D.geomdom.isradial ~= 0)
	set(axl,'xscale','log');
end
if (soln1D.geomdom.isradial == 0)
	xlabel('$x$','interpreter','latex');
else
	xlabel('$r$','interpreter','latex');
end
ylabel('$\phi$','interpreter','latex');
if (postproc.iscomparison1D == 1)
	phiCcurrentcomparison1D =	load([postproc.folderloccompare1D,'/phiC.dat']);
	phiAcurrentcomparison1D =	load([postproc.folderloccompare1D,'/phiA.dat']);
	phiBcurrentcomparison1D =	load([postproc.folderloccompare1D,'/phiB.dat']);
end
ifront.tsaved = 0;
%	plotting regime switches
plot(axl,[xliml,xlimr],[0,0],'k:','linewidth',1);
plot(axl,	...
min(xlimr,((soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^(1/double(soln1D.geomdom.isradial))))*[1,1],	...
[-1.25,1.25],'--','color',[0.2,0.9,0.9],'linewidth',2);
%	plotting on linear axis
if (soln1D.geomdom.isradial ~= 0)
	plot(axllin,[xliml,xlimr],[0,0],'k:','linewidth',1);
	plot(axllin,	...
	min(xlimr,((soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^(1/double(soln1D.geomdom.isradial))))*[1,1],...
	[-1.25,1.25],'--','color',[0.2,0.9,0.9],'linewidth',2);
end
itrangemov =		1+double(postproc.isexpected):length(soln1D.simul.tstepsave);
itrangeSS =			itrangemov(floor(linspace(1+double(postproc.isexpected),	...
					length(soln1D.simul.tstepsave),postproc.nTransportSS)));
if (postproc.iTransport == 0)
	itrange =		itrangeSS;
else
	itrange =		itrangemov;
end
for it = soln1D.simul.tstepsave(itrange)
% 	try
		ifront.tsaved =								ifront.tsaved+1;
%		analytical solutions
		if (postproc.obtainanalytical == 1)
			if (soln1D.geomdom.isradial == 0)
%				similarity solutions
				front.ana.x_simil =					soln1D.grd.t(it)^(1/6)*front.ana.z_simil+front.ana.r_f(it);
				front.ana.R_simil =					soln1D.solnsimil.R*(soln1D.grd.t(it)^(-2/3));
				front.ana.theta_simil =				-((front.ana.x_simil-front.ana.r_f(it))*soln1D.solnsimil.K)/	...
													sqrt(soln1D.grd.t(it));
				front.ana.phiA_simil =				soln1D.grd.t(it)^(-1/3)*front.ana.G_simil;
				front.ana.phiB_simil =				front.ana.phiA_simil-front.ana.theta_simil;
				front.ana.x_simil_ET =				soln1D.grd.t(it)^(1/6)*front.ana.z_simil_ET+front.ana.r_f(it);
				front.ana.R_simil_ET =				soln1D.solnsimil_ET.R*(soln1D.grd.t(it)^(-2/3));
				front.ana.theta_simil_ET =			-((front.ana.x_simil_ET-front.ana.r_f(it))*soln1D.solnsimil_ET.K)/		...
													sqrt(soln1D.grd.t(it));
				front.ana.phiA_simil_ET =			soln1D.grd.t(it)^(-1/3)*front.ana.G_simil_ET;
				front.ana.phiB_simil_ET =			front.ana.phiA_simil_ET-front.ana.theta_simil_ET;
			elseif (soln1D.geomdom.isradial == 1)
%				similarity solutions
				front.ana.x_simil =					soln1D.grd.t(it)^(1/6)*front.ana.z_simil+front.ana.r_f(it);
				front.ana.R_simil =					soln1D.solnsimil.R*(soln1D.grd.t(it)^(-2/3));
				front.ana.theta_simil =				-((front.ana.x_simil-front.ana.r_f(it))*soln1D.solnsimil.K)/	...
													sqrt(soln1D.grd.t(it));
				front.ana.phiA_simil =				soln1D.grd.t(it)^(-1/3)*front.ana.G_simil;
				front.ana.phiB_simil =				front.ana.phiA_simil-front.ana.theta_simil;
				front.ana.x_simil_ET =				soln1D.grd.t(it)^(1/2)*front.ana.z_simil_ET+front.ana.r_f(it);
				front.ana.R_simil_ET =				soln1D.solnsimil_ET.R;
				front.ana.theta_simil_ET =			-((front.ana.x_simil_ET-front.ana.r_f(it))*soln1D.solnsimil_ET.K)/	...
													sqrt(soln1D.grd.t(it));
				front.ana.phiA_simil_ET =			front.ana.G_simil_ET;
				front.ana.phiB_simil_ET =			front.ana.phiA_simil_ET-front.ana.theta_simil_ET;
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					front.ana.x_simil_disp =		front.ana.z_simil_disp+front.ana.r_f_disp(it);
					front.ana.R_simil_disp =		soln1D.solnsimil_disp.R*(soln1D.grd.t(it)^(-2/3));
					front.ana.theta_simil_disp =	-((front.ana.x_simil_disp-front.ana.r_f_disp(it))*soln1D.solnsimil_disp.K)/	...
													(soln1D.grd.t(it)^(1/3));
					front.ana.phiA_simil_disp =		soln1D.grd.t(it)^(-1/3)*front.ana.G_simil_disp;
					front.ana.phiB_simil_disp =		front.ana.phiA_simil_disp-front.ana.theta_simil_disp;
					front.ana.x_simil_disp_ET =		soln1D.grd.t(it)^(1/3)*front.ana.z_simil_disp_ET+front.ana.r_f_disp(it);
					front.ana.R_simil_disp_ET =		soln1D.solnsimil_disp_ET.R;
					front.ana.theta_simil_disp_ET =	-((front.ana.x_simil_disp_ET-front.ana.r_f_disp(it))*	...
													soln1D.solnsimil_disp_ET.K)/(soln1D.grd.t(it)^(1/3));
					front.ana.phiA_simil_disp_ET =	front.ana.G_simil_disp_ET;
					front.ana.phiB_simil_disp_ET =	front.ana.phiA_simil_disp_ET-front.ana.theta_simil_disp_ET;
				end
			else
%				similarity solutions
				front.ana.x_simil =					soln1D.solnsimil.K^(-1/3)*front.ana.z_simil+front.ana.r_f(it);
				front.ana.R_simil =					soln1D.solnsimil.RS;
				front.ana.theta_simil =				soln1D.solnsimil.thetaS;
				front.ana.phiA_simil =				soln1D.solnsimil.phiAS;
				front.ana.phiB_simil =				soln1D.solnsimil.phiBS;
				front.ana.x_simil_ET =				soln1D.grd.t(it)^(1/2)*front.ana.z_simil+front.ana.r_f_ET(it);
				front.ana.R_simil_ET =				soln1D.solnsimil_ET.R;
				front.ana.theta_simil_ET =			soln1D.solnsimil_ET.theta;
				front.ana.phiA_simil_ET =			soln1D.solnsimil_ET.phiA;
				front.ana.phiB_simil_ET =			soln1D.solnsimil_ET.phiB;
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))					
					front.ana.x_simil_disp =		front.ana.x_simil;
					front.ana.R_simil_disp =		front.ana.R_simil;
					front.ana.theta_simil_disp =	front.ana.theta_simil;
					front.ana.phiA_simil_disp =		front.ana.phiA_simil;
					front.ana.phiB_simil_disp =		front.ana.phiB_simil;
					front.ana.r_f_disp =			pchip(front.tsaved,front.num.r_f,soln1D.grd.t);
					front.ana.x_simil_disp_ET =		soln1D.grd.t(it)^(1/4)*front.ana.z_simil_disp_ET+front.ana.r_f_disp(it);
					front.ana.R_simil_disp_ET =		soln1D.solnsimil_disp_ET.R;
					front.ana.theta_simil_disp_ET =	soln1D.solnsimil_disp_ET.theta;
					front.ana.phiA_simil_disp_ET =	soln1D.solnsimil_disp_ET.phiA;
					front.ana.phiB_simil_disp_ET =	soln1D.solnsimil_disp_ET.phiB;
				end
			end
		end
%		expected solutions (corresponding to conservative solutions)
		if (postproc.isexpected == 1)
			phiAexpected =	((1.0/sqrt(4*pi*(1.0/soln1D.fltr.nondim.Pe)*(soln1D.grd.t(it))))/	...
							(1.0/sqrt(4*pi*(1.0/soln1D.fltr.nondim.Pe)*(soln1D.simul.t0))))*	...
							exp(-(((soln1D.grd.x-	...
							((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)/soln1D.geomdom.size.Lx)*soln1D.grd.t(it))-	...
							0.5*soln1D.geomdom.size.Lx).^2)/(4.0*(1.0/soln1D.fltr.nondim.Pe)*(soln1D.grd.t(it))));
			phiBexpected =	zeros(1,soln1D.simul.nx);
			for ifourier = 1:postproc.nfourier
				phiBexpected =	...
							phiBexpected+exp(-(2*ifourier-1)^2*pi^2*soln1D.grd.t(it))*	...
							(sin((2*ifourier-1)*pi*	...
							((soln1D.grd.x-	...
							((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)/soln1D.geomdom.size.Lx)*soln1D.grd.t(it))/	...
							soln1D.geomdom.size.Lx))/(2*ifourier-1));
			end
			phiBexpected =	(4.0/pi)*phiBexpected;
			phiCexpected =	0.5*(1.0-erf((0.5*soln1D.geomdom.size.Lx-	...
							((soln1D.grd.x-	...
							((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)/soln1D.geomdom.size.Lx)*	...
							soln1D.grd.t(it))))/	...
							(sqrt((4*(1.0/soln1D.fltr.nondim.Pe)*(soln1D.grd.t(it)))))));
		end
%		accumulating numerical solutions
		phiCcurrent =		load([postproc.folderloc1D,'/phiC_',num2str(it),'.dat']);
		phiAcurrent =		load([postproc.folderloc1D,'/phiA_',num2str(it),'.dat']);
		phiBcurrent =		load([postproc.folderloc1D,'/phiB_',num2str(it),'.dat']);
%		obtaining warped time numerical solutions
		if (postproc.iscomparison2D == 1)
			phiCcurrentcomparison2D =	load([postproc.folderloccompare2D,'/phiC_',num2str(it),'.dat']);
			phiAcurrentcomparison2D =	load([postproc.folderloccompare2D,'/phiA_',num2str(it),'.dat']);
			phiBcurrentcomparison2D =	load([postproc.folderloccompare2D,'/phiB_',num2str(it),'.dat']);
			if (soln1D.geomdom.isradial == 1)
				r_2D =			zeros(size(phiAcurrentcomparison2D,1)*size(phiAcurrentcomparison2D,2),1);
				phiA_2D =		zeros(size(phiAcurrentcomparison2D,1)*size(phiAcurrentcomparison2D,2),1);
				phiB_2D =		zeros(size(phiAcurrentcomparison2D,1)*size(phiAcurrentcomparison2D,2),1);
				phiC_2D =		zeros(size(phiAcurrentcomparison2D,1)*size(phiAcurrentcomparison2D,2),1);
				for ix = 1:size(phiAcurrentcomparison2D,2)
					for iy = 1:size(phiAcurrentcomparison2D,1)
						r_2D((ix-1)*size(phiAcurrentcomparison2D,1)+iy) =	...
														sqrt((solncompare2D.grd.x(ix)-0.5*solncompare2D.geomdom.size.Lx)^2+	...
														(solncompare2D.grd.y(iy)-0.5*solncompare2D.geomdom.size.Ly)^2);
						phiC_2D((ix-1)*size(phiAcurrentcomparison2D,1)+iy) =	...
														phiCcurrentcomparison2D(iy,ix);
						phiB_2D((ix-1)*size(phiAcurrentcomparison2D,1)+iy) =	...
														phiBcurrentcomparison2D(iy,ix);
						phiA_2D((ix-1)*size(phiAcurrentcomparison2D,1)+iy) =	...
														phiAcurrentcomparison2D(iy,ix);
					end
				end
			end
		end
%		numerical solutions for reaction rate and theta
		if ((postproc.iFront == 1) || (postproc.zoomin == 1) || (postproc.obtainanalytical == 1))
			front.num.R =							phiAcurrent.*phiBcurrent;
			front.num.theta =						phiAcurrent-phiBcurrent;
		end
%		analytical solution for theta
		if (soln1D.geomdom.isradial == 1)
			if (soln1D.fltr.nondim.Pe < postproc.thresPediff)
				front.ana.theta =					-soln1D.fltr.nondim.gamm+(1+soln1D.fltr.nondim.gamm)*					...
													gammainc(soln1D.grd.x.^2/(4*soln1D.grd.t(it)),soln1D.fltr.nondim.Pe/2,	...
													'upper');
			else
				YbyA =								(soln1D.grd.x.^2)/(2*soln1D.fltr.nondim.Pe*soln1D.grd.t(it));
				M =									YbyA-1;
				omga =								(M./abs(M)).*sqrt(2*(YbyA-1-log(YbyA)));
				front.ana.theta =					-soln1D.fltr.nondim.gamm+(1+soln1D.fltr.nondim.gamm)*					...
													0.5*erfc(omga*sqrt(soln1D.fltr.nondim.Pe/4));			
			end
			front.ana.theta_disp =					-soln1D.fltr.nondim.gamm+(1+soln1D.fltr.nondim.gamm)*					...
													gammainc(soln1D.grd.x.^3/(9*soln1D.fltr.eta*soln1D.grd.t(it)),1/3,		...
													'upper');
		end
%		postproc.zooming in
		if (postproc.zoomin == 1)
			leftswitch =				1;
			rightswitch =				1;
			ixleft =					1;
			ixright =					soln1D.simul.nx;
			for ix = 1:soln1D.simul.nx
				if ((front.num.R(ix) >= 1E-6*max(front.num.R)) && (leftswitch == 1))
					ixleft =			ix;
					leftswitch =		0;
				end
			end
			for ix = soln1D.simul.nx:-1:1
				if ((front.num.R(ix) >= 1E-6*max(front.num.R)) && (rightswitch == 1))
					ixright =			ix;
					rightswitch =		0;
				end
			end
			if (ixright == ixleft)
				ixright =				ixleft+1;
			end
			set(axl,'xlim',[soln1D.grd.x(ixleft),soln1D.grd.x(ixright)]);
			set(axllin,'xlim',[soln1D.grd.x(ixleft),soln1D.grd.x(ixright)]);
		end
%		including title
		title([postproc.casename,': time is ',num2str(soln1D.grd.t(it))],'Interpreter','none','FontSize',12);
%		plotting collected solutions
		npl = 0;
%		plotting advective element location
		plot(axl,front.r_adv(it)*[1,1],1.25*[-1,1],'-','color',[0.35,0.35,0.35],'linewidth',3); npl=npl+1;
%		plotting analytical front locations
		if (postproc.obtainanalytical == 1)
			plot(axl,front.ana.r_f(it)*[1,1],1.25*[-1,1],'-.','color',[0.35,0.35,0.35],'linewidth',3); npl=npl+1;
			if ((soln1D.fltr.eta ~= 0) && ...
			((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
				plot(axl,front.ana.r_f_disp(it)*[1,1],1.25*[-1,1],'--','color',[0.35,0.35,0.35],'linewidth',3); npl=npl+1;
			end
		end
%		filling front widths
		if (postproc.plotfrontwidth == 1)
%			filling analytical front widths
			if (postproc.obtainanalytical == 1)
				fill(axl,	...
				[max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f(it)),	...
				min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f(it)),		...
				min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f(it)),		...
				max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f(it))],	...
				[ylims(1),ylims(1),ylims(2),ylims(2)],	...
				[0.35,0.35,0.35],'LineStyle','-.','LineWidth',1.5,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); npl=npl+1;
				fill(axl,	...
				[max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f_approx(it)),	...
				min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f_approx(it)),	...
				min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f_approx(it)),	...
				max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f_approx(it))],	...
				[ylims(1),ylims(1),ylims(2),ylims(2)],	...
				[0.35,0.35,0.35],'LineStyle','-.','LineWidth',2,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); npl=npl+1;	
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					fill(axl,	...
					[max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp(it)),	...
					min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp(it)),	...
					min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp(it)),	...
					max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp(it))],	...
					[ylims(1),ylims(1),ylims(2),ylims(2)],	...
					[0.35,0.35,0.35],'LineStyle','--','LineWidth',1.5,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); npl=npl+1;
					fill(axl,	...
					[max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp_approx(it)),	...
					min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp_approx(it)),	...
					min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp_approx(it)),	...
					max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp_approx(it))],	...
					[ylims(1),ylims(1),ylims(2),ylims(2)],	...
					[0.35,0.35,0.35],'LineStyle','--','LineWidth',2,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); npl=npl+1;	
				end
				fill(axl,	...
				[max(xliml,front.num.r_f(ifront.tsaved)-0.5*front.num.w_f(ifront.tsaved)),		...
				min(xlimr,front.num.r_f(ifront.tsaved)+0.5*front.num.w_f(ifront.tsaved)),		...
				min(xlimr,front.num.r_f(ifront.tsaved)+0.5*front.num.w_f(ifront.tsaved)),		...
				max(xliml,front.num.r_f(ifront.tsaved)-0.5*front.num.w_f(ifront.tsaved))],		...
				[ylims(1),ylims(1),ylims(2),ylims(2)],	...
				[0.35,0.35,0.35],'LineStyle','-','LineWidth',0.5,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); npl=npl+1;	
			end
			fill(axl,	...
			[max(xliml,front.num.r_f_orig(ifront.tsaved)-0.5*front.num.w_f_orig(ifront.tsaved,2,1)),	...
			min(xlimr,front.num.r_f_orig(ifront.tsaved)+0.5*front.num.w_f_orig(ifront.tsaved,2,1)),		...
			min(xlimr,front.num.r_f_orig(ifront.tsaved)+0.5*front.num.w_f_orig(ifront.tsaved,2,1)),		...
			max(xliml,front.num.r_f_orig(ifront.tsaved)-0.5*front.num.w_f_orig(ifront.tsaved,2,1))],	...
			[ylims(1),ylims(1),ylims(2),ylims(2)],	...
			[0.35,0.35,0.35],'LineStyle','-','LineWidth',1.0,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); npl=npl+1;
		end
		npl2retain = npl;
%		plotting numerical species solution
		plot(axl,soln1D.grd.x,phiCcurrent/max(phiCcurrent),'-','color',[1.0,0.2,0.2],'linewidth',2); npl=npl+1;
		plot(axl,soln1D.grd.x,phiBcurrent/max(1,soln1D.fltr.nondim.gamm),'-','color',[0.2,0.2,1.0],'linewidth',2); npl=npl+1;
		plot(axl,soln1D.grd.x,phiAcurrent/max(1,soln1D.fltr.nondim.gamm),'-','color',[0.2,1.0,0.2],'linewidth',2); npl=npl+1;
		if (postproc.isexpected == 1)
			plot(axl,soln1D.grd.x,phiCexpected/max(phiCcurrent),'--','color',[1.0,0.8,0.8],'linewidth',3); npl=npl+1;
			plot(axl,soln1D.grd.x,phiBexpected/max(1,soln1D.fltr.nondim.gamm),'--','color',[0.8,0.8,1.0],'linewidth',3); npl=npl+1;
			plot(axl,soln1D.grd.x,phiAexpected/max(1,soln1D.fltr.nondim.gamm),'--','color',[0.8,1.0,0.8],'linewidth',3); npl=npl+1;
		end
		if (postproc.iscomparison1D == 1)
			plot(axl,solncompare1D.grd.x,phiCcurrentcomparison1D(it,:)/max(phiCcurrent),	...
			'--','color',[1.0,0.4,0.4],'linewidth',4); npl=npl+1;
			plot(axl,solncompare1D.grd.x,phiBcurrentcomparison1D(it,:)/max(1,soln1D.fltr.nondim.gamm),		...
			'--','color',[0.4,0.4,1.0],'linewidth',4); npl=npl+1;
			plot(axl,solncompare1D.grd.x,phiAcurrentcomparison1D(it,:)/max(1,soln1D.fltr.nondim.gamm),		...
			'--','color',[0.4,1.0,0.4],'linewidth',4); npl=npl+1;
		end
		if (postproc.iscomparison2D == 1)
			if (soln1D.geomdom.isradial == 0)
				plot(axl,solncompare2D.grd.x,	...
				phiCcurrentcomparison2D((solncompare2D.simul.ny-1)/2+1,:)/max(phiCcurrent),	...
				':','color',[1.0,0.4,0.4],'linewidth',4); npl=npl+1;
				plot(axl,solncompare2D.grd.x,	...
				phiBcurrentcomparison2D((solncompare2D.simul.ny-1)/2+1,:)/max(1,soln1D.fltr.nondim.gamm),	...
				':','color',[0.4,0.4,1.0],'linewidth',4); npl=npl+1;
				plot(axl,solncompare2D.grd.x,	...
				phiAcurrentcomparison2D((solncompare2D.simul.ny-1)/2+1,:)/max(1,soln1D.fltr.nondim.gamm),	...
				':','color',[0.4,1.0,0.4],'linewidth',4); npl=npl+1;
			else
				plot(axl,r_2D,phiC_2D/max(phiCcurrent),'.','color',[1.0,0.4,0.4],'linewidth',1,'markersize',1.25); 
				npl=npl+1;
				plot(axl,r_2D,phiB_2D/max(1,soln1D.fltr.nondim.gamm),'.','color',[0.4,0.4,1.0],'linewidth',1,'markersize',1.25);
				npl=npl+1;
				plot(axl,r_2D,phiA_2D/max(1,soln1D.fltr.nondim.gamm),'.','color',[0.4,1.0,0.4],'linewidth',1,'markersize',1.25);
				npl=npl+1;
			end
		end
%		plotting analytical species solution
		if (postproc.obtainanalytical == 1)
			if (soln1D.grd.t(it) < 1.0)
				plot(axl,front.ana.x_simil_ET,front.ana.phiA_simil_ET/max(abs(phiAcurrent)),		...
				'-.','color',[0.2,1.0,0.2],'linewidth',1);
				npl=npl+1;
				plot(axl,front.ana.x_simil_ET,front.ana.phiB_simil_ET/max(abs(phiBcurrent)),		...
				'-.','color',[0.2,0.2,1.0],'linewidth',1);
				npl=npl+1;
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					plot(axl,front.ana.x_simil_disp_ET,front.ana.phiA_simil_disp_ET/max(abs(phiAcurrent)),	...
					'--','color',[0.2,1.0,0.2],'linewidth',1);	
					npl=npl+1;
					plot(axl,front.ana.x_simil_disp_ET,front.ana.phiB_simil_disp_ET/max(abs(phiBcurrent)),	...
					'--','color',[0.2,0.2,1.0],'linewidth',1);
					npl=npl+1;
				end
			else
				plot(axl,front.ana.x_simil,front.ana.phiA_simil/max(abs(phiAcurrent)),				...
				'-.','color',[0.2,1.0,0.2],'linewidth',1);
				npl=npl+1;
				plot(axl,front.ana.x_simil,front.ana.phiB_simil/max(abs(phiBcurrent)),				...
				'-.','color',[0.2,0.2,1.0],'linewidth',1);
				npl=npl+1;
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					plot(axl,front.ana.x_simil_disp,front.ana.phiA_simil_disp/max(abs(phiAcurrent)),			...
					'--','color',[0.2,1.0,0.2],'linewidth',1);
					npl=npl+1;
					plot(axl,front.ana.x_simil_disp,front.ana.phiB_simil_disp/max(abs(phiBcurrent)),			...
					'--','color',[0.2,0.2,1.0],'linewidth',1);
					npl=npl+1;
				end
			end				
		end
%		plotting on linear axis
		if (soln1D.geomdom.isradial ~= 0)
			npllin = 0;
%			plotting advective element location
			plot(axllin,front.r_adv(it)*[1,1],1.25*[-1,1],'-','color',[0.35,0.35,0.35],'linewidth',3); npllin=npllin+1;
%			plotting analytical front locations
			if (postproc.obtainanalytical == 1)
				plot(axllin,front.ana.r_f(it)*[1,1],1.25*[-1,1],'-.','color',[0.35,0.35,0.35],'linewidth',3); npllin=npllin+1;
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					plot(axllin,front.ana.r_f_disp(it)*[1,1],1.25*[-1,1],'--','color',[0.35,0.35,0.35],'linewidth',3); 
					npllin=npllin+1;
				end
			end
%			filling front widths
			if (postproc.plotfrontwidth == 1)
%				filling analytical front widths
				if (postproc.obtainanalytical == 1)
					fill(axllin,	...
					[max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f(it)),	...
					min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f(it)),		...
					min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f(it)),		...
					max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f(it))],	...
					[ylims(1),ylims(1),ylims(2),ylims(2)],	...
					[0.35,0.35,0.35],'LineStyle','-.','LineWidth',1.5,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); 
					npllin=npllin+1;
					fill(axllin,	...
					[max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f_approx(it)),	...
					min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f_approx(it)),	...
					min(xlimr,front.ana.r_f(it)+0.5*front.ana.w_f_approx(it)),	...
					max(xliml,front.ana.r_f(it)-0.5*front.ana.w_f_approx(it))],	...
					[ylims(1),ylims(1),ylims(2),ylims(2)],	...
					[0.35,0.35,0.35],'LineStyle','-.','LineWidth',2,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); 
					npllin=npllin+1;	
					if ((soln1D.fltr.eta ~= 0) && ...
					((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && ...
					(soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
						fill(axllin,	...
						[max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp(it)),	...
						min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp(it)),	...
						min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp(it)),	...
						max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp(it))],	...
						[ylims(1),ylims(1),ylims(2),ylims(2)],	...
						[0.35,0.35,0.35],'LineStyle','--','LineWidth',1.5,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); 
						npllin=npllin+1;		
						fill(axllin,	...
						[max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp_approx(it)),	...
						min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp_approx(it)),	...
						min(xlimr,front.ana.r_f_disp(it)+0.5*front.ana.w_f_disp_approx(it)),	...
						max(xliml,front.ana.r_f_disp(it)-0.5*front.ana.w_f_disp_approx(it))],	...
						[ylims(1),ylims(1),ylims(2),ylims(2)],	...
						[0.35,0.35,0.35],'LineStyle','--','LineWidth',2,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); 
						npllin=npllin+1;		
					end
					fill(axllin,	...
					[max(xliml,front.num.r_f(ifront.tsaved)-0.5*front.num.w_f(ifront.tsaved)),		...
					min(xlimr,front.num.r_f(ifront.tsaved)+0.5*front.num.w_f(ifront.tsaved)),		...
					min(xlimr,front.num.r_f(ifront.tsaved)+0.5*front.num.w_f(ifront.tsaved)),		...
					max(xliml,front.num.r_f(ifront.tsaved)-0.5*front.num.w_f(ifront.tsaved))],		...
					[ylims(1),ylims(1),ylims(2),ylims(2)],	...
					[0.35,0.35,0.35],'LineStyle','-','LineWidth',0.5,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); 
					npllin=npllin+1;	
				end
				fill(axllin,	...
				[max(xliml,front.num.r_f_orig(ifront.tsaved)-0.5*front.num.w_f_orig(ifront.tsaved,2,1)),	...
				min(xlimr,front.num.r_f_orig(ifront.tsaved)+0.5*front.num.w_f_orig(ifront.tsaved,2,1)),		...
				min(xlimr,front.num.r_f_orig(ifront.tsaved)+0.5*front.num.w_f_orig(ifront.tsaved,2,1)),		...
				max(xliml,front.num.r_f_orig(ifront.tsaved)-0.5*front.num.w_f_orig(ifront.tsaved,2,1))],	...
				[ylims(1),ylims(1),ylims(2),ylims(2)],	...
				[0.35,0.35,0.35],'LineStyle','-','LineWidth',1.0,'EdgeColor',[0.35,0.35,0.35],'FaceAlpha',0.01); 
				npllin=npllin+1;	
			end
			npllin2retain = npllin;
%			plotting numerical species solution
			plot(axllin,soln1D.grd.x,phiCcurrent/max(phiCcurrent),'-','color',[1.0,0.2,0.2],'linewidth',2); 
			npllin=npllin+1;
			plot(axllin,soln1D.grd.x,phiBcurrent/max(1,soln1D.fltr.nondim.gamm),'-','color',[0.2,0.2,1.0],'linewidth',2); 
			npllin=npllin+1;
			plot(axllin,soln1D.grd.x,phiAcurrent/max(1,soln1D.fltr.nondim.gamm),'-','color',[0.2,1.0,0.2],'linewidth',2); 
			npllin=npllin+1;
			if (postproc.isexpected == 1)
				plot(axllin,soln1D.grd.x,phiCexpected/max(phiCcurrent),'--','color',[1.0,0.8,0.8],'linewidth',3); 
				npllin=npllin+1;
				plot(axllin,soln1D.grd.x,phiBexpected/max(1,soln1D.fltr.nondim.gamm),'--','color',[0.8,0.8,1.0],'linewidth',3); 
				npllin=npllin+1;
				plot(axllin,soln1D.grd.x,phiAexpected/max(1,soln1D.fltr.nondim.gamm),'--','color',[0.8,1.0,0.8],'linewidth',3); 
				npllin=npllin+1;
			end
			if (postproc.iscomparison1D == 1)
				plot(axllin,solncompare1D.grd.x,phiCcurrentcomparison1D(it,:)/max(phiCcurrent),	...
				'--','color',[1.0,0.4,0.4],'linewidth',4); 
				npllin=npllin+1;
				plot(axllin,solncompare1D.grd.x,phiBcurrentcomparison1D(it,:)/max(1,soln1D.fltr.nondim.gamm),		...
				'--','color',[0.4,0.4,1.0],'linewidth',4); 
				npllin=npllin+1;
				plot(axllin,solncompare1D.grd.x,phiAcurrentcomparison1D(it,:)/max(1,soln1D.fltr.nondim.gamm),		...
				'--','color',[0.4,1.0,0.4],'linewidth',4); 
				npllin=npllin+1;
			end
			if (postproc.iscomparison2D == 1)
				if (soln1D.geomdom.isradial == 0)
					plot(axllin,solncompare2D.grd.x,	...
					phiCcurrentcomparison2D((solncompare2D.simul.ny-1)/2+1,:)/max(phiCcurrent),	...
					':','color',[1.0,0.4,0.4],'linewidth',4); npllin=npllin+1;
					plot(axllin,solncompare2D.grd.x,	...
					phiBcurrentcomparison2D((solncompare2D.simul.ny-1)/2+1,:)/max(1,soln1D.fltr.nondim.gamm),	...
					':','color',[0.4,0.4,1.0],'linewidth',4); npllin=npllin+1;
					plot(axllin,solncompare2D.grd.x,	...
					phiAcurrentcomparison2D((solncompare2D.simul.ny-1)/2+1,:)/max(1,soln1D.fltr.nondim.gamm),	...
					':','color',[0.4,1.0,0.4],'linewidth',4); npllin=npllin+1;
				else
					plot(axllin,r_2D,phiC_2D/max(phiCcurrent),'.','color',[1.0,0.4,0.4],'linewidth',1,'markersize',1.25); 
					npllin=npllin+1;
					plot(axllin,r_2D,phiB_2D/max(1,soln1D.fltr.nondim.gamm),	...
					'.','color',[0.4,0.4,1.0],'linewidth',1,'markersize',1.25);
					npllin=npllin+1;
					plot(axllin,r_2D,phiA_2D/max(1,soln1D.fltr.nondim.gamm),	...
					'.','color',[0.4,1.0,0.4],'linewidth',1,'markersize',1.25);
					npllin=npllin+1;
				end
			end
%			plotting analytical species solution
			if (postproc.obtainanalytical == 1)
				if (soln1D.grd.t(it) < 1.0)
					plot(axllin,front.ana.x_simil_ET,front.ana.phiA_simil_ET/max(abs(phiAcurrent)),			...
					'-.','color',[0.2,1.0,0.2],'linewidth',1);
					npllin=npllin+1;
					plot(axllin,front.ana.x_simil_ET,front.ana.phiB_simil_ET/max(abs(phiBcurrent)),			...
					'-.','color',[0.2,0.2,1.0],'linewidth',1);
					npllin=npllin+1;
					if ((soln1D.fltr.eta ~= 0) && ...
					((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && ...
					(soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
						plot(axllin,front.ana.x_simil_disp_ET,front.ana.phiA_simil_disp_ET/max(abs(phiAcurrent)),	...
						'--','color',[0.2,1.0,0.2],'linewidth',1);
						npllin=npllin+1;
						plot(axllin,front.ana.x_simil_disp_ET,front.ana.phiB_simil_disp_ET/max(abs(phiBcurrent)),	...
						'--','color',[0.2,0.2,1.0],'linewidth',1);
						npllin=npllin+1;
					end
				else
					plot(axllin,front.ana.x_simil,front.ana.phiA_simil/max(abs(phiAcurrent)),		...
					'-.','color',[0.2,1.0,0.2],'linewidth',1);
					npllin=npllin+1;
					plot(axllin,front.ana.x_simil,front.ana.phiB_simil/max(abs(phiBcurrent)),		...
					'-.','color',[0.2,0.2,1.0],'linewidth',1);
					npllin=npllin+1;
					if ((soln1D.fltr.eta ~= 0) && ...
					((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && ...
					(soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
						plot(axllin,front.ana.x_simil_disp,front.ana.phiA_simil_disp/max(abs(phiAcurrent)),	...
						'--','color',[0.2,1.0,0.2],'linewidth',1);
						npllin=npllin+1;
						plot(axllin,front.ana.x_simil_disp,front.ana.phiB_simil_disp/max(abs(phiBcurrent)),	...
						'--','color',[0.2,0.2,1.0],'linewidth',1);
						npllin=npllin+1;
					end
				end					
			end
		end
		pause(postproc.timepausevid);
		movframecurr =     getframe(figl);
		writeVideo(phimovvid,movframecurr);
		if ((postproc.iTransportSS == 1) && (sum(soln1D.simul.tstepsave(itrangeSS)==it)~=0))
			saveas(figl,[postproc.folderloc1D,'/phiSS/phi_it_',num2str(it)],'fig');
		end
		for ii = 1:npl-npl2retain
			lastline = get(axl,'children');
			delete(lastline(1));
		end
		if (soln1D.geomdom.isradial ~= 0)
			for ii = 1:npllin-npllin2retain
				lastline = get(axllin,'children');
				delete(lastline(1));
			end
		end
%		resetting plotline counter
		npl = npl2retain;
		if (soln1D.geomdom.isradial ~= 0)
			npllin = npllin2retain;
		end
		if (postproc.obtainanalytical == 1)
%			plotting conservative species and reaction terms
			plot(axl,soln1D.grd.x,front.num.theta/max(abs(front.num.theta)),'-','color',[0.7,0.7,0.2],'linewidth',2);
			npl = npl+1;
			if (soln1D.geomdom.isradial == 1)
				plot(axl,soln1D.grd.x,front.ana.theta/max(abs(front.num.theta)),'--','color',[0.7,0.7,0.2],'linewidth',3);
				npl = npl+1;
				plot(axl,soln1D.grd.x,front.ana.theta_disp/max(abs(front.num.theta)),':','color',[0.7,0.7,0.2],'linewidth',2);
				npl = npl+1;
			end
			plot(axl,soln1D.grd.x,((front.num.R)/max(abs(front.num.R))),'-','color',[0.7,0.2,0.7],'linewidth',2);
			npl = npl+1;
			if (soln1D.grd.t(it) < 1.0)
				plot(axl,front.ana.x_simil_ET,front.ana.theta_simil_ET/max(abs(front.num.theta)),	...
				'-.','color',[0.7,0.7,0.2],'linewidth',1);
				npl = npl+1;
				plot(axl,front.ana.x_simil_ET,front.ana.R_simil_ET/max(abs(front.num.R)),'-.','color',[0.7,0.2,0.7],'linewidth',1);
				npl = npl+1;
			else
				plot(axl,front.ana.x_simil,front.ana.theta_simil/max(abs(front.num.theta)),			...
				'-.','color',[0.7,0.7,0.2],'linewidth',1);
				npl = npl+1;
				plot(axl,front.ana.x_simil,front.ana.R_simil/max(abs(front.num.R)),'-.','color',[0.7,0.2,0.7],'linewidth',1);
				npl = npl+1;
			end
			if ((soln1D.fltr.eta ~= 0) && ...
			((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
				if (soln1D.grd.t(it) < 1.0)
					plot(axl,front.ana.x_simil_disp_ET,front.ana.theta_simil_disp_ET/max(abs(front.num.theta)),		...
					'--','color',[0.7,0.7,0.2],'linewidth',1);
					npl = npl+1;
					plot(axl,front.ana.x_simil_disp_ET,front.ana.R_simil_disp_ET/max(abs(front.num.R)),				...
						'--','color',[0.7,0.2,0.7],'linewidth',1);
					npl = npl+1;
				else
					plot(axl,front.ana.x_simil_disp,front.ana.theta_simil_disp/max(abs(front.num.theta)),	...
					'--','color',[0.7,0.7,0.2],'linewidth',1);
					npl = npl+1;
					plot(axl,front.ana.x_simil_disp,front.ana.R_simil_disp/max(abs(front.num.R)),			...
						'--','color',[0.7,0.2,0.7],'linewidth',1);
					npl = npl+1;
				end
			end
%			plotting front and maximum reaction rates
			if (postproc.plotmaxfrontR == 1)
				if (soln1D.grd.t(it) < 1.0)
					plot(axl,front.ana.r_f(it),front.ana.R_f_ET(it)/max(abs(front.num.R)),		...
					's','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
					npl = npl+1;
					plot(axl,front.ana.r_f(it),front.ana.R_f_approx_ET(it)/max(abs(front.num.R)),						...
					'o','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
					npl = npl+1;
					plot(axl,front.ana.x_simil(front.ana.R_max_loc_ET),front.ana.R_max_ET(it)/max(abs(front.num.R)),	...
					'v','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
					npl = npl+1;
					if ((soln1D.fltr.eta ~= 0) && ...
					((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe>=1/postproc.thresPedispsimil)))
						plot(axl,front.ana.r_f_disp(it),front.ana.R_f_disp_ET(it)/max(abs(front.num.R)),				...
						'p','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npl = npl+1;
						plot(axl,front.ana.r_f_disp(it),front.ana.R_f_disp_approx_ET(it)/max(abs(front.num.R)),			...
						'h','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npl = npl+1;
						plot(axl,front.ana.x_simil_disp(front.ana.R_max_loc_disp_ET),	...
						front.ana.R_max_disp_ET(it)/max(abs(front.num.R)),				...
						'<','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npl = npl+1;
					end
				else
					plot(axl,front.ana.r_f(it),front.ana.R_f(it)/max(abs(front.num.R)),		...
					's','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
					npl = npl+1;
					plot(axl,front.ana.r_f(it),front.ana.R_f_approx(it)/max(abs(front.num.R)),		...
					'o','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
					npl = npl+1;
					plot(axl,front.ana.x_simil(front.ana.R_max_loc),front.ana.R_max(it)/max(abs(front.num.R)),	...
					'v','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
					npl = npl+1;
					if ((soln1D.fltr.eta ~= 0) && ...
					((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe>=1/postproc.thresPedispsimil)))
						plot(axl,front.ana.r_f_disp(it),front.ana.R_f_disp(it)/max(abs(front.num.R)),	...
						'p','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npl = npl+1;
						plot(axl,front.ana.r_f_disp(it),front.ana.R_f_disp_approx(it)/max(abs(front.num.R)),	...
						'h','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npl = npl+1;
						plot(axl,front.ana.x_simil_disp(front.ana.R_max_loc_disp),			...
						front.ana.R_max_disp(it)/max(abs(front.num.R)),						...
						'<','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npl = npl+1;
					end
				end
			end
			if (soln1D.geomdom.isradial ~= 0)
%				plotting conservative species and reaction terms
				plot(axllin,soln1D.grd.x,front.num.theta/max(abs(front.num.theta)),'-','color',[0.7,0.7,0.2],'linewidth',2);
				npllin = npllin+1;
				if (soln1D.geomdom.isradial == 1)
					plot(axllin,soln1D.grd.x,front.ana.theta/max(abs(front.num.theta)),'--','color',[0.7,0.7,0.2],'linewidth',3);
					npllin = npllin+1;
					plot(axllin,soln1D.grd.x,front.ana.theta_disp/max(abs(front.num.theta)),	...
					':','color',[0.7,0.7,0.2],'linewidth',2);
					npllin = npllin+1;
				end
				plot(axllin,soln1D.grd.x,((front.num.R)/max(abs(front.num.R))),	...
				'-','color',[0.7,0.2,0.7],'linewidth',2);
				npllin = npllin+1;
				plot(axllin,front.ana.x_simil,front.ana.theta_simil/max(abs(front.num.theta)),			...
				'-.','color',[0.7,0.7,0.2],'linewidth',1);
				npllin = npllin+1;
				if (soln1D.grd.t(it) < 1.0)
					plot(axllin,front.ana.x_simil_ET,front.ana.R_simil_ET/max(abs(front.num.R)),		...
					'-.','color',[0.7,0.2,0.7],'linewidth',1);
					npllin = npllin+1;
				else
					plot(axllin,front.ana.x_simil,front.ana.R_simil/max(abs(front.num.R)),				...
					'-.','color',[0.7,0.2,0.7],'linewidth',1);
					npllin = npllin+1;
				end
				if ((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					plot(axllin,front.ana.x_simil_disp,front.ana.theta_simil_disp/max(abs(front.num.theta)),	...
					'--','color',[0.7,0.7,0.2],'linewidth',1);
					npllin = npllin+1;
					if (soln1D.grd.t(it) < 1.0)
						plot(axllin,front.ana.x_simil_disp_ET,front.ana.R_simil_disp_ET/max(abs(front.num.R)),	...
						'--','color',[0.7,0.2,0.7],'linewidth',1);	
						npllin = npllin+1;
					else
						plot(axllin,front.ana.x_simil_disp,front.ana.R_simil_disp/max(abs(front.num.R)),		...	
						'--','color',[0.7,0.2,0.7],'linewidth',1);	
						npllin = npllin+1;
					end
				end
%				plotting front and maximum reaction rates
				if (postproc.plotmaxfrontR == 1)
					if (soln1D.grd.t(it) < 1.0)
						plot(axllin,front.ana.r_f(it),front.ana.R_f_ET(it)/max(abs(front.num.R)),		...
						's','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npllin = npllin+1;
						plot(axl,front.ana.r_f(it),front.ana.R_f_approx_ET(it)/max(abs(front.num.R)),		...
						'o','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npllin = npllin+1;
						plot(axllin,front.ana.x_simil(front.ana.R_max_loc_ET),front.ana.R_max_ET(it)/max(abs(front.num.R)),	...
						'v','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npllin = npllin+1;
						if ((soln1D.fltr.eta ~= 0) && ...
						((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && ...
						(soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
							plot(axllin,front.ana.r_f_disp(it),front.ana.R_f_disp_ET(it)/max(abs(front.num.R)),		...
							'p','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
							npllin = npllin+1;
							plot(axllin,front.ana.r_f_disp(it),front.ana.R_f_disp_approx_ET(it)/max(abs(front.num.R)),	...
							'h','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
							npllin = npllin+1;
							plot(axllin,front.ana.x_simil_disp(front.ana.R_max_loc_disp_ET),	...
								front.ana.R_max_disp_ET(it)/max(abs(front.num.R)),...
							'<','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
							npllin = npllin+1;
						end
					else
						plot(axllin,front.ana.r_f(it),front.ana.R_f(it)/max(abs(front.num.R)),			...
						's','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npllin = npllin+1;
						plot(axl,front.ana.r_f(it),front.ana.R_f_approx(it)/max(abs(front.num.R)),		...
						'o','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npllin = npllin+1;
						plot(axllin,front.ana.x_simil(front.ana.R_max_loc),front.ana.R_max(it)/max(abs(front.num.R)),	...
						'v','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
						npllin = npllin+1;
						if ((soln1D.fltr.eta ~= 0) && ...
							((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && ...
							(soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
							plot(axllin,front.ana.r_f_disp(it),front.ana.R_f_disp(it)/max(abs(front.num.R)),		...
							'p','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
							npllin = npllin+1;
							plot(axllin,front.ana.r_f_disp(it),front.ana.R_f_disp_approx(it)/max(abs(front.num.R)),	...
							'h','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
							npllin = npllin+1;
							plot(axllin,front.ana.x_simil_disp(front.ana.R_max_loc_disp),							...
							front.ana.R_max_disp(it)/max(abs(front.num.R)),...
							'<','color',[1.0,0.2,1.0],'linewidth',2,'markersize',3);
							npllin = npllin+1;
						end						
					end
				end
			end
			pause(postproc.timepausevid);
			movframecurr =     getframe(figl);
			writeVideo(phimovvidA,movframecurr);
		if ((postproc.iTransportSS == 1) && (sum(soln1D.simul.tstepsave(itrangeSS)==it)~=0))
			saveas(figl,[postproc.folderloc1D,'/phiSS/anasoln_it_',num2str(it)],'fig');
		end
		end
		for ii = 1:npl
			lastline = get(axl,'children');
			delete(lastline(1));
		end
		if (soln1D.geomdom.isradial ~= 0)
			for ii = 1:npllin
				lastline = get(axllin,'children');
				delete(lastline(1));
			end
		end
% 	catch
% 	end
end
close(figl);
close(phimovvid);
close(phimovvidA);
