% plotting mass of C

figl  =	figure('position',[100,100,750,500],'visible',vistog); axl = axes(figl,'position',[0.20,0.20,0.75,0.75]);
set(axl,'box','on'); grid(gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
set(axl,'xscale','log'); set(axl,'yscale','log');
xlabel('$t$','interpreter','latex');
ylabel('$\displaystyle M_{C}/M_{C(0)}$','interpreter','latex');
title(casename,'Interpreter','none','FontSize',12);
if ((soln1D.simul.compMassC == 0) || ((soln1D.simul.compMassC == 1) && (forceMassCnum == 1)))
	plot(axl,front.tsaved,soln_analyze.MassC,'-','color',[0.2,0.2,1.0],'linewidth',2);
	plot(axl,front.tsaved,soln_analyze.MassC_4mR,'-.','color',[0.2,0.2,1.0],'linewidth',2);
	plot(axl,front.tsaved,soln_analyze.MassC_4mR_twk,'--','color',[0.2,0.2,1.0],'linewidth',2);
	plot(axl,front.tsaved,soln_analyze.barR,'-','color',[0.2,1.0,0.2],'linewidth',2);
end
if (soln1D.simul.compMassC == 1)
	plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC,'-','color',[0.2,0.6,0.6],'linewidth',2);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC_4mR,'-.','color',[0.2,0.6,0.6],'linewidth',2);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC_4mR_twk,'--','color',[0.2,0.6,0.6],'linewidth',2);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.barR,'-','color',[0.5,0.8,0.1],'linewidth',2);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.dMassCdt,'--','color',[0.5,0.8,0.1],'linewidth',2.5);
end
legend({'$\displaystyle M_C$','$\displaystyle \int\limits_{\tau=0}^{\tau=t} \bar{R}~d\tau$',							...
		'$\displaystyle \int\limits_{\tau=0}^{\tau=t}\left(\bar{R}2\pi\phi_{C(r\rightarrow 0)}\right)~d\tau$',			...
		'$\displaystyle \bar{R} = \int\limits_{r\rightarrow 0}^{r\rightarrow\infty} \phi_C~dA$'},						...
		'interpreter','latex','fontsize',12,'location','eastoutside');
saveas(figl,[folderloc1D,'/MassC_Difference',namenormalize{inormalize+1}],'fig');
close(figl);