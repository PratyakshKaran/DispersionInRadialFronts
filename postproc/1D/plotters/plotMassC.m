% plotting mass of C

figl  =	figure('position',[100,100,750,500],'visible',postproc.vistog); axl = axes(figl,'position',[0.20,0.20,0.75,0.75]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
set(axl,'xscale','log'); set(axl,'yscale','log');
xlabel('$t$','interpreter','latex');
if (soln1D.geomdom.isradial == 0)
	ylblstring =	'$\displaystyle {M_{\rm C}}/{L_yL_z} = m_{\rm C} = \int_{x=0}^{x=L_x}\phi_C dx$';
elseif (soln1D.geomdom.isradial == 1)
	ylblstring =	'$\displaystyle {M_{\rm C}}/{L_z} = m_{\rm C} = 2\pi\int_{r=a}^{r=A}\phi_C rdr$';
else
	ylblstring =	'$\displaystyle M_{\rm C} = 4\pi\int_{r=a}^{r=A}\phi_C r^2dr$';
end
ylabel(ylblstring,'interpreter','latex');
title(postproc.casename,'Interpreter','none','FontSize',12);
if ((soln1D.simul.compMassC == 0) || ((soln1D.simul.compMassC == 1) && (postproc.forcemassC == 1)))
	plot(axl,front.tsaved,soln_analyze.MassC,'-','color',[0.0,0.0,1.0],'linewidth',2.5);
	plot(axl,front.tsaved,soln_analyze.MassC_4mR,'-.','color',[0.0,0.0,1.0],'linewidth',2.5);
	plot(axl,front.tsaved,soln_analyze.MassC_4mR_twk,'--','color',[0.0,0.0,1.0],'linewidth',2.5);
	plot(axl,front.tsaved,soln_analyze.barR,'-','color',[0.0,1.0,0.0],'linewidth',2.5);
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,2):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,2));
				coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
				patchline(front.tsaved(tspread),soln_analyze.MassC(tspread),				...
				'linestyle','-','linewidth',5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(front.tsaved(tspread),soln_analyze.MassC_4mR(tspread),			...
				'linestyle','-.','linewidth',5, ...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(front.tsaved(tspread),soln_analyze.MassC_4mR_twk(tspread),		...
				'linestyle','--','linewidth',5, ...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(front.tsaved(tspread),soln_analyze.barR(tspread),					...
				'linestyle','-','linewidth',5,	...
				'edgecolor',[0.0,1.0,0.0],'facecolor',[0.0,1.0,0.0],'edgealpha',coltr,'facealpha',coltr);
			end
		end
	end
end
if (soln1D.simul.compMassC == 1)
	plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC,'-','color',[0.0,0.0,1.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC_4mR,'-.','color',[0.0,0.0,1.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC_4mR_twk,'--','color',[0.0,0.0,1.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.barR,'-','color',[0.0,1.0,0.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.dMassCdt,'--','color',[0.0,1.0,0.0],'linewidth',1.25);
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,1):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,1));
				coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.MassC(tspread),				...
				'linestyle','-','linewidth',2.5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);				
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.MassC_4mR(tspread),			...
				'linestyle','-.','linewidth',2.5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.MassC_4mR_twk(tspread),		...
				'linestyle','--','linewidth',2.5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.barR(tspread),				...
				'linestyle','-','linewidth',2.5,	...
				'edgecolor',[0.0,1.0,0.0],'facecolor',[0.0,1.0,0.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.dMassCdt(tspread),			...
				'linestyle','--','linewidth',2.5,	...
				'edgecolor',[0.0,1.0,0.0],'facecolor',[0.0,1.0,0.0],'edgealpha',coltr,'facealpha',coltr);				
			end
		end
	end
end
if (postproc.obtainanalytical == 1)
	range_ET = 1:sum(soln1D.grd.t<1.0);
	range_LT = range_ET(end)+1:length(soln1D.grd.t);	
	plot(axl,soln1D.grd.t(range_LT),front.ana.Mc(range_LT),'-.','color',[0.0,0.2,0.8],'linewidth',1.25);
	plot(axl,soln1D.grd.t(range_LT),front.ana.dMcdt(range_LT),'-.','color',[0.2,0.8,0.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t(range_ET),front.ana.Mc_ET(range_ET),'-.','color',[0.0,0.2,0.8],'linewidth',1.25);
	plot(axl,soln1D.grd.t(range_ET),front.ana.dMcdt_ET(range_ET),'-.','color',[0.2,0.8,0.0],'linewidth',1.25);
	if ((soln1D.geomdom.isradial == 1) && ...
		((soln1D.fltr.eta ~= 0) && ...
		((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil))))
		plot(axl,soln1D.grd.t(range_LT),front.ana.Mc_disp(range_LT),':','color',[0.0,0.2,0.8],'linewidth',1.25);
		plot(axl,soln1D.grd.t(range_LT),front.ana.dMcdt_disp(range_LT),':','color',[0.2,0.8,0.0],'linewidth',1.25);
		plot(axl,soln1D.grd.t(range_ET),front.ana.Mc_disp_ET(range_ET),':','color',[0.0,0.2,0.8],'linewidth',1.25);
		plot(axl,soln1D.grd.t(range_ET),front.ana.dMcdt_disp_ET(range_ET),':','color',[0.2,0.8,0.0],'linewidth',1.25);
	end
end
legend({'$\displaystyle M_C$','$\displaystyle \int\limits_{\tau=0}^{\tau=t} \bar{R}~d\tau$',							...
		'$\displaystyle \int\limits_{\tau=0}^{\tau=t}\left(\bar{R}2\pi\phi_{C(r\rightarrow 0)}\right)~d\tau$',			...
		'$\displaystyle \bar{R} = \int\limits_{r\rightarrow 0}^{r\rightarrow\infty} \phi_C~dA$'},						...
		'interpreter','latex','fontsize',12,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/MassC'],'fig');
close(figl);

if (soln1D.geomdom.isradial == 2)

	figl  =	figure('position',[100,100,750,500],'visible',postproc.vistog); axl = axes(figl,'position',[0.20,0.20,0.75,0.75]);
	set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
	set(axl,'xscale','log'); set(axl,'yscale','log');
	xlabel('$t$','interpreter','latex');
	ylblstring =	'$\displaystyle M_{\rm C}/\eta^{\frac{3}{4}} = (4\pi\int_{r=a}^{r=A}\phi_C r^2dr) /\eta^{\frac{3}{4}}$';
	ylabel(ylblstring,'interpreter','latex');
	title(postproc.casename,'Interpreter','none','FontSize',12);
	if ((soln1D.simul.compMassC == 0) || ((soln1D.simul.compMassC == 1) && (postproc.forcemassC == 1)))
		plot(axl,front.tsaved,soln_analyze.MassC/(soln1D.fltr.eta^(3/4)),'-','color',[0.0,0.0,1.0],'linewidth',2.5);
		plot(axl,front.tsaved,soln_analyze.MassC_4mR/(soln1D.fltr.eta^(3/4)),'-.','color',[0.0,0.0,1.0],'linewidth',2.5);
		plot(axl,front.tsaved,soln_analyze.MassC_4mR_twk/(soln1D.fltr.eta^(3/4)),'--','color',[0.0,0.0,1.0],'linewidth',2.5);
		plot(axl,front.tsaved,soln_analyze.barR/(soln1D.fltr.eta^(3/4)),'-','color',[0.0,1.0,0.0],'linewidth',2.5);
		for isolvmode = 1:2
			for ilimdef = 1:2
				for ifracTransitn = 1:length(postproc.fracTransitn)
					tspread = ...
					(front.trnszon(isolvmode,ifracTransitn,ilimdef,1,2):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,2));
					coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
					patchline(front.tsaved(tspread),soln_analyze.MassC(tspread)/(soln1D.fltr.eta^(3/4)),				...
					'linestyle','-','linewidth',5,	...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
					patchline(front.tsaved(tspread),soln_analyze.MassC_4mR(tspread)/(soln1D.fltr.eta^(3/4)),			...
					'linestyle','-.','linewidth',5, ...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
					patchline(front.tsaved(tspread),soln_analyze.MassC_4mR_twk(tspread)/(soln1D.fltr.eta^(3/4)),		...
					'linestyle','--','linewidth',5, ...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
					patchline(front.tsaved(tspread),soln_analyze.barR(tspread)/(soln1D.fltr.eta^(3/4)),					...
					'linestyle','-','linewidth',5,	...
					'edgecolor',[0.0,1.0,0.0],'facecolor',[0.0,1.0,0.0],'edgealpha',coltr,'facealpha',coltr);
				end
			end
		end
	end
	if (soln1D.simul.compMassC == 1)
		plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC/(soln1D.fltr.eta^(3/4)),'-','color',[0.0,0.0,1.0],'linewidth',1.25);
		plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC_4mR/(soln1D.fltr.eta^(3/4)),'-.','color',[0.0,0.0,1.0],'linewidth',1.25);
		plot(axl,soln1D.grd.t,soln1D.solnreactive.MassC_4mR_twk/(soln1D.fltr.eta^(3/4)),'--','color',[0,0.0,1.0],'linewidth',1.25);
		plot(axl,soln1D.grd.t,soln1D.solnreactive.barR/(soln1D.fltr.eta^(3/4)),'-','color',[0.0,1.0,0.0],'linewidth',1.25);
		plot(axl,soln1D.grd.t,soln1D.solnreactive.dMassCdt/(soln1D.fltr.eta^(3/4)),'--','color',[0.0,1.0,0.0],'linewidth',1.25);
		for isolvmode = 1:2
			for ilimdef = 1:2
				for ifracTransitn = 1:length(postproc.fracTransitn)
					tspread = ...
					(front.trnszon(isolvmode,ifracTransitn,ilimdef,1,1):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,1));
					coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
					patchline(soln1D.grd.t(tspread),soln1D.solnreactive.MassC(tspread)/(soln1D.fltr.eta^(3/4)),				...
					'linestyle','-','linewidth',2.5,	...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);				
					patchline(soln1D.grd.t(tspread),soln1D.solnreactive.MassC_4mR(tspread)/(soln1D.fltr.eta^(3/4)),			...
					'linestyle','-.','linewidth',2.5,	...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
					patchline(soln1D.grd.t(tspread),soln1D.solnreactive.MassC_4mR_twk(tspread)/(soln1D.fltr.eta^(3/4)),		...
					'linestyle','--','linewidth',2.5,	...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
					patchline(soln1D.grd.t(tspread),soln1D.solnreactive.barR(tspread)/(soln1D.fltr.eta^(3/4)),				...
					'linestyle','-','linewidth',2.5,	...
					'edgecolor',[0.0,1.0,0.0],'facecolor',[0.0,1.0,0.0],'edgealpha',coltr,'facealpha',coltr);
					patchline(soln1D.grd.t(tspread),soln1D.solnreactive.dMassCdt(tspread)/(soln1D.fltr.eta^(3/4)),			...
					'linestyle','--','linewidth',2.5,	...
					'edgecolor',[0.0,1.0,0.0],'facecolor',[0.0,1.0,0.0],'edgealpha',coltr,'facealpha',coltr);				
				end
			end
		end
	end
	if (postproc.obtainanalytical == 1)
		range_ET = 1:sum(soln1D.grd.t<1.0);
		range_LT = range_ET(end)+1:length(soln1D.grd.t);	
		plot(axl,soln1D.grd.t(range_LT),front.ana.Mc(range_LT)/(soln1D.fltr.eta^(3/4)), ...
			'-.','color',[0.0,0.2,0.8],'linewidth',1.25);
		plot(axl,soln1D.grd.t(range_LT),front.ana.dMcdt(range_LT)/(soln1D.fltr.eta^(3/4)), ...
			'-.','color',[0.2,0.8,0.0],'linewidth',1.25);
		plot(axl,soln1D.grd.t(range_ET),front.ana.Mc_ET(range_ET)/(soln1D.fltr.eta^(3/4)), ...
			'-.','color',[0.0,0.2,0.8],'linewidth',1.25);
		plot(axl,soln1D.grd.t(range_ET),front.ana.dMcdt_ET(range_ET)/(soln1D.fltr.eta^(3/4)), ...
			'-.','color',[0.2,0.8,0.0],'linewidth',1.25);
		if ((soln1D.geomdom.isradial == 1) && ...
			((soln1D.fltr.eta ~= 0) && ...
			((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil))))
			plot(axl,soln1D.grd.t(range_LT),front.ana.Mc_disp(range_LT)/(soln1D.fltr.eta^(3/4)), ...
				':','color',[0.0,0.2,0.8],'linewidth',1.25);
			plot(axl,soln1D.grd.t(range_LT),front.ana.dMcdt_disp(range_LT)/(soln1D.fltr.eta^(3/4)), ...
				':','color',[0.2,0.8,0.0],'linewidth',1.25);
			plot(axl,soln1D.grd.t(range_ET),front.ana.Mc_disp_ET(range_ET)/(soln1D.fltr.eta^(3/4)), ...
				':','color',[0.0,0.2,0.8],'linewidth',1.25);
			plot(axl,soln1D.grd.t(range_ET),front.ana.dMcdt_disp_ET(range_ET)/(soln1D.fltr.eta^(3/4)), ...
				':','color',[0.2,0.8,0.0],'linewidth',1.25);
		end
	end
	legend({'$\displaystyle M_C$','$\displaystyle \int\limits_{\tau=0}^{\tau=t} \bar{R}~d\tau$',							...
			'$\displaystyle \int\limits_{\tau=0}^{\tau=t}\left(\bar{R}2\pi\phi_{C(r\rightarrow 0)}\right)~d\tau$',			...
			'$\displaystyle \bar{R} = \int\limits_{r\rightarrow 0}^{r\rightarrow\infty} \phi_C~dA$'},						...
			'interpreter','latex','fontsize',12,'location','eastoutside');
	saveas(figl,[postproc.folderloc1D,'/MassC_by_etapow0d75'],'fig');
	close(figl);

end