% plotting front location
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
set(axl,'xscale','log'); set(axl,'yscale','log');
xlabel('$t$','interpreter','latex');
ylabel('$x_f-x_0$','interpreter','latex');
title(postproc.casename,'Interpreter','none','FontSize',12);
plot(axl,soln1D.grd.t,front.r_adv-front.r_adv(1),'-','color',[0.0,0.0,0.0],'linewidth',0.5);
plot(axl,front.tsaved,front.num.r_f-front.num.r_f(1),'--','color',[0.0,0.0,1.0],'linewidth',2.5);
hold on;
plot(axl,front.tsaved,front.num.r_f_orig-front.num.r_f_orig(1),'-','color',[0.0,0.0,1.0],'linewidth',2.5);
for isolvmode = 1:2
	for ilimdef = 1:2
		for ifracTransitn = 1:length(postproc.fracTransitn)
			tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,2):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,2));
			coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
			patchline(front.tsaved(tspread),front.num.r_f(tspread)-front.num.r_f(1),				...
			'linestyle','--','linewidth',5,	...
			'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
			patchline(front.tsaved(tspread),front.num.r_f_orig(tspread)-front.num.r_f_orig(1),		...
			'linestyle','-','linewidth',5, ...
			'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
		end
	end
end
if (soln1D.simul.compfront == 1)
	plot(axl,soln1D.grd.t,soln1D.solnreactive.xf-soln1D.solnreactive.xf(1),'--','color',[0.0,0.0,1.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.xforig-soln1D.solnreactive.xforig(1),'-','color',[0.0,0.0,1.0],'linewidth',1.25);
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,1):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,1));
				coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.xf(tspread)-soln1D.solnreactive.xf(1),				...
				'linestyle','--','linewidth',5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.xforig(tspread)-soln1D.solnreactive.xforig(1),		...
				'linestyle','-','linewidth',5, ...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
			end
		end
	end
end
if (postproc.obtainanalytical == 1)
	plot(axl,soln1D.grd.t,front.ana.r_f-front.ana.r_f(1),'-','color',[0.0,0.3,0.7],'linewidth',1);
	if ((soln1D.geomdom.isradial == 1) && ...
	((soln1D.fltr.eta ~= 0) && ...
	((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil))))
		plot(axl,soln1D.grd.t,front.ana.r_f_disp-front.ana.r_f_disp(1),'-','color',[0.3,0.7,0.0],'linewidth',1);
	end
end
legend({'advective front','$\max(R)$','$\min(\theta)$'},'interpreter','latex','fontsize',12,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Location_Front'],'fig');
close(figl);

% plotting front reaction rate
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
set(axl,'xscale','log'); set(axl,'yscale','log');
xlabel('$t$','interpreter','latex');
ylabel('$R_f$','interpreter','latex');
title(postproc.casename,'Interpreter','none','FontSize',12);
plot(axl,front.tsaved,front.num.R_f,'--','color',[0.0,0.0,1.0],'linewidth',2.5);
hold on;
plot(axl,front.tsaved,front.num.R_f_orig,'-','color',[0.0,0.0,1.0],'linewidth',2.5);
for isolvmode = 1:2
	for ilimdef = 1:2
		for ifracTransitn = 1:length(postproc.fracTransitn)
			tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,2):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,2));
			coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
			patchline(front.tsaved(tspread),front.num.R_f(tspread),				...
			'linestyle','--','linewidth',5,	...
			'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
			patchline(front.tsaved(tspread),front.num.R_f_orig(tspread),		...
			'linestyle','-','linewidth',5, ...
			'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
		end
	end
end
if (soln1D.simul.compfront == 1)
	plot(axl,soln1D.grd.t,soln1D.solnreactive.Rf,'--','color',[0.0,0.0,1.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.Rforig,'-','color',[0.0,0.0,1.0],'linewidth',1.25);
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,1):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,1));
				coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.Rf(tspread),				...
				'linestyle','--','linewidth',5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.Rforig(tspread),			...
				'linestyle','-','linewidth',5, ...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
			end
		end
	end
end
if (postproc.obtainanalytical == 1)
	range_ET = 1:sum(soln1D.grd.t<1.0);
	range_LT = range_ET(end)+1:length(soln1D.grd.t);	
	plot(axl,soln1D.grd.t(range_LT),front.ana.R_f(range_LT),'-','color',[0.0,0.3,0.7],'linewidth',1);
	plot(axl,soln1D.grd.t(range_LT),front.ana.R_max(range_LT),'--','color',[0.0,0.3,0.7],'linewidth',1);
	plot(axl,soln1D.grd.t(range_ET),front.ana.R_f_ET(range_ET),'-','color',[0.0,0.3,0.7],'linewidth',1);
	plot(axl,soln1D.grd.t(range_ET),front.ana.R_max_ET(range_ET),'--','color',[0.0,0.3,0.7],'linewidth',1);
	if ((soln1D.geomdom.isradial == 1) && ...
	((soln1D.fltr.eta ~= 0) && ...
	((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil))))
		plot(axl,soln1D.grd.t(range_LT),front.ana.R_f_disp(range_LT),'-','color',[0.3,0.7,0.0],'linewidth',1);
		plot(axl,soln1D.grd.t(range_LT),front.ana.R_max_disp(range_LT),'--','color',[0.3,0.7,0.0],'linewidth',1);
		plot(axl,soln1D.grd.t(range_ET),front.ana.R_f_disp_ET(range_ET),'-','color',[0.3,0.7,0.0],'linewidth',1);
		plot(axl,soln1D.grd.t(range_ET),front.ana.R_max_disp_ET(range_ET),'--','color',[0.3,0.7,0.0],'linewidth',1);
	end
end
legend({'at $\max(R)$','at $\min(\theta)$'},'interpreter','latex','fontsize',12,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/ReactionRate_Front'],'fig');
close(figl);

% plotting reaction zone width
linetypewidth = {':','--';'-','-.'};
linethickwidth = {1.5,1.5;2.5,1.5};
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
set(axl,'xscale','log'); set(axl,'yscale','log');
xlabel('$t$','interpreter','latex');
ylabel('$w_f$','interpreter','latex');
title(postproc.casename,'Interpreter','none','FontSize',12);
plot(axl,front.tsaved,front.num.w_f,'--','color',[0.0,0.0,1.0],'linewidth',2.5);
hold on;
for iwidthdef = 1:2
	for iwidthbasefunc = 1:2
		plot(axl,front.tsaved,front.num.w_f_orig(:,iwidthdef,iwidthbasefunc),	...
			linetypewidth{iwidthdef,iwidthbasefunc},'color',[0.0,0.0,1.0],'linewidth',linethickwidth{iwidthdef,iwidthbasefunc});
	end
end
for isolvmode = 1:2
	for ilimdef = 1:2
		for ifracTransitn = 1:length(postproc.fracTransitn)
			tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,2):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,2));
			coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
			patchline(front.tsaved(tspread),front.num.w_f(tspread),				...
			'linestyle','--','linewidth',5,	...
			'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
			for iwidthdef = 1:2
				for iwidthbasefunc = 1:2
					patchline(front.tsaved(tspread),front.num.w_f_orig(tspread,iwidthdef,iwidthbasefunc),		...
					'linestyle',linetypewidth{iwidthdef,iwidthbasefunc},'linewidth',linethickwidth{iwidthdef,iwidthbasefunc}, ...
					'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				end
			end
		end
	end
end
if (soln1D.simul.compfront == 1)
	plot(axl,soln1D.grd.t,soln1D.solnreactive.wf,'--','color',[0.0,0.0,1.0],'linewidth',1.25);
	plot(axl,soln1D.grd.t,soln1D.solnreactive.wforig,'-','color',[0.0,0.0,1.0],'linewidth',1.25);
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				tspread = (front.trnszon(isolvmode,ifracTransitn,ilimdef,1,1):front.trnszon(isolvmode,ifracTransitn,ilimdef,2,1));
				coltr = 0.20*(ifracTransitn/length(postproc.fracTransitn));
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.wf(tspread),				...
				'linestyle','--','linewidth',5,	...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
				patchline(soln1D.grd.t(tspread),soln1D.solnreactive.wforig(tspread),			...
				'linestyle','-','linewidth',5, ...
				'edgecolor',[0.0,0.0,1.0],'facecolor',[0.0,0.0,1.0],'edgealpha',coltr,'facealpha',coltr);
			end
		end
	end
end
if (postproc.obtainanalytical == 1)
	range_ET = 1:sum(soln1D.grd.t<1.0);
	range_LT = range_ET(end)+1:length(soln1D.grd.t);	
	plot(axl,soln1D.grd.t(range_LT),front.ana.w_f(range_LT),'-','color',[0.0,0.3,0.7],'linewidth',1);
	plot(axl,soln1D.grd.t(range_LT),front.ana.w_f_approx(range_LT),'--','color',[0.0,0.3,0.7],'linewidth',1);
	plot(axl,soln1D.grd.t(range_ET),front.ana.w_f_ET(range_ET),'-','color',[0.0,0.3,0.7],'linewidth',1);
	plot(axl,soln1D.grd.t(range_ET),front.ana.w_f_approx_ET(range_ET),'--','color',[0.0,0.3,0.7],'linewidth',1);
	if ((soln1D.geomdom.isradial == 1) && ...
		((soln1D.fltr.eta ~= 0) && ...
		((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil))))
		plot(axl,soln1D.grd.t(range_LT),front.ana.w_f_disp(range_LT),'-','color',[0.3,0.7,0.0]	,'linewidth',1);
		plot(axl,soln1D.grd.t(range_LT),front.ana.w_f_disp_approx(range_LT),'--','color',[0.0,0.3,0.7],'linewidth',1);
		plot(axl,soln1D.grd.t(range_ET),front.ana.w_f_disp_ET(range_ET),'-','color',[0.3,0.7,0.0]	,'linewidth',1);
		plot(axl,soln1D.grd.t(range_ET),front.ana.w_f_disp_approx_ET(range_ET),'--','color',[0.0,0.3,0.7],'linewidth',1);
	end
end
legend({'half concentration','second moment'},'interpreter','latex','fontsize',12,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Width_Front'],'fig');
close(figl);
