% plotting flow solution

if (soln1D.geomdom.isradial == 0)
	Q =								-(soln1D.fltr.flowbc.Hright0-soln1D.fltr.flowbc.Hleft0)/soln1D.geomdom.size.Lx;
	H_ana =							soln1D.fltr.flowbc.Hleft0-Q*soln1D.grd.x;
	v_ana =							Q*ones(1,soln1D.simul.nx);
elseif (soln1D.geomdom.isradial == 1)
	Q =								-((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)/...
									(log(soln1D.geomdom.size.a)-log(soln1D.geomdom.size.A)));
	H_ana =							((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)*log(soln1D.grd.x)-	...
									(soln1D.fltr.flowbc.Hleft0*log(soln1D.geomdom.size.A)-	...
									soln1D.fltr.flowbc.Hright0*log(soln1D.geomdom.size.a)))/	...
									(log(soln1D.geomdom.size.a)-log(soln1D.geomdom.size.A));
	v_ana =							Q./soln1D.grd.x;
else
	Q =								(((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)*		...
									soln1D.geomdom.size.A*soln1D.geomdom.size.a)/(soln1D.geomdom.size.A-soln1D.geomdom.size.a));
	H_ana =							soln1D.fltr.flowbc.Hleft0-	...
									(((soln1D.fltr.flowbc.Hleft0-soln1D.fltr.flowbc.Hright0)*		...
									(soln1D.grd.x-soln1D.geomdom.size.a))/	...
									(soln1D.geomdom.size.A-soln1D.geomdom.size.a)).*(soln1D.geomdom.size.A./soln1D.grd.x);
	v_ana =							Q./(soln1D.grd.x.^2);
end
figl  =	figure('position',[100,35,1250,700],'visible',postproc.vistog); axl = axes(figl,'position',[0.125,0.15,0.85,0.800]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
if (soln1D.geomdom.isradial == 0)
	xliml =			0.0;
	xlimr =			soln1D.geomdom.size.Lx;
else
	xliml =			max([soln1D.geomdom.size.a,postproc.maxxrange_log*soln1D.geomdom.size.A]);
	xlimr =			soln1D.geomdom.size.A;
end
lowbound =	min([min(min(soln1D.soln.flow.H)),min(min(soln1D.soln.flow.vel.v))]); 
highbound = max([max(max(soln1D.soln.flow.H)),max(max(soln1D.soln.flow.vel.v))]);
set(axl,'xlim',[xliml,xlimr]);
set(axl,'ylim',[0.85*double(lowbound>0)*lowbound+1.25*double(lowbound<0)*lowbound,	...
				1.25*double(highbound>0)*highbound+0.85*double(highbound<0)*highbound+eps]);
if (soln1D.geomdom.isradial ~= 0)
	set(axl,'xscale','log');
	set(axl,'yscale','log');
end
title(postproc.casename,'Interpreter','none','FontSize',12);
if (soln1D.geomdom.isradial == 0)
	xlabel('$x$','interpreter','latex');
else
	xlabel('$r$','interpreter','latex');
end
ylabel('$H,v$','interpreter','latex');
if (soln1D.simul.givenvel == 0)
	plot(axl,soln1D.grd.x,soln1D.soln.flow.H,'-','color',[0.05,0.35,0.50],'linewidth',1.5);
	plot(axl,soln1D.grd.x,soln1D.soln.flow.vel.v,'--','color',[0.2,0.2,1.0],'linewidth',1.5);
	plot(axl,soln1D.grd.x,H_ana,'c-','linewidth',0.2);
	plot(axl,soln1D.grd.x,v_ana,'c--','linewidth',0.2);
	plot(axl,soln1D.grd.x,soln1D.soln.flow.vel.v.*(soln1D.grd.x'.^soln1D.geomdom.isradial),...
		'-.','color',[0.2,0.2,1.0],'linewidth',1.5);
	plot(axl,soln1D.grd.x,v_ana.*(soln1D.grd.x.^soln1D.geomdom.isradial),'c-.','linewidth',0.2);
else
	plot(axl,soln1D.grd.x,soln1D.soln.flow.H,'-','color',[0.05,0.35,0.50],'linewidth',1.5);
	plot(axl,soln1D.grd.x,soln1D.soln.flow.vel.v','--','color',[0.2,0.2,1.0],'linewidth',1.5);
	plot(axl,soln1D.grd.x,H_ana,'c-','linewidth',0.2);
	plot(axl,soln1D.grd.x,v_ana,'c--','linewidth',0.2);
	plot(axl,soln1D.grd.x,soln1D.soln.flow.vel.v'.*(soln1D.grd.x'.^soln1D.geomdom.isradial),...
		'-.','color',[0.2,0.2,1.0],'linewidth',1.5);
	plot(axl,soln1D.grd.x,v_ana.*(soln1D.grd.x.^soln1D.geomdom.isradial),'c-.','linewidth',0.2);
end
lgnd = {'$H$','$v$'};
if (soln1D.geomdom.isradial == 0)
	lgnd{3} =	'$\displaystyle H_I-\frac{(H_I-H_O)x}{L_x}$ ';
	lgnd{4} =	'$\displaystyle \frac{H_I-H_O}{L_x}$';
	lgnd{5} =	'$\displaystyle v_x$';
	lgnd{6} =	'$\displaystyle v_{x \rm{(ana)}}$';
elseif (soln1D.geomdom.isradial == 1)
	lgnd{3} =	'$\displaystyle \frac{(H_I-H_O)\log(r)-(H_I\log(A)-H_O\log(a))}{\log(a)-\log(A)}$ ';
	lgnd{4} =	'$\displaystyle -\frac{(H_I-H_O)}{r}$';
	lgnd{5} =	'$\displaystyle vr$';
	lgnd{6} =	'$\displaystyle v_{\rm{(ana)}}r$';
else
	lgnd{3} =	'$\displaystyle H_I-\frac{(H_I-H_O)(r-a)}{A-a}\frac{A}{r}$ ';
	lgnd{4} =	'$\displaystyle -\frac{(H_I-H_O)aA}{(A-a)r^2}$';
	lgnd{5} =	'$\displaystyle vr^2$';
	lgnd{6} =	'$\displaystyle v_{\rm{(ana)}}r^2$';
end
legend(lgnd,'interpreter','latex','fontsize',12,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Flow'],'fig');
close(figl);