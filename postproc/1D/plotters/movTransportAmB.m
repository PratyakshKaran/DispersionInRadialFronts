% generating AmB movie

phimovvid = VideoWriter([postproc.folderloc1D,'/Transport_A_minus_B_Equation.avi']);
open(phimovvid);
figl  =	figure('position',[100,35,1250,700],'visible',postproc.vistogmov);
if (soln1D.geomdom.isradial ~= 0)
	axl = axes(figl,'position',[0.125,0.125,0.80,0.325]); hold on;
	if (soln1D.geomdom.isradial ~= 0)
		axllin = axes(figl,'position',[0.125,0.600,0.80,0.35]); hold on;
	end
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
end
set(axl,'xlim',[xliml,xlimr]);
if (soln1D.geomdom.isradial ~= 0)
	set(axl,'xscale','log');
end
if (soln1D.geomdom.isradial == 0)
	xlabel('$x$','interpreter','latex');
else
	xlabel('$r$','interpreter','latex');
end
ylabel('Equation($\theta = \phi_A-\phi_B$)','interpreter','latex');
if (soln1D.geomdom.isradial == 0)
	LHSstring =	'LHS = $\displaystyle \frac{\partial \theta}{\partial t} + v\frac{\partial \theta}{\partial x}$';
	RHSstring = 'RHS = $\displaystyle \frac{1}{Pe}\left[(\eta Pe v+c_{df})\frac{\partial^2 \theta}{\partial x^2}\right]$';
elseif (soln1D.geomdom.isradial == 1)
	LHSstring =	'LHS = $\displaystyle \frac{\partial \theta}{\partial t} + v\frac{\partial \theta}{\partial r}$';
	RHSstring = 'RHS';
else
	LHSstring =	'LHS = $\displaystyle \frac{\partial \theta}{\partial t} + v\frac{\partial \theta}{\partial r}$';
	RHSstring = 'RHS';
end
itprevprev =	soln1D.simul.tstepsave(1);
itprev =		soln1D.simul.tstepsave(2);
for it = soln1D.simul.tstepsave(3:end)
	try
		title([postproc.casename,': time is ',num2str(soln1D.grd.t(it))],'Interpreter','none','FontSize',12);
		phiAcurrent =		load([postproc.folderloc1D,'/phiA_',num2str(it),'.dat']);
		phiBcurrent =		load([postproc.folderloc1D,'/phiB_',num2str(it),'.dat']);
		phiAprev =			load([postproc.folderloc1D,'/phiA_',num2str(itprev),'.dat']);
		phiBprev =			load([postproc.folderloc1D,'/phiB_',num2str(itprev),'.dat']);
		phiAprevprev =		load([postproc.folderloc1D,'/phiA_',num2str(itprevprev),'.dat']);
		phiBprevprev =		load([postproc.folderloc1D,'/phiB_',num2str(itprevprev),'.dat']);
		phiAmBcurrent =		phiAcurrent-phiBcurrent;
		phiAmBprev =		phiAprev-phiBprev;
		phiAmBprevprev =	phiAprevprev-phiBprevprev;
		LHS =				zeros(1,soln1D.simul.nx);
		RHS =				zeros(1,soln1D.simul.nx);
		if (soln1D.simul.tstepprevs == 1)
			for ix = 1:soln1D.simul.nx
				LHS(ix) =	(phiAmBcurrent(ix)-phiAmBprev(ix))/(soln1D.grd.t(it)-soln1D.grd.t(itprev))+	...
							soln1D.soln.flow.vel.v(ix)*		...
							differentiation(phiAmBcurrent,soln1D.grd.x,ix,1,dirx{ix},postproc.derivwarn);
			end
		else
			for ix = 1:soln1D.simul.nx
				LHS(ix) =	differentiation([phiAmBprevprev(ix),phiAmBprev(ix),phiAmBcurrent(ix)],	...
							[soln1D.grd.t(itprevprev),soln1D.grd.t(itprev),soln1D.grd.t(it)],3,1,'bd',postproc.derivwarn) +	...
							soln1D.soln.flow.vel.v(ix)*		...
							differentiation(phiAmBcurrent,soln1D.grd.x,ix,1,dirx{ix},postproc.derivwarn);
			end
		end
		for ix = 2:soln1D.simul.nx-1
			RHS(ix) =		(soln1D.fltr.eta*soln1D.fltr.nondim.Pe*		...
							(soln1D.soln.flow.vel.v(ix)+soln1D.fltr.transpbc.frcdispers)+	...
							(1.0-double(soln1D.simul.nodiffusion)))*	...
							differentiation(phiAmBcurrent,soln1D.grd.x,ix,2,dirx{ix},postproc.derivwarn);
		end
		if (soln1D.geomdom.isradial ~= 0)
			for ix = 2:soln1D.simul.nx-1
				RHS(ix) =	RHS(ix)+	...
							(soln1D.fltr.eta*	...
							differentiation(soln1D.soln.flow.vel.v,soln1D.grd.x,ix,1,dirx{ix},postproc.derivwarn)+	...
							((soln1D.fltr.eta*soln1D.fltr.nondim.Pe*		...
							(soln1D.soln.flow.vel.v(ix)+soln1D.fltr.transpbc.frcdispers)+...
							(1.0-double(soln1D.simul.nodiffusion)))/(soln1D.grd.x(ix)/double(soln1D.geomdom.isradial))))*	...
							differentiation(phiAmBcurrent,soln1D.grd.x,ix,1,dirx{ix},postproc.derivwarn);
			end
		end
		plot(axl,soln1D.grd.x(2:soln1D.simul.nx-1),LHS(2:soln1D.simul.nx-1),':','color',[0.2,0.2,1.0],'linewidth',1.5);
		plot(axl,soln1D.grd.x(2:soln1D.simul.nx-1),(1.0/soln1D.fltr.nondim.Pe)*RHS(2:soln1D.simul.nx-1),	...
			'--','color',[0.2,1.0,0.2],'linewidth',0.2);
		if (soln1D.geomdom.isradial ~= 0)
			plot(axllin,soln1D.grd.x(2:soln1D.simul.nx-1),LHS(2:soln1D.simul.nx-1),':','color',[0.2,0.2,1.0],'linewidth',1.5);
			plot(axllin,soln1D.grd.x(2:soln1D.simul.nx-1),(1.0/soln1D.fltr.nondim.Pe)*RHS(2:soln1D.simul.nx-1),	...
				'--','color',[0.2,1.0,0.2],'linewidth',0.2);
		end
		legend({LHSstring,RHSstring},'interpreter','latex','fontsize',12,'location','eastoutside');
		pause(0.01);
		movframecurr =     getframe(figl);
		writeVideo(phimovvid,movframecurr);
		for ii = 1:2
			lastline = get(axl,'children');
			delete(lastline(1));
			if (soln1D.geomdom.isradial ~= 0)
				lastline = get(axllin,'children');
				delete(lastline(1));
			end
		end
		itprevprev =	itprev;
		itprev =		it;
	catch
	end
end
close(figl);
close(phimovvid);
