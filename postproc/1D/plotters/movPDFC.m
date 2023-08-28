% generating PDF of product with time
phiCmax =					-Inf;
phiCmin =					Inf;
soln_analyze.Cpdfbucket =	zeros(soln1D.simul.nt,nCpdfbucket+1);
for it = soln1D.simul.tstepsave
	phiCcurrent =			load([postproc.folderloc1D,'/phiC_',num2str(it),'.dat']);
	phiCcurrent =			phiCcurrent';
	phiCmax =				max([phiCmax,max(max(phiCcurrent))]);
	phiCmin =				min([phiCmin,min(min(phiCcurrent))]);
end
for it = soln1D.simul.tstepsave
	phiCcurrent =			load([postproc.folderloc1D,'/phiC_',num2str(it),'.dat']);
	phiCcurrent =			phiCcurrent';
	phiOS =					phiCcurrent-phiCmin;
	phiOSmax =				max(max(phiOS));
	phiOSnorm =				phiOS/(phiOSmax+eps);
	dc =					1.0/nCpdfbucket;
	blowupfac =				1.0/dc;
	phiOSnormre =			phiOSnorm*blowupfac;
	soln_analyze.Cpdfbucketlbl =						zeros(1,nCpdfbucket+1);
	for iCpdfbucket = 1:nCpdfbucket+1
		soln_analyze.Cpdfbucketlbl(iCpdfbucket) =		phiCmin+(phiCmax-phiCmin)*((iCpdfbucket-1)/nCpdfbucket);
	end
	for ix = 1:soln1D.simul.nx
		soln_analyze.Cpdfbucket(it,floor(phiOSnormre(ix))+1) =	...
														soln_analyze.Cpdfbucket(it,floor(phiOSnormre(ix))+1)+1.0;
	end
end
phimovvid =				VideoWriter([postproc.folderloc1D,'/PDFC.avi']);
open(phimovvid);
figl  =	figure('position',[100,100,750,500],'visible',postproc.vistogmov);
axl = axes(figl,'position',[0.20,0.20,0.75,0.75]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
set(axl,'xlim',[phiCmin,phiCmax]); set(axl,'ylim',[1.0E-5,1.25]);
set(axl,'yscale','log'); set(axl,'xscale','log');
xlabel('$\phi_C$','interpreter','latex');
ylabel('PDF($\phi_C$)','interpreter','latex');
for it = soln1D.simul.tstepsave
	plot(axl,soln_analyze.Cpdfbucketlbl,soln_analyze.Cpdfbucket(it,:)/sum(abs(soln_analyze.Cpdfbucket(it,:))),	...
		'k-','linewidth',1);
	title([postproc.casename,': time step is ',soln1D.grd.t(num2str(it))],'Interpreter','none','FontSize',12);
	pause(0.5);
	movframecurr =     getframe(figl);
	writeVideo(phimovvid,movframecurr);
	for ii = 1:1
		lastline = get(axl,'children');
		delete(lastline(1));
	end
end
close(figl);
close(phimovvid);
