% plotting terms of GDE

% plotting GDE terms at injection line/point
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'xscale','log'); set(axl,'yscale','log'); set(gca,'ylim',[1E-10,1E10]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
xlabel('$t$','interpreter','latex');
ylabel('GDE terms','interpreter','latex');
for iend = 1:2
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				ittrim = front.trnszon(isolvmode,ifracTransitn,ilimdef,iend,2);
				for ispecies = 1:4
					for iterm = 1:2
						colline = [0.0+1.0*(ispecies==1),0.0+1.0*(ispecies==2),0.0+1.0*(ispecies==3)];
						thickterm = iterm;
						plot(axl,front.tsaved(1:ittrim),abs(dphidr_centre(1:ittrim,ispecies,iterm)'),		...
						'-','color',colline,'linewidth',thickterm);
					end
				end
			end
		end
	end
end
legend({'$\displaystyle \left|\frac{\partial \phi_{A}}{\partial r}\right|$',		...
		'$\displaystyle \left|\frac{\partial \phi_{B}}{\partial r}\right|$',		...
		'$\displaystyle \left|\frac{\partial \phi_{C}}{\partial r}\right|$',		...
		'$\displaystyle \left|\frac{\partial \theta}{\partial r}\right|$',         ...
		'$\displaystyle \left|\frac{\partial^2 \phi_{A}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\frac{\partial^2 \phi_{B}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\frac{\partial^2 \phi_{C}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\frac{\partial^2 \theta}{\partial r^2}\right|$'},	...			
		'interpreter','latex','fontsize',18,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Terms_GDE_Centre'],'fig');
close(figl);
% plotting GDE terms at front
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'xscale','log'); set(axl,'yscale','log'); set(gca,'ylim',[1E-10,1E10]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
xlabel('$t$','interpreter','latex');
ylabel('GDE terms','interpreter','latex');
for iend = 1:2
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				ittrim = front.trnszon(isolvmode,ifracTransitn,ilimdef,iend,2);
				for ispecies = 1:4
					for iterm = 1:2
						colline = [0.0+1.0*(ispecies==1),0.0+1.0*(ispecies==2),0.0+1.0*(ispecies==3)];
						thickterm = iterm;
						plot(axl,front.tsaved(1:ittrim),abs(dphidr_f(1:ittrim,ispecies,iterm)'),		...
						'-','color',colline,'linewidth',thickterm);
					end
				end
				for iterm = 1:2
					ispecies = 4;
					colline = [0.1+0.5*(ispecies==1),0.1+0.5*(ispecies==2),0.1+0.5*(ispecies==3)];
					thickterm = iterm;
					plot(axl,front.tsaved(1:ittrim),abs(1./(front.num.w_f_orig(1:ittrim,2,1).^iterm)'),		...
					'--','color',colline,'linewidth',thickterm);
				end
			end
		end
	end
end
legend({'$\displaystyle \left|\frac{\partial \phi_{A}}{\partial r}\right|$',		...
		'$\displaystyle \left|\frac{\partial \phi_{B}}{\partial r}\right|$',		...
		'$\displaystyle \left|\frac{\partial \phi_{C}}{\partial r}\right|$',		...
		'$\displaystyle \left|\frac{\partial \theta}{\partial r}\right|$',         ...
		'$\displaystyle \left|\frac{\partial^2 \phi_{A}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\frac{\partial^2 \phi_{B}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\frac{\partial^2 \phi_{C}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\frac{\partial^2 \theta}{\partial r^2}\right|$'},	...			
		'interpreter','latex','fontsize',18,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Terms_GDE_Front'],'fig');
close(figl);
% plotting (full) GDE terms at injection line/point
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'xscale','log'); set(axl,'yscale','log'); set(gca,'ylim',[1E-10,1E10]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
xlabel('$t$','interpreter','latex');
ylabel('GDE terms','interpreter','latex');
for iend = 1:2
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				ittrim = front.trnszon(isolvmode,ifracTransitn,ilimdef,iend,2);
				for ispecies = 1:4
					for iterm = 1:2
						colline = [0.0+1.0*(ispecies==1),0.0+1.0*(ispecies==2),0.2+0.8*(ispecies==3)];
						thickterm = iterm;
						plot(axl,front.tsaved(1:ittrim),abs(dphidr_centre_full(1:ittrim,ispecies,iterm)'),		...
						'-','color',colline,'linewidth',thickterm);
					end
				end
			end
		end
	end
end
legend({'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \phi_{A}}{\partial r}\right|$',		...
		'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \phi_{B}}{\partial r}\right|$',		...
		'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \phi_{C}}{\partial r}\right|$',		...
		'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \theta}{\partial r}\right|$',         ...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \phi_{A}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \phi_{B}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \phi_{C}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \theta}{\partial r^2}\right|$'},	...			
		'interpreter','latex','fontsize',18,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Terms_Full_GDE_Centre'],'fig');
close(figl);
% plotting (full) GDE terms at front
figl  =	figure('position',[100,100,1200,500],'visible',postproc.vistog);
axl = axes(figl,'position',[0.10,0.20,0.75,0.75]);
set(axl,'xscale','log'); set(axl,'yscale','log'); set(gca,'ylim',[1E-10,1E10]);
set(axl,'box','on'); grid(postproc.gridtog); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); hold on;
xlabel('$t$','interpreter','latex');
ylabel('GDE terms','interpreter','latex');
for iend = 1:2
	for isolvmode = 1:2
		for ilimdef = 1:2
			for ifracTransitn = 1:length(postproc.fracTransitn)
				ittrim = front.trnszon(isolvmode,ifracTransitn,ilimdef,iend,2);
				for ispecies = 1:4
					for iterm = 1:2
						colline = [0.0+1.0*(ispecies==1),0.0+1.0*(ispecies==2),0.0+1.0*(ispecies==3)];
						thickterm = iterm;
						plot(axl,front.tsaved(1:ittrim),abs(dphidr_f_full(1:ittrim,ispecies,iterm)'),		...
						'-','color',colline,'linewidth',thickterm);
					end
				end
				for iterm = 1:2
					ispecies = 4;
					colline = [0.1+0.5*(ispecies==1),0.1+0.5*(ispecies==2),0.1+0.5*(ispecies==3)];
					thickterm = iterm;
					plot(axl,front.tsaved(1:ittrim),	...
					abs(multiplier_f(1:ittrim)'.*(1./(front.num.w_f_orig(1:ittrim,2,1).^iterm))'),...
					'--','color',colline,'linewidth',thickterm);
				end
			end
		end
	end
end
legend({'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \phi_{A}}{\partial r}\right|$',		...
		'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \phi_{B}}{\partial r}\right|$',		...
		'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \phi_{C}}{\partial r}\right|$',		...
		'$\displaystyle \left|\left(\frac{Pe-1}{r}\right)\frac{\partial \theta}{\partial r}\right|$',         ...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \phi_{A}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \phi_{B}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \phi_{C}}{\partial r^2}\right|$',	...
		'$\displaystyle \left|\left(\frac{\eta Pe}{r}+1\right)\frac{\partial^2 \theta}{\partial r^2}\right|$'},	...			
		'interpreter','latex','fontsize',18,'location','eastoutside');
saveas(figl,[postproc.folderloc1D,'/Terms_Full_GDE_Front'],'fig');
close(figl);