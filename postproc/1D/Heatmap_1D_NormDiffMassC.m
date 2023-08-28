% initial clean-up and set-up
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');

% whether to additionally post-process front solution
iFront = 1;							%	plot front properties
% post-processing option toggles
vistog = 'on';						%	plot windows visibility toggle
gridtog = 'on';						%	grid toggle
derivwarn = 'true';					%	derivative computation warnings display toggle
forcenumfront = 1;					%	should the front properties be numerically computed even if post-proc solution is saved
forceMassCnum =	1;					%	should the mass of product be numerically computed even if post-proc solution is saved
% similarity solution (approximate analytical solution) computation parameters

% overriding front computation as per requirement
speciesname = ['A','B','C'];

% simulation output data structure
namegeom =						{'CI'};
eta =							[0.0E+0,	...
								1.00E-3,	3.33E-3,	1.00E-2,	3.33E-2,	1.00E-1,	3.33E-1,	...
								1.00E+0,	3.33E+0,	1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	...
								1.00E+3		];
etaname =						{'0d0EP0',	...
								'1d00EM3',	'3d33EM3',	'1d00EM2',	'3d33EM2',	'1d00EM1',	'3d33EM1',	...
								'1d00EP0',	'3d33EP0',	'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	...
								'1d00EP3'	};
Pe =							[1.00E+0,	3.33E+0,	1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	1.00E+3,	];
Pename =						{'1d00EP0',	'3d33EP0',	'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	'1d00EP3',	};
folderloc1Dparent =				'D:\Pratyaksh\LocalOnlyRawData_Level_0\1D\1D_Production_B';

namenormalize =					{'','_Normalized','_Ratio'};

for iGeom = 1:1
casename =						[namegeom{iGeom},'_Pe_',Pename{1},'_eta_',etaname{1}];
% location of specific case data
folderloc1D =					[folderloc1Dparent,'/',casename];

for iPe = 1:length(Pe)	

for ieta = 1:length(eta)
%	try to clear previous solution left-overs
	try
		clear soln1D;
		close soln_analyze;
		clear front;
	catch
	end
%	name of specific case
	casename =						[namegeom{iGeom},'_Pe_',Pename{iPe},'_eta_',etaname{ieta}];
%	location of specific case data
	folderloc1D =					[folderloc1Dparent,'/',casename];
%	analyzing data
%	loading solution and pre-setting differentiation directions
	if (isfile([folderloc1D,'/soln_posproc.mat']) ~= 0)
		load([folderloc1D,'/soln_posproc.mat']);
	else
		soln1D =								load([folderloc1D,'/soln.mat']);
		soln1D.simul.compMassC =				0;
		if (soln1D.simul.compMassC == 1)
			soln1D.solnreactive =				load([folderloc1D,'/solnreactive.mat']);	
			soln1D.solnreactive =				soln1D.solnreactive.solnreactive;
		end
	end
%	trimming step save array to unique values
	soln1D.simul.tstepsave =					unique(soln1D.simul.tstepsave);
	front.tsaved=								soln1D.grd.t(soln1D.simul.tstepsave);
%	obtaining mass of product C
	if ((exist('soln_analyze','var') == 0) || (forceMassCnum == 1))
		needmassCsoln =							1;
	else
		needmassCsoln =							0;
	end
	if (needmassCsoln == 1)
		[~,inod_larget] = min((front.tsaved-10).^2);
		it = soln1D.simul.tstepsave(inod_larget);
			phiCcurrent =						load([folderloc1D,'/phiC_',num2str(it),'.dat']);
			soln_analyze.MassC_HM(iPe,ieta) =	((double(soln1D.geomdom.isradial)*2*pi)^double(soln1D.geomdom.isradial~=0))*		...
												trapz(soln1D.grd.x,phiCcurrent.*(soln1D.grd.x.^soln1D.geomdom.isradial));
	end
end

end

for iPe = 1:length(Pe)
	for ieta = 2:length(eta)
		soln_analyze.MassC_HM(iPe,ieta) =	soln_analyze.MassC_HM(iPe,ieta)/soln_analyze.MassC_HM(iPe,1)-1;
	end
end

figl  =	figure('position',[100,35,1250,700],'visible','on'); axl = axes(figl,'position',[0.125,0.15,0.85,0.800]);
set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
set(gca,'linewidth',3); hold on;
ylabel('$Pe$','interpreter','latex');
xlabel('$\eta$','interpreter','latex');
set(gca,'xlim',[min(log10(-eps+min(eta(2:end)))),max(log10(eps+max(eta(2:end))))]); 
set(gca,'ylim',[min(log10(-eps+min(Pe))),max(log10(eps+max(Pe)))]); 
surface(log10(eta(2:end)),log10(Pe),log10(eps+soln_analyze.MassC_HM(1:end,2:end)),'LineStyle','none');
% plot(axl,log10(eps+postproc.eta_in),log10(eps+t_mid_in(iPe_nonsave,:)),':','color',[0.3,0.3,0.3],'linewidth',3);
xticks([-3, -1, 1, 3]);
yticks([0, 1, 2, 3]);
xticklabels({'$10^-3}$', '$10^{-1}$', '$10^1$', '$10^3$'});
yticklabels({'$10^{0}$', '$10^{1}$', '$10^2$', '$10^{3}$'});
saveas(figl,[folderloc1Dparent,'/',namegeom{iGeom},'_NormDiffMassCHM'],'fig');
close(figl);

end
