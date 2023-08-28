% initial clean-up and set-up [default definitions taken for all variables]
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:MKDIR:DirectoryExists');

% specifying parent folder location
folderloc1Dparent = 'D:\LocalOnlyRawData_Level_0\1D\1D_Production_B';

% specifying etaeter values and names
etalabel =						'$\log{\eta}$';
Pelabel =						'$\log{Pe}$';
contrastname =					{'ratio', 'diff', 'normdiff','baseval'};
propname =						{'rfSET', 'wfSET', 'barRET', 'MCET'};
eta =							[ ...
								0.0E+0,		...
								1.00E-3,	3.33E-3,	1.00E-2,	3.33E-2,	...
								1.00E-1];
etaname =						{	...
								'0d0EP0',		...
								'1d00EM3',	'3d33EM3',	'1d00EM2',	'3d33EM2',	...
								'1d00EM1'};
Pe =							flip([1.00E+1,		3.33E+1,	1.00E+2,	3.33E+2,	1.00E+3		]);
Pename =						flip({'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	'1d00EP3'	});

% initializing time grids
contrast_FrProp =				zeros(length(Pe),length(eta),4,4);
t_endtime =						zeros(length(Pe),length(eta));

% looping over etaeter values
for iPe = 1:length(Pe)
	for ieta = 1:length(eta)
%		setting output folder
		clear params;
		casename = ['CI_Pe_',Pename{iPe},'_eta_',etaname{ieta}]; 
		folderloc1D = [folderloc1Dparent,'/',casename];
		params = load([folderloc1D,'/soln.mat']);
		t_endtime(iPe,ieta) =	params.grd.t(params.simul.tstepsave(end));
%		obtaining last time-step solution
		phiAET =								load([folderloc1D,'/phiA_',num2str(params.simul.tstepsave(end)),'.dat']);
		phiBET =								load([folderloc1D,'/phiB_',num2str(params.simul.tstepsave(end)),'.dat']);
		phiCET =								load([folderloc1D,'/phiC_',num2str(params.simul.tstepsave(end)),'.dat']);
		thetaET =								phiAET-phiBET;
		RET =									phiAET.*phiBET;
		[~,frontnode] =							min(thetaET.^2);
		rfET =									params.grd.x(frontnode);
		barRET =								trapz(params.grd.x,params.grd.x.^2.*RET);
		wfET =									trapz(params.grd.x,params.grd.x.^2.*(params.grd.x-rfET).*RET)/barRET;
		MCET =									trapz(params.grd.x,params.grd.x.^2.*phiCET);
%		obtaining contrast data
		if (ieta == 1)
			rfET_base = rfET; barRET_base = barRET; wfET_base = wfET; MCET_base = MCET;
		end
		contrast_FrProp(iPe,ieta,1,1) =				rfET/rfET_base;
		contrast_FrProp(iPe,ieta,1,2) =				wfET/wfET_base;
		contrast_FrProp(iPe,ieta,1,3) =				barRET/barRET_base;
		contrast_FrProp(iPe,ieta,1,4) =				MCET/MCET_base;
		contrast_FrProp(iPe,ieta,2,1) =				rfET-rfET_base;
		contrast_FrProp(iPe,ieta,2,3) =				wfET-wfET_base;
		contrast_FrProp(iPe,ieta,2,3) =				barRET-barRET_base;
		contrast_FrProp(iPe,ieta,2,4) =				MCET-MCET_base;
		contrast_FrProp(iPe,ieta,3,1) =				rfET/rfET_base-1;
		contrast_FrProp(iPe,ieta,3,2) =				wfET/wfET_base-1;
		contrast_FrProp(iPe,ieta,3,3) =				barRET/barRET_base-1;
		contrast_FrProp(iPe,ieta,3,4) =				MCET/MCET_base-1;
		contrast_FrProp(iPe,ieta,4,1) =				rfET;
		contrast_FrProp(iPe,ieta,4,2) =				wfET;
		contrast_FrProp(iPe,ieta,4,3) =				barRET;
		contrast_FrProp(iPe,ieta,4,4) =				MCET;
		disp(['Done with Pe = ',num2str(Pe(iPe)),', eta = ',num2str(eta(ieta))]);
	end
end
% generating heatmaps
[logetaeta, logPePe] = meshgrid(log10(eta(2:end)),log10(Pe));
for icontrast = 1:4
	for iprop = 1:4
		pcolor(logetaeta,logPePe,real(log10(contrast_FrProp(:,2:end,icontrast,iprop)))); shading interp;
		saveas(gcf,[folderloc1Dparent,'/CI_HM_',contrastname{icontrast},'_',propname{iprop},'_Log_Real'],'fig');
		close(gcf);
		pcolor(logetaeta,logPePe,imag(log10(contrast_FrProp(:,2:end,icontrast,iprop)))); shading interp;
		saveas(gcf,[folderloc1Dparent,'/CI_HM_',contrastname{icontrast},'_',propname{iprop},'_Log_Imag'],'fig');
		close(gcf);
		pcolor(logetaeta,logPePe,contrast_FrProp(:,2:end,icontrast,iprop)); shading interp;
		saveas(gcf,[folderloc1Dparent,'/CI_HM_',contrastname{icontrast},'_',propname{iprop}],'fig');
		close(gcf);
	end
end
% saving final output
save([folderloc1Dparent,'/CI_HM_data.mat']);