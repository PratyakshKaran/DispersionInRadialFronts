% initial clean-up and set-up
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');
warning('off','MATLAB:nearlySingularMatrix');

% postprocessing parameters (see PostProc_1D for definition)
postproc.obtainfrontnumforced =	1;
postproc.facnearstat =			0.995;
postproc.factrans =				1.0;

% specifying parent folder location
%{
postproc.folderloc1Dparent =	'P:\LocalOnlyRawData_Level_0\1D\1D_Production_B';
%}
postproc.folderloc1Dparent =	'E:\LocalOnlyRawData_Level_0\1D\1D_Production_C [onlySI]';
%
% specifying values and names
%{
postproc.geom =					'CI';
postproc.etalabel =				'$\eta$';
postproc.eta =					[	...
								1.00E-3,	3.33E-3,	1.00E-2,	3.33E-2,	...
								1.00E-1,	3.33E-1,	1.00E+0,	3.33E+0,	...
								1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	...
								1.00E+3];
postproc.etaname =				{	...
								'1d00EM3',	'3d33EM3',	'1d00EM2',	'3d33EM2',	...
								'1d00EM1',	'3d33EM1',	'1d00EP0',	'3d33EP0',	...
								'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	...
								'1d00EP3'};
postproc.Pe =					[1.00E+2,	];
postproc.Pename =				{'1d00EP2',	};
%}
%
postproc.geom =					'SI';
postproc.etalabel =				'$\eta$';
postproc.eta =					[	...
								1.00E-3,	3.33E-3,	1.00E-2,	3.33E-2,	...
								1.00E-1,	3.33E-1,	1.00E+0,	3.33E+0,	...
								1.00E+1,	3.33E+1,	1.00E+2,	3.33E+2,	...
								1.00E+3,	3.33E+3,	1.00E+4,	3.33E+4,	...
								1.00E+5,	3.33E+5,							...
								1.00E+6];
postproc.etaname =				{	...
								'1d00EM3',	'3d33EM3',	'1d00EM2',	'3d33EM2',	...
								'1d00EM1',	'3d33EM1',	'1d00EP0',	'3d33EP0',	...
								'1d00EP1',	'3d33EP1',	'1d00EP2',	'3d33EP2',	...
								'1d00EP3',	'3d33EP3',	'1d00EP4',	'3d33EP4',	...
								'1d00EP5',	'3d33EP5',							...
								'1d00EP3'};
postproc.Pe =					[1.00E+0,	1.00E+1,	1.00E+2,	1.00E+3		];
postproc.Pename =				{'1d00EP0',	'1d00EP1',	'1d00EP2',	'1d00EP3'	};
%

% output smoothening parameters
postproc.dosmooth =				1;
postproc.n_eta_in =				1001;
postproc.eta_in =				exp(linspace(log(postproc.eta(1)),log(postproc.eta(end)),postproc.n_eta_in));

% initializing time grids
nod_transA =					ones(length(postproc.Pename),length(postproc.eta));
nod_transB =					ones(length(postproc.Pename),length(postproc.eta));
t_transA =						zeros(length(postproc.Pename),length(postproc.eta));
t_transB =						zeros(length(postproc.Pename),length(postproc.eta));
t_transB_orig =					zeros(length(postproc.Pename),length(postproc.eta));
rfS =							zeros(length(postproc.Pename),length(postproc.eta));
t_transA_in =					zeros(length(postproc.Pename),postproc.n_eta_in);
t_transB_in =					zeros(length(postproc.Pename),postproc.n_eta_in);

% looping over parameter values
for iPe_nonsave = 1:length(postproc.Pename)
%	try
%		postproctemp2 =		postproc;
		% iPetemp2 =			iPe_nonsave;
	%	load([postproc.folderloc1Dparent,'/',postproc.geom,'_RegimePlotData_Pe_',postproc.Pename{iPe_nonsave},'.mat']);
	%	postproc =			postproctemp2;
	% 	iPe_nonsave =				iPetemp2;
%	catch
		for ieta = 1:length(postproc.etaname)
			try
				clear soln1D;
				close soln_analyze;
				clear front;
			catch
			end
%			specifying casename
			casename =				[postproc.geom,'_Pe_',postproc.Pename{iPe_nonsave},'_eta_',postproc.etaname{ieta}];
			folderloc1D =			[postproc.folderloc1Dparent,'/',casename];
			disp(casename);
%			loading solution and pre-setting differentiation directions
			postproctemp =			postproc;
			soln1D =				load([folderloc1D,'/soln.mat']);
			postproc =				postproctemp;
			soln1D.simul.tstepsave=	unique(soln1D.simul.tstepsave);
%			obtaining time grid that is saved to files
			front.tsaved =			soln1D.grd.t(soln1D.simul.tstepsave);
			if ((iPe_nonsave == 1) && (ieta == 1))
				front_num_r_f_orig=	zeros(length(postproc.Pename),length(postproc.eta),length(soln1D.simul.tstepsave));
			end
%			obtaining front location
			ifront.tsaved =			0;				
			for it = soln1D.simul.tstepsave
				ifront.tsaved =		ifront.tsaved+1;
				phiCcurrent =		load([folderloc1D,'/phiC_',num2str(it),'.dat']);
				phiAcurrent =		load([folderloc1D,'/phiA_',num2str(it),'.dat']);
				phiBcurrent =		load([folderloc1D,'/phiB_',num2str(it),'.dat']);
				front.num.R =		phiAcurrent.*phiBcurrent;
				front.num.theta =	phiAcurrent-phiBcurrent;
				[~,inodrf_orig] =	min(front.num.theta.^2);
				front.num.r_f_orig_node(ifront.tsaved) =		inodrf_orig;
				front.num.r_f_orig(ifront.tsaved) =				soln1D.grd.x(inodrf_orig);
				front_num_r_f_orig(iPe_nonsave,ieta,ifront.tsaved) =	soln1D.grd.x(inodrf_orig);
			end
%			obtaining stationary solution
			soln1D.simul.STLtol =						1E-5;
			soln1D.simul.STLrelax =						0.99;
			if (soln1D.geomdom.isradial == 2)
				phiS =										...
				stationary1D_reactive(soln1D.geomdom,soln1D.simul,soln1D.fltr,soln1D.grd,soln1D.soln.flow);
				if (soln1D.fltr.eta == 0)
					rfS(iPe_nonsave,ieta) =						soln1D.fltr.nondim.Pe/log(1+soln1D.fltr.nondim.gamm);
				else
					rfS(iPe_nonsave,ieta) =						...
					sqrt(soln1D.fltr.eta*soln1D.fltr.nondim.Pe)*			...
					tan(sqrt(soln1D.fltr.eta/(soln1D.fltr.nondim.Pe))*	...
					log(1-(1/(1+soln1D.fltr.nondim.gamm))*									...
					(1-exp((pi/2)*sqrt(soln1D.fltr.nondim.Pe/soln1D.fltr.eta)))));
				end
			end
		end
		avoidVariable = 'postproc, iPe, ieta, postproctemp, postproctemp1, postproctemp2 iPetemp2, iPe_nonsave';
		save([postproc.folderloc1Dparent,'/',postproc.geom,'_RegimePlotData_Pe_',postproc.Pename{iPe_nonsave},'.mat'], ...
		'-regexp', ['^(?!', avoidVariable,'$).']);
	end

%	obtaining transition zones
	for ieta = 1:length(postproc.etaname)
		if (soln1D.geomdom.isradial == 2)
			[~,nod_transA(iPe_nonsave,ieta)] =			...
			min(((front_num_r_f_orig(iPe_nonsave,ieta,:).^2)-postproc.factrans*postproc.eta(ieta)*postproc.Pe(iPe_nonsave)).^2);
			nod_transA(iPe_nonsave,ieta) =				soln1D.simul.tstepsave(nod_transA(iPe_nonsave,ieta));
			t_transA(iPe_nonsave,ieta) =				soln1D.grd.t(nod_transA(iPe_nonsave,ieta));
			if (front_num_r_f_orig(iPe_nonsave,ieta,end)^2 < (1/postproc.factrans)*postproc.eta(ieta)*postproc.Pe(iPe_nonsave))
				t_transA(iPe_nonsave,ieta) =			1E100;
			end
%			postproc.facnearstat
			[~,nod_transB(iPe_nonsave,ieta)] =			min((front_num_r_f_orig(iPe_nonsave,ieta,:)- ...
														postproc.facnearstat*rfS(iPe_nonsave,ieta)).^2);
%			nod_transB(iPe_nonsave,ieta)
%			nod_transB(iPe_nonsave,ieta) =				soln1D.simul.tstepsave(nod_transB(iPe_nonsave,ieta));
			t_transB(iPe_nonsave,ieta) =				soln1D.grd.t(nod_transB(iPe_nonsave,ieta));
		else
			[~,nod_transA(iPe_nonsave,ieta)] =			...
			min(((front_num_r_f_orig(iPe_nonsave,ieta,:))-postproc.eta(ieta)*postproc.Pe(iPe_nonsave)).^2);
			nod_transA(iPe_nonsave,ieta) =				soln1D.simul.tstepsave(nod_transA(iPe_nonsave,ieta));
			t_transA(iPe_nonsave,ieta) =				soln1D.grd.t(nod_transA(iPe_nonsave,ieta));
			if (nod_transA(iPe_nonsave,ieta) == length(soln1D.grd.t))
				t_transA(iPe_nonsave,ieta) =			pchip(postproc.eta(max([1,ieta-3]):ieta-1), ...
														t_transA(iPe_nonsave,max([1,ieta-3]):ieta-1), ...
														postproc.eta(ieta));
			end
			alph =										0.475417;
			bet =										alph^6;
			t_transB(iPe_nonsave,ieta) =				bet*(postproc.eta(ieta)^2)*(postproc.Pe(iPe_nonsave)^3)*	...
														((gammaincinv(1/2,postproc.Pe(iPe_nonsave)/2,'upper')).^(-3));
		end
	end
%	flipping the transition times if large Pe
%{
	if (postproc.Pe(iPe_nonsave) >= 1.00E+01)
		t_transB_orig(iPe_nonsave,:) =		t_transB(iPe_nonsave,:);
		t_transB(iPe_nonsave,:) =			t_transA(iPe_nonsave,:);
		t_transA(iPe_nonsave,:) =			t_transB_orig(iPe_nonsave,:);
	end
%}
%	obtaining profile to generate regime plot
	if (postproc.dosmooth == 1)
		t_transA_in(iPe_nonsave,:) =	exp(pchip(log(eps+postproc.eta),log(eps+t_transA(iPe_nonsave,:)),log(postproc.eta_in)));
		t_transB_in(iPe_nonsave,:) =	exp(pchip(log(eps+postproc.eta),log(eps+t_transB(iPe_nonsave,:)),log(postproc.eta_in)));
	else
		postproc.n_eta_in =		length(postproc.eta);
		postproc.eta_in =		postproc.eta;
		t_transA_in =			t_transA;
		t_transB_in =			t_transB;
		t_bot_in =				t_bot;
		t_mid_in =				t_mid;
		t_top_in =				t_top;
	end

%	plotting regime
	[~,inodeta_smaller] = min((postproc.eta_in-0.1*postproc.Pe(iPe_nonsave)).^2);
	[~,inodeta] = min((postproc.eta_in-postproc.Pe(iPe_nonsave)).^2);
	[~,inodeta_larger] = min((postproc.eta_in-10*postproc.Pe(iPe_nonsave)).^2);
	[~,inodcrossA] =	min((t_transA_in(iPe_nonsave,:)-1).^2);
	[~,inodcrossB] =	min((t_transA_in(iPe_nonsave,:)-t_transB_in(iPe_nonsave,:)).^2);
	inodcrossA =		min([length(postproc.eta_in)-1,max([2,inodcrossA])]);
	inodcrossB =		min([length(postproc.eta_in)-1,max([2,inodcrossB])]);
	if (inodcrossA == inodcrossB)
		inodcrossB =	inodcrossB-1;
	end
	figl  =	figure('position',[100,35,1250,700],'visible','on'); axl = axes(figl,'position',[0.125,0.15,0.85,0.800]);
	set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
	set(gca,'linewidth',3); hold on;
	ylabel('$t$','interpreter','latex');
	xlabel(postproc.etalabel,'interpreter','latex');
	set(gca,'xlim',[min(log10(-eps+postproc.eta_in)),max(log10(eps+postproc.eta_in))]); 
	if (soln1D.geomdom.isradial == 2)
		set(gca,'ylim',[-7,11]);
	else
		set(gca,'ylim',[-7,5]);
	end
	etanodrange =		sort([1,inodeta_smaller,inodeta,inodeta_larger,length(postproc.eta_in)]);
	t_layer =			zeros(8,length(postproc.eta_in));
	for ieta = 1:length(postproc.eta_in)
		t_trans_temp =		sort([log10(eps+t_transA_in(iPe_nonsave,ieta)),log10(eps+t_transB_in(iPe_nonsave,ieta))]);
		t_layer(1,ieta) =	-7.0;
		t_layer(4,ieta) =	0.0;
		t_layer(7,ieta) =	11.0;
		if ((t_trans_temp(1) <= 0.0) && (t_trans_temp(2) <= 0.0))
			t_layer(2,ieta) =	t_trans_temp(1);
			t_layer(3,ieta) =	t_trans_temp(2);
			t_layer(5,ieta) =	0.0;
			t_layer(6,ieta) =	0.0;
		elseif ((t_trans_temp(1) <= 0.0) && (t_trans_temp(2) >= 0.0))
			t_layer(2,ieta) =	t_trans_temp(1);
			t_layer(3,ieta) =	0.0;
			t_layer(5,ieta) =	0.0;
			t_layer(6,ieta) =	t_trans_temp(2);
		else
			t_layer(2,ieta) =	0.0;
			t_layer(3,ieta) =	0.0;
			t_layer(5,ieta) =	t_trans_temp(1);
			t_layer(6,ieta) =	t_trans_temp(2);
		end
	end
	for ietanod = 1:length(etanodrange)-1
%		etanodrange(4) = 1000;
		area(axl,log10(eps+postproc.eta_in(etanodrange(ietanod):etanodrange(ietanod+1))),	...
		[t_layer(1,etanodrange(ietanod):etanodrange(ietanod+1));		...
		t_layer(2,etanodrange(ietanod):etanodrange(ietanod+1))-t_layer(1,etanodrange(ietanod):etanodrange(ietanod+1));		...
		t_layer(3,etanodrange(ietanod):etanodrange(ietanod+1))-t_layer(2,etanodrange(ietanod):etanodrange(ietanod+1));		...
		t_layer(4,etanodrange(ietanod):etanodrange(ietanod+1))-t_layer(3,etanodrange(ietanod):etanodrange(ietanod+1));		...
		t_layer(5,etanodrange(ietanod):etanodrange(ietanod+1))-t_layer(4,etanodrange(ietanod):etanodrange(ietanod+1));		...
		t_layer(6,etanodrange(ietanod):etanodrange(ietanod+1))-t_layer(5,etanodrange(ietanod):etanodrange(ietanod+1));		...
		t_layer(7,etanodrange(ietanod):etanodrange(ietanod+1))-t_layer(6,etanodrange(ietanod):etanodrange(ietanod+1));		...
		]',	...
		'linestyle','none');
	end
%	plot(axl,log10(eps+postproc.eta_in),log10(eps+t_mid_in(iPe_nonsave,:)),':','color',[0.3,0.3,0.3],'linewidth',3);
	if ((soln1D.geomdom.isradial == 1) && (postproc.Pe(iPe_nonsave) >= 1.00E+01))
		plot(axl,log10(eps+postproc.eta_in),log10(eps+t_transB_in(iPe_nonsave,:)),'-','color',[0.9,0.1,0.1],'linewidth',3);
		plot(axl,log10(eps+postproc.eta_in),log10(eps+t_transA_in(iPe_nonsave,:)),'-','color',[0.1,0.9,0.3],'linewidth',3);
	else
		plot(axl,log10(eps+postproc.eta_in),log10(eps+t_transA_in(iPe_nonsave,:)),'-','color',[0.9,0.1,0.1],'linewidth',3);
		plot(axl,log10(eps+postproc.eta_in),log10(eps+t_transB_in(iPe_nonsave,:)),'-','color',[0.1,0.9,0.3],'linewidth',3);
	end
	if (soln1D.geomdom.isradial == 2)
		plot(axl,log10(eps+postproc.eta_in(inodeta_smaller))*[1,1],[-7,11],':','color',[0.3,0.3,0.3],'linewidth',3);
		plot(axl,log10(eps+postproc.eta_in(inodeta))*[1,1],[-7,11],':','color',[0.3,0.3,0.3],'linewidth',3);
		plot(axl,log10(eps+postproc.eta_in(inodeta_larger))*[1,1],[-7,11],':','color',[0.3,0.3,0.3],'linewidth',3);
		plot(axl,log10(eps+postproc.eta_in(1:inodeta_larger)), ...
		log10(eps+(1/3)*postproc.eta_in(1:inodeta_larger).^(3/2)*postproc.Pe(iPe_nonsave)^(3/2)),	...
		'-','color',[1,0.2,0.2],'linewidth',1.5);
	else
		plot(axl,log10(postproc.eta_in),log10(eps+(0.500000000000000*postproc.eta_in.^2*postproc.Pe(iPe_nonsave)^2)), ...
		'-.','color',[0.9,0.1,0.1],'linewidth',1.5);
		plot(axl,log10(postproc.eta_in),log10(eps+(1.163284365916267*postproc.eta_in.^2*postproc.Pe(iPe_nonsave)^3)), ...
		'--','color',[0.9,0.1,0.1],'linewidth',1.5);
		plot(axl,log10(postproc.eta_in), ...
		log10(eps+(0.0115465*(postproc.eta_in.^2)*(postproc.Pe(iPe_nonsave)^3)* ...
		((gammaincinv(1/2,postproc.Pe(iPe_nonsave)/2,'upper'))^(-3)))), ...
		'-.','color',[0.1,0.9,0.3],'linewidth',1.5);
		plot(axl,log10(postproc.eta_in),log10(eps+(8*0.0115465*postproc.eta_in.^2)), ...
		'--','color',[0.1,0.9,0.3],'linewidth',1.5);
	end
	if (soln1D.geomdom.isradial == 2)
		xticks([-3, 0, 3, 6]);
		yticks([-7, -1, 5, 11]);
		xticklabels({'$10^{-3}$', '$10^{0}$', '$10^3$', '$10^6$'});
		yticklabels({'$10^{-7}$', '$10^{-1}$', '$10^5$', '$10^{11}$'});
	else
		xticks([-3, -1, 1, 3]);
		yticks([-7, -4, -1, 2, 5]);
		xticklabels({'$10^{-3}$', '$10^{-1}$', '$10^{1}$', '$10^3$'});
		yticklabels({'$10^{-7}$', '$10^{-4}$', '$10^{-1}$', '$10^2$', '$10^5$'});
	end
	view(90,-90);
	saveas(figl,[postproc.folderloc1Dparent,'/',postproc.geom,'_RegimePlot_Pe_',postproc.Pename{iPe_nonsave}],'fig');
	close(figl);
	
% end
