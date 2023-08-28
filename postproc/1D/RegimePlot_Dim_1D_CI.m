% initial clean-up and set-up
clear;
close all;
clc;
addpath('../simroutines/functions');
addpath('./anasolvers');
addpath('./plotters');
warning('off','MATLAB:nearlySingularMatrix');

% specifying parent folder location
%
postproc.folderloc1Dparent =	'B:\RegimePlotDim_CI';

% setting parameter values
Name =		'Pe_27778_lL_10cm';
Q0 =		2.7778E-4;
DM =		1E-8;
lL =		1E-1;
Pe =		Q0/DM;
logtst =	0;
logtend =	11;
logtRst =	0;
logtRend =	11;
nt =		101;
ntR =		101;
nhatch =	51;

% obtaining variables
t_R =		10.^linspace(logtst,logtend,nt);
t =			10.^linspace(logtRst,logtRend,ntR);
R =			sqrt(Q0*t);
eta =		lL./sqrt(Q0*t_R);
ttrL =		1.16*eta.^2.*t_R;
ttrP =		((eta.^2*Pe^3)/(4*gammaincinv(1/2,Pe/2,'upper'))).*t_R;
r_f_Disp =	sqrt((4*gammaincinv(1/2,Pe/2,'upper'))/Pe)*sqrt(t_R./t_R).*sqrt(Q0*t_R);
r_f_Diff =	(9*gammaincinv(1/2,1/3,'upper')*eta).^(1/3).*sqrt(Q0*t_R);
[~,ind1sttrans] =		min((r_f_Disp-eta).^2);
[~,ind2ndtrans] =		min((r_f_Diff-eta*Pe).^2);

% generating sqrt Q0tbylL regime diagram
figl  =	figure('position',[50,21.8,800,731.2],'visible','on'); 
axl = axes(figl,'position',[0.14675,0.119116238449763,0.7,0.818285714285715]);
set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
set(gca,'linewidth',6); set(gca,'xscale','log'); set(gca,'yscale','log'); hold on;
ylabel('$t_R/l_L^2 = \left(k_Rc_0\right)^{-1}/l_L^2$','interpreter','latex');
xlabel("$\bar{r}'/l_L = \sqrt{Q_0 t'}/l_L$",'interpreter','latex');
xx = [sqrt(Q0*t(1))/lL, sqrt(Q0*ttrP(1))/lL, sqrt(Q0*ttrP(1))/lL];
yy = [t_R(1)/(lL^2), ttrP(1)/(lL^2), t_R(1)/(lL^2)]; 
patch(xx, yy, [1.00,0.70,0.70],'linestyle','none');	
xx = [sqrt(Q0*t(1))/lL, sqrt(Q0*ttrP(1))/lL, sqrt(Q0*ttrP(1))/lL, sqrt(Q0*t(1))/lL]; 
yy = [t_R(1)/(lL^2), ttrP(1)/(lL^2), t_R(end)/(lL^2), t_R(end)/(lL^2)]; 
patch(xx, yy, [0.80,0.50,0.50],'linestyle','none');
xx = [sqrt(Q0*ttrP(1))/lL, sqrt(Q0*ttrP(1))/lL, sqrt(Q0*t(end))/lL, sqrt(Q0*t(end))/lL]; 
yy = [t_R(1)/(lL^2), ttrP(1)/(lL^2), t_R(end)/(lL^2), t_R(1)/(lL^2)]; 
patch(xx, yy, [0.70,0.70,1.00],'linestyle','none');
xx = [sqrt(Q0*ttrP(1))/lL, sqrt(Q0*t(end))/lL, sqrt(Q0*ttrP(1))/lL]; 
yy = [ttrP(1)/(lL^2), t_R(end)/(lL^2), t_R(end)/(lL^2)]; 
patch(xx, yy, [0.50,0.50,0.80],'linestyle','none');
plot(axl,sqrt(Q0*ttrL)/lL,t_R/(lL^2),':','color',[0.0,0.0,0.0],'linewidth',2);
plot(axl,sqrt(Q0*ttrP)/lL,t_R/(lL^2),'-.','color',[0.0,0.0,0.0],'linewidth',2);
hatchht = 10.^linspace(log10(t_R(1)/(lL^2)),log10(t_R(end)/(lL^2)),nhatch);
for ihatch = 1:nhatch
	plot([sqrt(Q0*mean(ttrL))/lL,sqrt(Q0*mean(ttrP))/lL],hatchht(ihatch)*[1,1],'k-','linewidth',1);
end
plot(axl,R/lL,t_R/(lL^2),'-','color',[0.0,0.0,0.0],'linewidth',2);
yticks([1E-9, 1E-6, 1E-3, 1, 3600*24, 3600*24*365*100, 3600*24*365*1E6]);
yticklabels({'ns/m^2', '$\mu$s/m^2','ms/m^2', 's/m^2', 'day/m^2', 'cent/m^2', 'Myr/m^2'});
xticks([1E-6, 1E-4, 1E-2, 1E0, 1E2, 1E4, 1E6]);
xticklabels({'$10^{-6}$', '$10^{-4}$', '$10^{-2}$', '$10^{0}$', '$10^{2}$', '$10^{4}$', '$10^{6}$'});
ytickangle(90);
set(gca,'xlim',sqrt(Q0*[t(1),t(end)])/lL);
set(gca,'ylim',[t_R(1),t_R(end)]);
saveas(figl,[postproc.folderloc1Dparent,'/RegimePlot_Dim_sqrtQ0tbylL_',Name],'fig');
close(figl);

% generating sqrt Q0t regime diagram
figl  =	figure('position',[50,21.8,800,731.2],'visible','on'); 
axl = axes(figl,'position',[0.14675,0.119116238449763,0.7,0.818285714285715]);
set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
set(gca,'linewidth',6); set(gca,'xscale','log'); set(gca,'yscale','log'); hold on;
ylabel('$t_R = \left(k_Rc_0\right)^{-1}$','interpreter','latex');
xlabel("$\bar{r}' = \sqrt{Q_0 t'}$",'interpreter','latex');
xx = [sqrt(Q0*t(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*ttrP(1))]; yy = [t_R(1), ttrP(1), t_R(1)]; 
patch(xx, yy, [1.00,0.70,0.70],'linestyle','none');	
xx = [sqrt(Q0*t(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*t(1))]; yy = [t_R(1), ttrP(1), t_R(end), t_R(end)]; 
patch(xx, yy, [0.80,0.50,0.50],'linestyle','none');
xx = [sqrt(Q0*ttrP(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*t(end)), sqrt(Q0*t(end))]; yy = [t_R(1), ttrP(1), t_R(end), t_R(1)]; 
patch(xx, yy, [0.70,0.70,1.00],'linestyle','none');
xx = [sqrt(Q0*ttrP(1)), sqrt(Q0*t(end)), sqrt(Q0*ttrP(1))]; yy = [ttrP(1), t_R(end), t_R(end)]; 
patch(xx, yy, [0.50,0.50,0.80],'linestyle','none');
plot(axl,sqrt(Q0*ttrL),t_R,':','color',[0.0,0.0,0.0],'linewidth',2);
plot(axl,sqrt(Q0*ttrP),t_R,'-.','color',[0.0,0.0,0.0],'linewidth',2);
hatchht = 10.^linspace(log10(t_R(1)),log10(t_R(end)),nhatch);
for ihatch = 1:nhatch
	plot([sqrt(Q0*mean(ttrL)),sqrt(Q0*mean(ttrP))],hatchht(ihatch)*[1,1],'k-','linewidth',1);
end
plot(axl,R,t_R,'-','color',[0.0,0.0,0.0],'linewidth',2);
yticks([1E-9, 1E-6, 1E-3, 1, 3600*24, 3600*24*365*100, 3600*24*365*1E6]);
yticklabels({'ns', '$\mu$s','ms', 's', 'day', 'cent', 'Myr'});
xticks([1E-6, 1E-4, 1E-2, 1E0, 1E2, 1E4, 1E6]);
xticklabels({'1 $\mu$', '0.1 mm', '1 cm', '1 m', '0.1 km', '10 km', '1000 km'});
ytickangle(90);
set(gca,'xlim',sqrt(Q0*[t(1),t(end)]));
set(gca,'ylim',[t_R(1),t_R(end)]);
saveas(figl,[postproc.folderloc1Dparent,'/RegimePlot_Dim_sqrtQ0t_',Name],'fig');
close(figl);

% generating sqrt t regime diagram
figl  =	figure('position',[50,21.8,800,731.2],'visible','on'); 
axl = axes(figl,'position',[0.14675,0.119116238449763,0.7,0.818285714285715]);
set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
set(gca,'linewidth',6); set(gca,'xscale','log'); set(gca,'yscale','log'); hold on;
ylabel('$t_R = \left(k_Rc_0\right)^{-1}$','interpreter','latex');
xlabel("$t'$",'interpreter','latex');
xx = [(t(1)), (ttrP(1)), (ttrP(1))]; yy = [t_R(1), ttrP(1), t_R(1)]; 
patch(xx, yy, [1.00,0.70,0.70],'linestyle','none');	
xx = [(t(1)), (ttrP(1)), (ttrP(1)), (t(1))]; yy = [t_R(1), ttrP(1), t_R(end), t_R(end)]; 
patch(xx, yy, [0.80,0.50,0.50],'linestyle','none');
xx = [(ttrP(1)), (ttrP(1)), (t(end)), (t(end))]; yy = [t_R(1), ttrP(1), t_R(end), t_R(1)]; 
patch(xx, yy, [0.70,0.70,1.00],'linestyle','none');
xx = [(ttrP(1)), (t(end)), (ttrP(1))]; yy = [ttrP(1), t_R(end), t_R(end)]; 
patch(xx, yy, [0.50,0.50,0.80],'linestyle','none');
plot(axl,(ttrL),t_R,':','color',[0.0,0.0,0.0],'linewidth',2);
plot(axl,(ttrP),t_R,'-.','color',[0.0,0.0,0.0],'linewidth',2);
hatchht = 10.^linspace(log10(t_R(1)),log10(t_R(end)),nhatch);
for ihatch = 1:nhatch
	plot([(mean(ttrL)),(mean(ttrP))],hatchht(ihatch)*[1,1],'k-','linewidth',1);
end
plot(axl,t,t_R,'-','color',[0.0,0.0,0.0],'linewidth',2);
yticks([1E-9, 1E-6, 1E-3, 1, 60, 3600*24, 3600*24*365, 3600*24*365*1000, 3600*24*365*1E6]);
yticklabels({'ns', '$\mu$s','ms', 's', 'min', 'day', 'yr', 'mil', 'Myr'});
xticks([1E-9, 1E-6, 1E-3, 1, 60, 3600*24, 3600*24*365, 3600*24*365*1000, 3600*24*365*1E6]);
xticklabels({'ns', '$\mu$s','ms', 's', 'min', 'day', 'yr', 'mil', 'Myr'});
ytickangle(90);
set(gca,'xlim',([t(1),t(end)]));
set(gca,'ylim',[t_R(1),t_R(end)]);
saveas(figl,[postproc.folderloc1Dparent,'/RegimePlot_Dim_t_',Name],'fig');
close(figl);

% generating r_f regime diagram
figl  =	figure('position',[50,21.8,800,731.2],'visible','on'); 
axl = axes(figl,'position',[0.14675,0.119116238449763,0.7,0.818285714285715]);
set(axl,'box','on'); grid('off'); set(axl,'fontsize',25); set(axl,'ticklabelinterpreter','latex'); 
set(gca,'linewidth',6); set(gca,'xscale','log'); set(gca,'yscale','log'); hold on;
ylabel('$t_R = \left(k_Rc_0\right)^{-1}$','interpreter','latex');
xlabel("${r_f}'$",'interpreter','latex');
%{
xx = [sqrt(Q0*t(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*ttrP(1))]; yy = [t_R(1), ttrP(1), t_R(1)]; 
patch(xx, yy, [1.00,0.70,0.70],'linestyle','none');	
xx = [sqrt(Q0*t(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*t(1))]; yy = [t_R(1), ttrP(1), t_R(end), t_R(end)]; 
patch(xx, yy, [0.80,0.50,0.50],'linestyle','none');
xx = [sqrt(Q0*ttrP(1)), sqrt(Q0*ttrP(1)), sqrt(Q0*t(end)), sqrt(Q0*t(end))]; yy = [t_R(1), ttrP(1), t_R(end), t_R(1)]; 
patch(xx, yy, [0.70,0.70,1.00],'linestyle','none');
xx = [sqrt(Q0*ttrP(1)), sqrt(Q0*t(end)), sqrt(Q0*ttrP(1))]; yy = [ttrP(1), t_R(end), t_R(end)]; 
patch(xx, yy, [0.50,0.50,0.80],'linestyle','none');
%}
plot(axl,eta,t_R,':','color',[0.0,0.0,0.0],'linewidth',2);
plot(axl,eta*Pe,t_R,'-.','color',[0.0,0.0,0.0],'linewidth',2);
%{
hatchht = 10.^linspace(log10(t_R(1)),log10(t_R(end)),nhatch);
for ihatch = 1:nhatch
	plot([sqrt(Q0*mean(ttrL)),sqrt(Q0*mean(ttrP))],hatchht(ihatch)*[1,1],'k-','linewidth',1);
end
%}
rfnumend =		load('rf_for_CI_Dim.mat');
indPe =			1+double(Pe>1000);
if (Pe > 1000)
	indL = 45;
	indR = 69;
else
	indL = 15;
	indR = 69;
end

firstpart =			r_f_Diff(1:indL-1);
secondpart =		pchip(rfnumend.eta(1:24),rfnumend.rf_end_num(indPe,1:24),eta);
secondpart =		secondpart(indL:indR).*sqrt(Q0*t_R(indL:indR));
thirdpart =			r_f_Disp(indR+1:end);
r_f_tR =			[firstpart,secondpart,thirdpart];
plot(axl,r_f_tR,t_R,'-','color',[0.0,0.0,0.0],'linewidth',2);
yticks([1E-9, 1E-6, 1E-3, 1, 3600*24, 3600*24*365*100, 3600*24*365*1E6]);
yticklabels({'ns', '$\mu$s','ms', 's', 'day', 'cent', 'Myr'});
%{
xticks([1E-10, 1E-6, 1E-2, 1E2, 1E6]);
xticklabels({'1 \angstrom', '1 $\mu$', '1 cm', '0.1 km', '1000 km'});
%}
xticks([1E-8, 1E-5, 1E-2, 1E1, 1E4, 1E8, 1E12]);
xticklabels({'10 nm', '10 $\mu$m', '1 cm', '10 m', '$10$ km',  '$10^5$ km', '$1$ AU'});
ytickangle(90);
% set(gca,'xlim',sqrt(Q0*[t(1),t(end)]));
set(gca,'ylim',[t_R(1),t_R(end)]);
saveas(figl,[postproc.folderloc1Dparent,'/RegimePlot_Dim_rf_',Name],'fig');
% close(figl);

%{
ticks([1E-9, 1E-6, 1E-3, 1, 60, 3600, 3600*24, 3600*24*365, 3600*24*365*100, 3600*24*365*1000*10, 3600*24*365*1E6]);
ticklabels({'ns', '$\mu$s','ms', 's', 'min', 'hr', 'day', 'yr', 'cent', 'myria', 'Myr'});
%}