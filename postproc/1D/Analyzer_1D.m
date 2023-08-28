% loading solution and pre-setting differentiation directions
if ((isfile([postproc.folderloc1D,'/soln_postproc.mat']) == 1) && (1 == 0))
	postproc_original = postproc;
	load([postproc.folderloc1D,'/soln_postproc.mat']);
	postproc = postproc_original;
else
	soln1D =									load([postproc.folderloc1D,'/soln.mat']);
	soln1D.simul.compMassC =					0;
	soln1D.simul.compfront =					0;
	if ((soln1D.simul.compMassC == 1) || (soln1D.simul.compfront == 1))
		soln1D.solnreactive =					load([postproc.folderloc1D,'/solnreactive.mat']);	
		soln1D.solnreactive =					soln1D.solnreactive.solnreactive;
	end
end
% trimming step save array to unique values
soln1D.simul.tstepsave =						unique(soln1D.simul.tstepsave);
% obtaining analytical solution
if (postproc.obtainanalytical == 1)
	if ((isfield(soln1D,'solnsimil') == 0) || (postproc.forcesolnsimil == 1))
		soln1D.simul.simil =					postproc.simil;
		if (soln1D.geomdom.isradial == 0)
			soln1D.soln.flow.vel.Px =			soln1D.soln.flow.vel.v;
			soln1D.fltr.flowbc.v =				soln1D.soln.flow.vel.Px(1);
			soln1D.solnsimil =					solver_S0(soln1D.simul,soln1D.fltr);
			soln1D.solnsimil_ET =				soln1D.solnsimil;
		elseif (soln1D.geomdom.isradial == 1)
			soln1D.solnsimil =					solver_S1(soln1D.simul,soln1D.fltr);
			soln1D.solnsimil_ET =				solver_S1_ET(soln1D.simul,soln1D.fltr);
			if	((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
				soln1D.solnsimil_disp =			solver_S2(soln1D.simul,soln1D.fltr);
				soln1D.solnsimil_disp_ET =		solver_S2_ET(soln1D.simul,soln1D.fltr);
			end
		else
			soln1D.solnsimil =					solver_S3(soln1D.simul,soln1D.fltr);
			soln1D.solnsimil_ET =				solver_S4(soln1D.simul,soln1D.fltr);
			if	((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
				soln1D.solnsimil_disp =			soln1D.solnsimil;
				soln1D.solnsimil_disp_ET =		soln1D.solnsimil_ET;
				eta_orig =						soln1D.fltr.eta;
				soln1D.fltr.eta =				0.0;
				soln1D.solnsimil =				solver_S3(soln1D.simul,soln1D.fltr);
				soln1D.solnsimil_ET =			solver_S4(soln1D.simul,soln1D.fltr);
				soln1D.fltr.eta =				eta_orig;
			end
		end
	end
end
dirx    =										cell(1,soln1D.simul.nx);
dirx(:) =										{'cd'};
dirx(1) =										{'fd'};
dirx(soln1D.simul.nx) =							{'bd'};
dirt    =										cell(1,soln1D.simul.nt);
dirt(:) =										{'cd'};
dirt(1) =										{'fd'};
dirt(soln1D.simul.nt) =							{'bd'};
% resetting grid to pre-conformed space if solved in conformed space
try
	if (soln1D.simul.confmap.space == 1)
		soln1D.grd.x =							soln1D.grd.x./(1-soln1D.grd.x);
	end
	if (soln1D.simul.confmap.time == 1)
		soln1D.grd.t =							soln1D.grd.t./(1-soln1D.grd.t);
	end
catch
end

% loading comparison solutions
if (postproc.iscomparison2D == 1)
	solncompare2D =								load([postproc.folderloccompare2D,'/soln.mat']);
	solncompare2D.simul.tstepsave =				unique(solncompare2D.simul.tstepsave);
end
if (postproc.iscomparison1D == 1)
	solncompare1D =								load([postproc.folderloccompare1D,'/soln.mat']);
	solncompare1D.simul.tstepsave =				unique(solncompare1D.simul.tstepsave);
end

% obtaining time grid that is saved to files
front.tsaved =									soln1D.grd.t(soln1D.simul.tstepsave);

% obtaining advective front location
front.barq_adv =								double(soln1D.geomdom.isradial~=0)*	...
												(double(soln1D.geomdom.isradial+1)^(1/double(soln1D.geomdom.isradial+1)));
front.r_adv =									front.barq_adv*soln1D.grd.t.^(1/double(soln1D.geomdom.isradial+1));

% obtaining analytical solution
if (postproc.obtainanalytical == 1)
% checking if analytical solution has to be obtained
	if ((exist('front','var') == 0) || (postproc.forceanalytical == 1))
		needanasoln =							1;
	else
		if (isfield(front,'ana') == 0)
			needanasoln =						1;
		else
			needanasoln =						0;
		end
	end
	if (needanasoln == 1)
%		insights from analytical approximations
		if (soln1D.geomdom.isradial == 0)
%			front properties
			front.ana.r_f =						0.5*soln1D.geomdom.size.Lx+soln1D.solnsimil.barq*sqrt(soln1D.grd.t);
			front.ana.w_f =						soln1D.solnsimil.baralpha*(soln1D.grd.t.^(1/6));
			front.ana.w_f_approx =				soln1D.solnsimil.baralphaapprox*(soln1D.grd.t.^(1/6));
			front.ana.R_f =						soln1D.solnsimil.barbeta*(soln1D.grd.t.^(-2/3));
			front.ana.R_f_approx =				soln1D.solnsimil.barbetaapprox*(soln1D.grd.t.^(-2/3));
			front.ana.R_max =					soln1D.solnsimil.Rmax*(soln1D.grd.t.^(-2/3));
			front.ana.R_max_loc =				soln1D.solnsimil.locRmax.*soln1D.grd.t.^0;
			front.ana.Mc =						soln1D.solnsimil.barZ_0*soln1D.grd.t.^(1/2)+		...
												soln1D.solnsimil.barZ_1*soln1D.grd.t.^(1/2);
			front.ana.dMcdt =					soln1D.solnsimil.barzeta_0*soln1D.grd.t.^(-1/2)+	...
												soln1D.solnsimil.barzeta_1*soln1D.grd.t.^(-1/2);
			front.ana.w_f_ET =					soln1D.solnsimil_ET.baralpha*(soln1D.grd.t.^(1/6));
			front.ana.w_f_approx_ET =			soln1D.solnsimil_ET.baralphaapprox*(soln1D.grd.t.^(1/6));
			front.ana.R_f_ET =					soln1D.solnsimil_ET.barbeta*(soln1D.grd.t.^(-2/3));
			front.ana.R_f_approx_ET =			soln1D.solnsimil_ET.barbetaapprox*(soln1D.grd.t.^(-2/3));
			front.ana.R_max_ET =				soln1D.solnsimil_ET.Rmax*(soln1D.grd.t.^(-2/3));
			front.ana.R_max_loc_ET =			soln1D.solnsimil_ET.locRmax.*soln1D.grd.t.^0;
			front.ana.Mc_ET =					soln1D.solnsimil_ET.barZ_0*soln1D.grd.t.^(1/2)+		...
												soln1D.solnsimil_ET.barZ_1*soln1D.grd.t.^(1/2);
			front.ana.dMcdt_ET =				soln1D.solnsimil_ET.barzeta_0*soln1D.grd.t.^(-1/2)+	...
												soln1D.solnsimil_ET.barzeta_1*soln1D.grd.t.^(-1/2);
%			similarity solutions
			front.ana.z_simil =					soln1D.solnsimil.x;
			front.ana.G_simil =					soln1D.solnsimil.G;
			front.ana.z_simil_ET =				soln1D.solnsimil_ET.x;
			front.ana.G_simil_ET =				soln1D.solnsimil_ET.G;
		elseif (soln1D.geomdom.isradial == 1)
%			front properties
			front.ana.r_f =						soln1D.solnsimil.barq*sqrt(soln1D.grd.t);
			front.ana.w_f =						soln1D.solnsimil.baralpha*(soln1D.grd.t.^(1/6));
			front.ana.w_f_approx =				soln1D.solnsimil.baralphaapprox*(soln1D.grd.t.^(1/6));
			front.ana.R_f =						soln1D.solnsimil.barbeta*(soln1D.grd.t.^(-2/3));
			front.ana.R_f_approx =				soln1D.solnsimil.barbetaapprox*(soln1D.grd.t.^(-2/3));
			front.ana.R_max =					soln1D.solnsimil.Rmax*(soln1D.grd.t.^(-2/3));
			front.ana.R_max_loc =				soln1D.solnsimil.locRmax.*soln1D.grd.t.^0;
			front.ana.Mc =						soln1D.solnsimil.barZ_0*soln1D.grd.t;
			front.ana.dMcdt =					soln1D.solnsimil.barzeta_0*soln1D.grd.t.^0;
			front.ana.w_f_ET =					soln1D.solnsimil_ET.baralpha*(soln1D.grd.t.^(1/2));
			front.ana.w_f_approx_ET =			soln1D.solnsimil_ET.baralphaapprox*(soln1D.grd.t.^(1/2));
			front.ana.R_f_ET =					soln1D.solnsimil_ET.barbeta*(soln1D.grd.t.^0);
			front.ana.R_f_approx_ET =			soln1D.solnsimil_ET.barbetaapprox*(soln1D.grd.t.^0);
			front.ana.R_max_ET =				soln1D.solnsimil_ET.Rmax*(soln1D.grd.t.^(0));
			front.ana.R_max_loc_ET =			soln1D.solnsimil_ET.locRmax;
			front.ana.Mc_ET =					soln1D.solnsimil_ET.barZ_0*soln1D.grd.t.^2;
			front.ana.dMcdt_ET =				soln1D.solnsimil_ET.barzeta_0*soln1D.grd.t;
			if	((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
				front.ana.r_f_disp =			soln1D.solnsimil_disp.barq*(soln1D.grd.t).^(1/3);
				front.ana.w_f_disp =			soln1D.solnsimil_disp.baralpha*(soln1D.grd.t).^(0);
				front.ana.w_f_disp_approx =		soln1D.solnsimil_disp.baralphaapprox*(soln1D.grd.t).^(0);
				front.ana.R_f_disp =			soln1D.solnsimil_disp.barbeta*(soln1D.grd.t.^(-2/3));
				front.ana.R_f_disp_approx =		soln1D.solnsimil_disp.barbetaapprox*(soln1D.grd.t.^(-2/3));
				front.ana.R_max_disp =			soln1D.solnsimil_disp.Rmax*(soln1D.grd.t.^(-2/3));
				front.ana.R_max_loc_disp =		soln1D.solnsimil_disp.locRmax.*soln1D.grd.t.^0;
				front.ana.Mc_disp =				soln1D.solnsimil_disp.barZ_0*soln1D.grd.t.^(2/3);
				front.ana.dMcdt_disp =			soln1D.solnsimil_disp.barzeta_0*soln1D.grd.t.^(-1/3);
				front.ana.w_f_disp_ET =			soln1D.solnsimil_disp.baralpha*(soln1D.grd.t).^(1/3);
				front.ana.w_f_disp_approx_ET =	soln1D.solnsimil_disp_ET.baralphaapprox*(soln1D.grd.t).^(1/3);
				front.ana.R_f_disp_ET =			soln1D.solnsimil_disp_ET.barbeta*(soln1D.grd.t.^0);
				front.ana.R_f_disp_approx_ET =	soln1D.solnsimil_disp_ET.barbetaapprox*(soln1D.grd.t.^0);
				front.ana.R_max_disp_ET =		soln1D.solnsimil_disp_ET.Rmax*(soln1D.grd.t.^0);
				front.ana.R_max_loc_disp_ET =	soln1D.solnsimil_disp_ET.locRmax;
				front.ana.Mc_disp_ET =			soln1D.solnsimil_disp_ET.barZ_0*soln1D.grd.t.^(5/3);
				front.ana.dMcdt_disp_ET =		soln1D.solnsimil_disp_ET.barzeta_0*soln1D.grd.t.^(2/3);
			end
%			similarity solutions
			front.ana.z_simil =					soln1D.solnsimil.x;
			front.ana.G_simil =					soln1D.solnsimil.G;
			front.ana.z_simil_ET =				soln1D.solnsimil_ET.x;
			front.ana.G_simil_ET =				soln1D.solnsimil_ET.G;
			if	((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
				front.ana.z_simil_disp =		soln1D.solnsimil_disp.x;
				front.ana.G_simil_disp =		soln1D.solnsimil_disp.G;
				front.ana.z_simil_disp_ET =		soln1D.solnsimil_disp_ET.x;
				front.ana.G_simil_disp_ET =		soln1D.solnsimil_disp_ET.G;
			end
		else
%			front properties
			front.ana.r_f =						soln1D.solnsimil.r_fS*soln1D.grd.t.^0;
			front.ana.w_f =						soln1D.solnsimil.w_fS*soln1D.grd.t.^0;
			front.ana.w_f_approx =				soln1D.solnsimil.w_fS_approx*soln1D.grd.t.^0;
			front.ana.R_f =						soln1D.solnsimil.RFS*soln1D.grd.t.^0;
			front.ana.R_f_approx =				soln1D.solnsimil.RFS*soln1D.grd.t.^0;
			front.ana.R_max =					soln1D.solnsimil.RFS*soln1D.grd.t.^0;
			front.ana.R_max_loc =				soln1D.solnsimil.locRmax.*soln1D.grd.t.^0;
			front.ana.Mc =						soln1D.solnsimil.barZ_0*soln1D.grd.t;
			front.ana.dMcdt =					soln1D.solnsimil.barzeta_0*soln1D.grd.t.^0;
			front.ana.r_f_ET =					front.r_adv;
			front.ana.w_f_ET =					soln1D.solnsimil_ET.baralpha*(soln1D.grd.t.^(1/2));
			front.ana.w_f_approx_ET =			soln1D.solnsimil_ET.baralpha_half*(soln1D.grd.t.^(1/2));
			front.ana.R_f_ET =					soln1D.solnsimil_ET.Rf.*soln1D.grd.t.^0;
			front.ana.R_f_approx_ET =			soln1D.solnsimil_ET.Rf.*soln1D.grd.t.^0;
			front.ana.R_max_ET =				soln1D.solnsimil_ET.Rmax.*soln1D.grd.t.^0;
			front.ana.R_max_loc_ET =			soln1D.solnsimil_ET.locRmax.*soln1D.grd.t.^0;
			front.ana.Mc_ET =					soln1D.solnsimil_ET.barZ_0*soln1D.grd.t.^(5/3) +	...
												soln1D.solnsimil_ET.barZ_1*soln1D.grd.t.^2;
			front.ana.dMcdt_ET =				soln1D.solnsimil_ET.barzeta_0*soln1D.grd.t.^(2/3) +	...
												soln1D.solnsimil_ET.barzeta_1*soln1D.grd.t;
			front.ana.z_simil =					soln1D.solnsimil.x;
			front.ana.G_simil =					soln1D.solnsimil.G;
			front.ana.z_simil_ET =				soln1D.solnsimil_ET.x;
			front.ana.G_simil_ET =				soln1D.solnsimil_ET.G;
			if (soln1D.fltr.eta == 0)
			else
				front.ana.r_f_ET =				front.r_adv.^(double(soln1D.fltr.eta~=0)+1);
				front.ana.w_f_ET =				...
				soln1D.solnsimil_ET.baralpha*(soln1D.grd.t.^(1/(2*(double(soln1D.fltr.eta~=0)+1))));
				front.ana.w_f_approx_ET =		...
				soln1D.solnsimil_ET.baralpha*(soln1D.grd.t.^(1/(2*(double(soln1D.fltr.eta~=0)+1))));
			end
			if (soln1D.fltr.eta == 0)
			else
				front.ana.Mc_ET =				soln1D.solnsimil_ET.barZ_0*soln1D.grd.t.^(5/3) +	...
												soln1D.solnsimil_ET.barZ_1*soln1D.grd.t.^(3/2);
				front.ana.dMcdt_ET =			soln1D.solnsimil_ET.barzeta_0*soln1D.grd.t.^(2/3) +	...
												soln1D.solnsimil_ET.barzeta_1*soln1D.grd.t.^(1/2);
				front.ana.z_simil =				soln1D.solnsimil.x;
				front.ana.G_simil =				soln1D.solnsimil.G;
				front.ana.z_simil_ET =			soln1D.solnsimil_ET.x;
				front.ana.G_simil_ET =			soln1D.solnsimil_ET.G;
			end
			if	((soln1D.fltr.eta ~= 0) && ...
				((soln1D.fltr.nondim.Pe <= postproc.thresPedispsimil) && (soln1D.fltr.nondim.Pe >= 1/postproc.thresPedispsimil)))
					front.ana.r_f_disp =				soln1D.solnsimil_disp.r_fS*soln1D.grd.t.^0;
					front.ana.w_f_disp =				soln1D.solnsimil_disp.w_fS*soln1D.grd.t.^0;
					front.ana.w_f_approx_disp =			soln1D.solnsimil_disp.w_fS_approx*soln1D.grd.t.^0;
					front.ana.R_f_disp =				soln1D.solnsimil_disp.RFS*soln1D.grd.t.^0;
					front.ana.R_f_approx_disp =			soln1D.solnsimil_disp.RFS*soln1D.grd.t.^0;
					front.ana.R_max_disp =				soln1D.solnsimil_disp.RFS*soln1D.grd.t.^0;
					front.ana.R_max_loc_disp =			soln1D.solnsimil_disp.locRmax.*soln1D.grd.t.^0;
					front.ana.Mc_disp =					soln1D.solnsimil_disp.barZ_0*soln1D.grd.t;
					front.ana.dMcdt_disp =				soln1D.solnsimil_disp.barzeta_0*soln1D.grd.t.^0;
					front.ana.r_f_disp_ET =				front.r_adv;
					front.ana.w_f_disp_ET =				soln1D.solnsimil_disp_ET.baralpha*(soln1D.grd.t.^(1/4));
					front.ana.w_f_approx_disp_ET =		soln1D.solnsimil_disp_ET.baralpha_half*(soln1D.grd.t.^(1/4));
					front.ana.R_f_disp_ET =				soln1D.solnsimil_disp_ET.Rf.*soln1D.grd.t.^0;
					front.ana.R_f_approx_disp_ET =		soln1D.solnsimil_disp_ET.Rf.*soln1D.grd.t.^0;
					front.ana.R_max_disp_ET =			soln1D.solnsimil_disp_ET.Rmax.*soln1D.grd.t.^0;
					front.ana.R_max_loc_disp_ET =		soln1D.solnsimil_disp_ET.locRmax.*soln1D.grd.t.^0;
					front.ana.Mc_disp_ET =				soln1D.solnsimil_disp_ET.barZ_0*soln1D.grd.t.^(5/3) +		...
														soln1D.solnsimil_disp_ET.barZ_1*soln1D.grd.t.^(3/2);
					front.ana.dMcdt_ET =				soln1D.solnsimil_disp_ET.barzeta_0*soln1D.grd.t.^(2/3) +	...
														soln1D.solnsimil_disp_ET.barzeta_1*soln1D.grd.t.^(1/2);
					front.ana.z_simil_disp =			soln1D.solnsimil_disp.x;
					front.ana.G_simil_disp =			soln1D.solnsimil_disp.G;
					front.ana.z_simil_disp_ET =			soln1D.solnsimil_disp_ET.x;
					front.ana.G_simil_disp_ET =			soln1D.solnsimil_disp_ET.G;
			end
		end
	end
end
% obtaining numerical reaction rate, reaction zone width
if (postproc.iFront == 1)
	if ((exist('front','var') == 0) || (postproc.forcenumfront == 1))
		neednumfrontsoln =						1;
	else
		if (isfield(front,'num') == 0)
			neednumfrontsoln =					1;
		else
			neednumfrontsoln =					0;
		end
	end
	if (neednumfrontsoln == 1)
		if ((soln1D.simul.compfront == 0) || ((postproc.forcenumfront == 1) && (soln1D.simul.compfront == 1)))
			ifront.tsaved =						0;
			for it = soln1D.simul.tstepsave
				ifront.tsaved =					ifront.tsaved+1;
				phiCcurrent =					load([postproc.folderloc1D,'/phiC_',num2str(it),'.dat']);
				phiAcurrent =					load([postproc.folderloc1D,'/phiA_',num2str(it),'.dat']);
				phiBcurrent =					load([postproc.folderloc1D,'/phiB_',num2str(it),'.dat']);
				front.num.R =					phiAcurrent.*phiBcurrent;
				front.num.theta =				phiAcurrent-phiBcurrent;
				[front.num.R_f(ifront.tsaved),front.num.inodrf(ifront.tsaved)] =		...
												max(front.num.R);
				[~,front.num.inodrf_orig(ifront.tsaved)] =								...
												min(front.num.theta.^2);
				front.num.r_f(ifront.tsaved)=	soln1D.grd.x(front.num.inodrf(ifront.tsaved));
				front.num.r_f_orig(ifront.tsaved) =			...
												soln1D.grd.x(front.num.inodrf_orig(ifront.tsaved));
				front.num.R_f_orig(ifront.tsaved) =			...
												front.num.R(front.num.inodrf_orig(ifront.tsaved));
				for iwidthdef = 1:2
					for iwidthbasefunc = 1:2
						IntNr =					trapz(soln1D.grd.x,		...
												(soln1D.grd.x-front.num.r_f_orig(ifront.tsaved)).^2.*		...
												front.num.R.*...
												(soln1D.grd.x.^(double((iwidthdef-1)*soln1D.geomdom.isradial))));
						IntDr =					trapz(soln1D.grd.x,front.num.R.*...
												(soln1D.grd.x.^(double((iwidthdef-1)*soln1D.geomdom.isradial))));
						front.num.w_f_orig(ifront.tsaved,iwidthdef,iwidthbasefunc) =						...
												sqrt(IntNr/IntDr);
					end
				end
				leftswitch =					1;
				rightswitch =					1;
				ixleft =						1;
				ixright =						soln1D.simul.nx;
				for ix = 1:soln1D.simul.nx
					if ((front.num.R(ix) > 0.5*(front.num.R_f(ifront.tsaved))) && (leftswitch == 1))
						ixleft =				ix;
						leftswitch =			0;
					end
					if ((front.num.R(ix) < 0.5*(front.num.R_f(ifront.tsaved))) && (leftswitch == 0) && (rightswitch == 1))
						ixright =				ix-1;
						rightswitch =			0;
					end
				end
				front.num.w_f(ifront.tsaved) =	soln1D.grd.x(ixright)-soln1D.grd.x(ixleft);	
			end
		end
	end
end

% obtaining transition zones
front.trnszon =									ones(2,length(postproc.fracTransitn),2,2,2);
if (soln1D.simul.compfront == 1)
	for ifracTransitn = 1:length(postproc.fracTransitn)
		[~,front.trnszon(1,ifracTransitn,1,1,1)]=	min((soln1D.solnreactive.xforig-	...
											((postproc.fracTransitn(ifracTransitn)*soln1D.fltr.eta)^		...
											(1/double(soln1D.geomdom.isradial)))).^2);
		[~,front.trnszon(1,ifracTransitn,1,2,1)]=	min((soln1D.solnreactive.xforig-	...
											(((1/postproc.fracTransitn(ifracTransitn))*soln1D.fltr.eta)^	...
											(1/double(soln1D.geomdom.isradial)))).^2);
		[~,front.trnszon(1,ifracTransitn,2,1,1)]=	min((soln1D.solnreactive.xforig-	...
											((postproc.fracTransitn(ifracTransitn)*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^		...
											(1/double(soln1D.geomdom.isradial)))).^2);
		[~,front.trnszon(1,ifracTransitn,2,2,1)]=	min((soln1D.solnreactive.xforig-	...
											(((1/postproc.fracTransitn(ifracTransitn))*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^	...
											(1/double(soln1D.geomdom.isradial)))).^2);
	end
	for ilimdef = 1:2
		for iside = 1:2
			[~,front.trnszon(1,ifracTransitn,ilimdef,iside,2)] =	...
											min((front.tsaved-soln1D.grd.t(front.trnszon(1,ifracTransitn,ilimdef,iside,1))).^2);
		end
	end
end
if (postproc.iFront == 1)
	if ((soln1D.simul.compfront == 0) || ((postproc.forcenumfront == 1) && (soln1D.simul.compfront == 1)))
		for ifracTransitn = 1:length(postproc.fracTransitn)
			[~,front.trnszon(2,ifracTransitn,1,1,2)] =			...
											min((front.num.r_f-	...
											((postproc.fracTransitn(ifracTransitn)*soln1D.fltr.eta)^							...
											(1/double(soln1D.geomdom.isradial)))).^2);
			[~,front.trnszon(2,ifracTransitn,1,2,2)] =		...
											min((front.num.r_f-	...
											(((1/postproc.fracTransitn(ifracTransitn))*soln1D.fltr.eta)^						...
											(1/double(soln1D.geomdom.isradial)))).^2);
			[~,front.trnszon(2,ifracTransitn,2,1,2)] =		...
											min((front.num.r_f-	...
											((postproc.fracTransitn(ifracTransitn)*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^		...
											(1/double(soln1D.geomdom.isradial)))).^2);
			[~,front.trnszon(2,ifracTransitn,2,2,2)] =		...
											min((front.num.r_f-	...
											(((1/postproc.fracTransitn(ifracTransitn))*soln1D.fltr.eta*soln1D.fltr.nondim.Pe)^	...
											(1/double(soln1D.geomdom.isradial)))).^2);
			for ilimdef = 1:2
				for iside = 1:2
					front.trnszon(2,ifracTransitn,ilimdef,iside,1) =	...
											soln1D.simul.tstepsave(front.trnszon(2,ifracTransitn,ilimdef,iside,2));
				end
			end
		end
	end
end

% obtaining mass of product C
if (postproc.iMassC == 1)
	if ((exist('soln_analyze','var') == 0) || (postproc.forcemassC == 1))
		needmassCsoln =						1;
	else
		needmassCsoln =						0;
	end
	if (needmassCsoln == 1)
		soln_analyze.MassC =				zeros(length(soln1D.simul.tstepsave),1);
		soln_analyze.MassC_4mR =			zeros(length(soln1D.simul.tstepsave),1);
		soln_analyze.MassC_4mR_twk=			zeros(length(soln1D.simul.tstepsave),1);
		soln_analyze.barR =					zeros(length(soln1D.simul.tstepsave),1);
		ifront.tsaved = 0;
		for it = soln1D.simul.tstepsave
			ifront.tsaved =					ifront.tsaved+1;
			phiCcurrent =					load([postproc.folderloc1D,'/phiC_',num2str(it),'.dat']);
			phiBcurrent =					load([postproc.folderloc1D,'/phiB_',num2str(it),'.dat']);
			phiAcurrent =					load([postproc.folderloc1D,'/phiA_',num2str(it),'.dat']);
			soln_analyze.MassC(ifront.tsaved) =		...
											((double(soln1D.geomdom.isradial)*2*pi)^double(soln1D.geomdom.isradial~=0))*		...
											trapz(soln1D.grd.x,phiCcurrent.*(soln1D.grd.x.^soln1D.geomdom.isradial));
			soln_analyze.barR(ifront.tsaved) =		...
											((double(soln1D.geomdom.isradial)*2*pi)^double(soln1D.geomdom.isradial~=0))*		...
											trapz(soln1D.grd.x,phiAcurrent.*phiBcurrent.*(soln1D.grd.x.^soln1D.geomdom.isradial));
			if (ifront.tsaved > 1)
				soln_analyze.MassC_4mR(ifront.tsaved) =		...
											0.5*(soln_analyze.barR(ifront.tsaved)+soln_analyze.barR(ifront.tsaved-1))*	...
											(soln1D.grd.t(ifront.tsaved)-soln1D.grd.t(ifront.tsaved));
				soln_analyze.MassC_4mR_twk(ifront.tsaved) =		...
											soln_analyze.MassC_4mR(ifront.tsaved) +			...
											2*pi*phiCcurrent(1)*(soln1D.grd.t(ifront.tsaved)-soln1D.grd.t(ifront.tsaved));
			end
		end
	end
end

% obtaining terms of GDE at center and front
if (postproc.iTermsGDE == 1)
	if ((exist('dphidr_centre','var') == 0) || (exist('dphidr_f','var') == 0) ||				...
		(exist('dphidr_centre_full','var') == 0) || (exist('dphidr_f_full','var') == 0)	) % ||		...
	%	(postproc.forcetermsGDE == 1))
	%	initializing matrices
		dphidr_centre = zeros(length(soln1D.simul.tstepsave),4,2);
		dphidr_f = zeros(length(soln1D.simul.tstepsave),4,2);
		dphidr_centre_full = zeros(length(soln1D.simul.tstepsave),4,2);
		dphidr_f_full = zeros(length(soln1D.simul.tstepsave),4,2);
		multiplier_centre = zeros(length(soln1D.simul.tstepsave),1);
		multiplier_f = zeros(length(soln1D.simul.tstepsave),1);
	%	obtaining values for terms of GDE
		ifront.tsaved = 0;
		for it = soln1D.simul.tstepsave
			ifront.tsaved =								ifront.tsaved+1;
	%		accumulating numerical solutions
			phiCcurrent =		load([postproc.folderloc1D,'/phiC_',num2str(it),'.dat']);
			phiAcurrent =		load([postproc.folderloc1D,'/phiA_',num2str(it),'.dat']);
			phiBcurrent =		load([postproc.folderloc1D,'/phiB_',num2str(it),'.dat']);
			phiAmBcurrent =		phiAcurrent-phiBcurrent;
			for iorder = 1:2
				dphidr_centre(ifront.tsaved,1,iorder) =	differentiation(phiAcurrent,soln1D.grd.x,2,iorder,'cd',postproc.derivwarn);
				dphidr_centre(ifront.tsaved,2,iorder) =	differentiation(phiBcurrent,soln1D.grd.x,2,iorder,'cd',postproc.derivwarn);
				dphidr_centre(ifront.tsaved,3,iorder) =	differentiation(phiCcurrent,soln1D.grd.x,2,iorder,'cd',postproc.derivwarn);
				dphidr_centre(ifront.tsaved,4,iorder) =	...
				differentiation(phiAmBcurrent,soln1D.grd.x,2,iorder,'cd',postproc.derivwarn);
				multiplier_centre(ifront.tsaved) =	...
				(double(iorder==2)*(((soln1D.fltr.eta*soln1D.fltr.nondim.Pe)/soln1D.grd.x(1))+1)+...
				double(iorder==1)*((soln1D.fltr.nondim.Pe-1)/soln1D.grd.x(1)));
				dphidr_centre_full(ifront.tsaved,1,iorder) =	...
				multiplier_centre(ifront.tsaved)*dphidr_centre(ifront.tsaved,1,iorder);
				dphidr_centre_full(ifront.tsaved,2,iorder) =	...
				multiplier_centre(ifront.tsaved)*dphidr_centre(ifront.tsaved,2,iorder);
				dphidr_centre_full(ifront.tsaved,3,iorder) =	...
				multiplier_centre(ifront.tsaved)*dphidr_centre(ifront.tsaved,3,iorder);
				dphidr_centre_full(ifront.tsaved,4,iorder) =	...
				multiplier_centre(ifront.tsaved)*dphidr_centre(ifront.tsaved,4,iorder);
			end
			for iorder = 1:2
				inode = min([max([2,front.num.inodrf_orig(ifront.tsaved)]),soln1D.simul.nx]);
				dphidr_f(ifront.tsaved,1,iorder) =	differentiation(phiAcurrent,soln1D.grd.x,inode,iorder,'cd',postproc.derivwarn);
				dphidr_f(ifront.tsaved,2,iorder) =	differentiation(phiBcurrent,soln1D.grd.x,inode,iorder,'cd',postproc.derivwarn);
				dphidr_f(ifront.tsaved,3,iorder) =	differentiation(phiCcurrent,soln1D.grd.x,inode,iorder,'cd',postproc.derivwarn);
				dphidr_f(ifront.tsaved,4,iorder) =	...
				differentiation(phiAmBcurrent,soln1D.grd.x,inode,iorder,'cd',postproc.derivwarn);
				multiplier_f(ifront.tsaved) =	...
				(double(iorder==2)*(((soln1D.fltr.eta*soln1D.fltr.nondim.Pe)/soln1D.grd.x(inode))+1)+...
				double(iorder==1)*((soln1D.fltr.nondim.Pe-1)/soln1D.grd.x(inode)));
				dphidr_f_full(ifront.tsaved,1,iorder) =	multiplier_f(ifront.tsaved)*dphidr_f(ifront.tsaved,1,iorder);
				dphidr_f_full(ifront.tsaved,2,iorder) =	multiplier_f(ifront.tsaved)*dphidr_f(ifront.tsaved,2,iorder);
				dphidr_f_full(ifront.tsaved,3,iorder) =	multiplier_f(ifront.tsaved)*dphidr_f(ifront.tsaved,3,iorder);
				dphidr_f_full(ifront.tsaved,4,iorder) =	multiplier_f(ifront.tsaved)*dphidr_f(ifront.tsaved,4,iorder);
			end
		end
	end
end

% saving postprocessing output to file
igeom_orig = igeom; iPe_orig = iPe; ieta_orig = ieta; close all; clear postproc_original;
clear igeom iPe ieta;
save([postproc.folderloc1D,'/soln_postproc.mat']);
igeom = igeom_orig; iPe = iPe_orig; ieta = ieta_orig;