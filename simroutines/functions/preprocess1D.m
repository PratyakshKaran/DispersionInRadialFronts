%-----------------------------------------------------------------------------------------------------------------------------------
% Pre-Processor for
% 1D Transient Reactive Transport Solver for Planar/Cylindrical/Spherical Injection using Finite Difference Method (FDM) &
% 1D Steady State Darcy Flow Solver corresponding for Planar/Cylindrical/Spherical Injection using Finite Difference Method (FDM)
% [Non-Function Program]
% 
%	written by --
%	uddipta ghosh (uddipta.ghosh@iitgn.ac.in)
%	pratyaksh karan (pratyakshkaran@gmail.com)
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Identifying Start Time and End Time
	starttime =						0.0;
	if ((simul.confmap.time == 0) || ((simul.confmap.time == 1) && (simul.confmap.timegridfromreal == 1)))
		endtime =					simul.dt*(simul.nt-1);
	else
		endtime =					1.0;
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Override Time Step Count to Satisfy Time Step Size Constraints	
	if (simul.limitdt == 1)
		simul.nt =					max([simul.nt,	...
									ceil(1+max([	...
									((endtime-starttime)/simul.dtstartmax)^(1/simul.tgridskew),	...
									1/(1.0-((1.0-(simul.dtendmax/(endtime-starttime)))^(1/simul.tgridskew)))	...
									]))]);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Time Grid Generation
	grd.t =							starttime+(endtime-starttime)*linspace(0.0,1.0,simul.nt).^simul.tgridskew;
	if ((simul.confmap.time == 1) && (simul.confmap.timegridfromreal == 1))
		grd.t =						grd.t./(1+grd.t);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Generating array of time steps to write solution to file
	if (strcmp(simul.t2savemode,'uniform'))
		simul.tstepsave =			floor(linspace(1,simul.nt,simul.nt2write));	
	elseif (strcmp(simul.t2savemode,'log'))
		simul.tstepsave =			floor(exp(linspace(log(1),log(simul.nt),simul.nt2write)));
	elseif (strcmp(simul.t2savemode,'revlog'))
		simul.tstepsave =			floor(exp(linspace(log(1),log(simul.nt),simul.nt2write)));
		simul.tstepsave =			simul.nt-simul.tstepsave;
		simul.tstepsave =			flip(simul.tstepsave);
	elseif (strcmp(simul.t2savemode,'bothlog'))
		simul.tstepsave =			floor(exp(linspace(log(1),log(simul.nt),floor(simul.nt2write/2))));
		simul.tstepsave1 =			floor(exp(linspace(log(1),log(simul.nt),ceil(simul.nt2write/2))));
		simul.tstepsave1 =			simul.nt-simul.tstepsave1;
		simul.tstepsave1 =			flip(simul.tstepsave1);
		simul.tstepsave = 			[simul.tstepsave,simul.tstepsave1];
		simul.tstepsave = 			sort(simul.tstepsave);
	elseif (strcmp(simul.t2savemode,'manual'))
		simul.tstepsave =			linspace(1,100,floor(simul.nt/3));
		simul.tstepsave =			[simul.tstepsave,linspace(101,1000,floor(simul.nt/3))];
		simul.tstepsave =			[simul.tstepsave,linspace(1001,10000,floor(simul.nt/3))];
		simul.nt2write =			length(simul.tstepsave);
	end
	if (simul.tstepsave(1) < 1)
		simul.tstepsave(1) =		1;
	end
	if (simul.tstepsave(end) > simul.nt)
		simul.tstepsave(end) =		simul.nt;
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Reset Geometry for Minimizing Computation
	if ((simul.compactgeom == 1) && ((simul.confmap.space == 0) || (geomdom.isradial == 0)))
		if (geomdom.isradial == 0)
			tildPe =				fltr.nondim.Pe/(1+fltr.eta*fltr.nondim.Pe*fltr.transpbc.frcdispers);
			alph =					(2/sqrt(tildPe))*erfinv((1-fltr.nondim.gamm)/(1+fltr.nondim.gamm));
			delt =					pi*sqrt(tildPe)*((1+fltr.nondim.gamm)/(2*sqrt(pi)))*exp(-0.25*(tildPe*alph)^2);
			geomdom.size.Lx =		min([geomdom.size.Lx,simul.compactgeomfactotal*(abs(alph)*sqrt(grd.t(end))+		...
									simul.compactgeomsubfacwidth*delt*grd.t(end)^(1/6))]);
		elseif (geomdom.isradial == 1)
			tildQ =					1/3;
			tildq =					((9*fltr.eta)^(1/3))*((gammaincinv(fltr.nondim.gamm/(1+fltr.nondim.gamm),tildQ,'upper'))^(1/3));
			tildK =					(1+fltr.nondim.gamm)*(tildq/2)*(1/gamma(tildQ))*(((tildq^2)/4)^(tildQ-1)*exp(-((tildq^2)/4)));
			barQ =					fltr.nondim.Pe/2;
			if (barQ > 10)
				z1 =				sqrt(2*barQ+erfinv((1-fltr.nondim.gamm)/(1+fltr.nondim.gamm)));
				barq =				2*z1;
				if (fltr.nondim.gamm ~= 1)
					barK =			((1+fltr.nondim.gamm)/sqrt(pi))*(exp(-(erfinv((1-fltr.nondim.gamm)/(1+fltr.nondim.gamm)))^2)*...
									(z1/erfinv((1-fltr.nondim.gamm)/(1+fltr.nondim.gamm)))*(1-(barQ/(z1^2))));
				else
					barK =			1.0;
				end
			else
				barq =				2*sqrt(gammaincinv(fltr.nondim.gamm/(fltr.nondim.gamm+1),barQ,'upper'));
				barK =				(1+fltr.nondim.gamm)*(barq/2)*(1/gamma(barQ))*(((barq^2)/4)^(barQ-1)*exp(-((barq^2)/4)));
			end
			baralpha =				0.979197710209806*(barK^(-1/3));
			tildalpha =				0.979197710209806*(tildK^(-1/3));
			geomdom.size.Lx =		min([geomdom.size.Lx,simul.compactgeomfactotal*	...
									max([tildq*(grd.t(end)).^(1/3)+simul.compactgeomsubfacwidth*tildalpha,	...
									barq*(grd.t(end)).^(1/2)+simul.compactgeomsubfacwidth*baralpha*(grd.t(end)).^(1/6)])]);
		else
			if (fltr.eta ~= 0)
				geomdom.size.Lx =	min([geomdom.size.Lx,simul.compactgeomfactotal*	...
									sqrt(fltr.eta*fltr.nondim.Pe)*tan(sqrt(fltr.eta/fltr.nondim.Pe)*		...
									log(1-(1/(1+fltr.nondim.gamm))*(1-exp((pi/2)*sqrt(fltr.nondim.Pe/fltr.eta)))))]);
			else
				geomdom.size.Lx =	min([geomdom.size.Lx,	...
									simul.compactgeomfactotal*(fltr.nondim.Pe/log(1+fltr.nondim.gamm))]);
			end
		end
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Override Neumann BC when conformally-mapped space is being used
	if (simul.confmap.space == 1)
		geomdom.transpbc.right =	[1,1,1];
		geomdom.flowbc.right =		1;
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Override Space Step Count to Satisfy Space Step Size Constraints	
	if (simul.limitdx == 1)
		if (geomdom.isradial == 0)
			geomdom.size.a =			0.0;
		else
			if ((simul.confmap.space == 0) || ((simul.confmap.space == 1) && (simul.confmap.timegridfromreal == 1)))
				geomdom.size.a =		fltr.transpbc.inject.k1;
			else
				geomdom.size.a =		fltr.transpbc.inject.k1/(1+fltr.transpbc.inject.k1);
			end
		end
		if ((simul.confmap.space == 0) || ((simul.confmap.space == 1) && (simul.confmap.timegridfromreal == 1)))
			geomdom.size.A =			geomdom.size.Lx;
		else
			geomdom.size.A =			1.0;
		end
		simul.dxRmax =					min([simul.dxRmax,simul.dxRmax_comprtoA*geomdom.size.A]);
		simul.nx =						max([simul.nx,	...
										ceil(1+max([	...
										((geomdom.size.A-geomdom.size.a)/simul.dxLmax)^(1/simul.xgridskew),	...
										1/(1.0-((1.0-(simul.dxRmax/(geomdom.size.A-geomdom.size.a)))^(1/simul.xgridskew)))	...
										]))]);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Computed Geometric & Domain Parameters
	if ((simul.confmap.space == 0) || ((simul.confmap.space == 1) && (simul.confmap.timegridfromreal == 1)))
		geomdom.size.a =				fltr.transpbc.inject.k1;
		geomdom.size.A =				geomdom.size.Lx;
	else
		geomdom.size.a =				fltr.transpbc.inject.k1/(1+fltr.transpbc.inject.k1);
		geomdom.size.A =				1.0;
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Setting End Pressure Values to obtain Desired Flow Rate
	if (simul.isflowrate == 1)
		if (geomdom.isradial == 0)
			fltr.flowbc.Hleft0 =	fltr.flowbc.Hright0+geomdom.size.Lx;
		elseif (geomdom.isradial == 1)
			fltr.flowbc.Hleft0 =	fltr.flowbc.Hright0-log(geomdom.size.a+1E-100)-log(geomdom.size.A);
		else
			fltr.flowbc.Hleft0 =	fltr.flowbc.Hright0+	...
									((geomdom.size.A-(geomdom.size.a+1E-100))/(geomdom.size.A*(geomdom.size.a+1E-100)));
		end
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Obtaining Absolute Front and Injection Parameters when provided values are Relative
if ((simul.frontparamsrel == 1) && (simul.confmap.space == 0))
 	fltr.init.front.k1 =			fltr.init.front.k1*geomdom.size.Lx;
 	fltr.init.front.k2 =			fltr.init.front.k2*geomdom.size.Lx;
end
if ((simul.injectparamsrel == 1) && (simul.confmap.space == 0))
 	fltr.transpbc.inject.k1 =		fltr.transpbc.inject.k1*geomdom.size.Lx;
end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Space Grid Generation
	if (geomdom.isradial == 0)
		leftend =					0.0;
		rightend =					geomdom.size.Lx;
	else
		if ((simul.confmap.space == 0) || ((simul.confmap.space == 1) && (simul.confmap.timegridfromreal == 1)))
			leftend =				geomdom.size.a;
			rightend =				geomdom.size.A;
		else
			leftend =				geomdom.size.a/(1+geomdom.size.a);
			rightend =				1.0;
		end
	end
	grd.x =							leftend+(rightend-leftend)*linspace(0.0,1.0,simul.nx).^simul.xgridskew;
	if ((simul.confmap.space == 1) && (simul.confmap.spacegridfromreal == 1))
		grd.x =						grd.x./(grd.x+1);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Permeability Generation
	fltr.Kperm =					fltr.perm.Kperm0*zeros(simul.nx,1);
	if (geomdom.perm == 0)			% uniform permeability
		fltr.Kperm =				fltr.perm.Kperm0*ones(simul.nx,1);
	elseif (geomdom.perm == 1)		% permeability based on a function
		for ix = 1:nx
			fltr.Kperm(ix) =		fltr.perm.Kperm0*(1.0+fltr.perm.Kvarrang*sin(4.0*pi*(grd.x(ix)/geomdom.size.Lx)));
		end
	elseif (geomdom.perm == 2)		% permeability generated using a stochasticity function
		fltr.Kperm =				permgen(fltr.perm,grd,geomdom);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Initial Reaction Front Generation
	fltr.init.phiA =				zeros(simul.nx,1);
	fltr.init.phiB =				zeros(simul.nx,1);
	fltr.init.phiC =				zeros(simul.nx,1);
	if (simul.justinject == 1)
		fltr.init.phiA(1) =			1.0;
		fltr.init.phiB(1:simul.nx)=	fltr.nondim.gamm;
	else
		fltr.init.phiA(1:simul.nx)=	0.5*(1.0-tanh((grd.x(1:simul.nx)-fltr.init.front.k1)/fltr.init.front.k2));
		fltr.init.phiB(1:simul.nx)=	fltr.nondim.gamm*0.5*(1.0+tanh((grd.x(1:simul.nx)-fltr.init.front.k1)/fltr.init.front.k2));
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	Creating Output Folder
	if (simul.localonly == 0)
		simul.outfoldrname =		['../../outdata/1D/',simul.casename];
	else
		simul.outfoldrname =		['../../../../LocalOnlyRawData_Level_',num2str(simul.localonlylevel),'/1D/',simul.casename];
	end
	mkdir('../../outdata/1D');
	mkdir(simul.outfoldrname);
%-----------------------------------------------------------------------------------------------------------------------------------
	