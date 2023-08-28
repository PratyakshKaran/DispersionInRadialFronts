% tests the front-attached frame solution for phi_A corresponding to the error-function solution
% [planar front]

function solnS0 =		solver_S0(simul,fltr)

%	obtaining computational parameters
	nx =				simul.simil.nx;
	absmaxx =			simul.simil.absmaxx;
	power =				simul.simil.power;
	eta =				fltr.eta;
	v =					fltr.flowbc.v;
	v_dispers =			fltr.transpbc.frcdispers;
	gamm =				fltr.nondim.gamm;
	Pe =				fltr.nondim.Pe;
	derivwarn =			simul.derivwarn;
	tolrootfind =		simul.simil.tolrootfind;
	rootfindrelax =		simul.simil.rootfindrelax;
	leftBCtuned =		simul.simil.leftBCtuned;
	rootfindmeth =		simul.simil.rootfindmeth;
	maxrootfinditer =	simul.simil.maxrootfinditer;

%	setting computation parameters
	tildPe =			Pe/(1+eta*Pe*(v+v_dispers));
	K =					sqrt(tildPe)*((1+gamm)/(2*sqrt(pi)))*exp(-(erfinv((1-gamm)/(1+gamm)))^2);
	barq =				(2/sqrt(tildPe))*erfinv((1-gamm)/(1+gamm));

%	generating x-grid
	xl =				flip(-absmaxx*linspace(0,1,(nx-1)/2+1).^power);
	xr =				absmaxx*linspace(0,1,(nx-1)/2+1).^power;
	x =					[xl,xr(2:end)];
	
%	initiating equation solver matrices
	row =				zeros(3*nx-4,1);
	col =				zeros(3*nx-4,1);
	jaco =				zeros(3*nx-4,1);
	B =					zeros(nx,1);
	
%	populating constant jacobian and residual terms
	icounter =			1;
	row(icounter) =		1;
	col(icounter) =		1;
	jaco(icounter) =	1.0;
	icounter =			icounter+1;
	row(icounter) =		nx;
	col(icounter) =		nx;
	jaco(icounter) =	1.0;
	icounter =			icounter+1;
	icounterresum =		icounter;
	if (strcmp(rootfindmeth,'STL') == 1)
		if (leftBCtuned == 1)
			B(1) =		-K*x(1);
		else
			B(1) =		1.0;
		end
	end

%	loop initial sanitization
	if (leftBCtuned == 1)
		G =				-K*x(1)+K*x(1)*((x-x(1))/(x(nx)-x(1)));
	else
		G =				1-((x-x(1))/(x(nx)-x(1)));
	end
	IntNr =		0.0;

%	populating loop-varying jacobian and residual entries
	errrootfind =		10*tolrootfind;
	iter = 0;
	while ((errrootfind > tolrootfind) && (iter <= maxrootfinditer))
		iter = iter+1;
		if (strcmp(rootfindmeth,'STL') == 1)
			Gprev =			G;
		end
		icounter =			icounterresum;
		for ix = 2:nx-1
			acoeffx2 =		differentcoeff(x,ix,2,'cd',derivwarn);
			for jx = -1:1
				row(icounter) =		ix;
				col(icounter) =		ix+jx;
				jaco(icounter) =	acoeffx2(3+jx)-	...
									tildPe*(K*x(ix)+(2-double(strcmp(rootfindmeth,'STL')==1))*G(ix))*double(jx==0);
				icounter =			icounter+1;
			end
			if (strcmp(rootfindmeth,'STL') == 0)
				B(ix) =				differentiation(G,x,ix,2,'cd',derivwarn)-tildPe*G(ix)*(K*x(ix)+G(ix));
			end
		end
		icounter =			icounter-1;
		if (strcmp(rootfindmeth,'STL') == 0)
			if (leftBCtuned == 1)
				B(1) =		G(1)+K*x(1);
			else
				B(1) =		G(1)-1.0;
			end
			B(nx) =			G(nx);
		end

%		numerical solution of G
		if (strcmp(rootfindmeth,'STL') == 1)
			G =			(1-rootfindrelax)*G+(rootfindrelax*sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B)';	
			errrootfind=sqrt(sum((G-Gprev).^2));
		else
			G =			G-(rootfindrelax*sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B)';
			errrootfind=sqrt(sum((B).^2));
		end
		disp(['solver_S0 Root finding error is ',num2str(errrootfind)]);
	end

%	obtaining front width integral
	IntNr =				trapz(x,x.^2.*(G.^2+K*G.*x));
	IntDr =				trapz(x,G.^2+K*G.*x);

%	containing output in output data structure
	solnS0.x =								x;
	solnS0.G =								G;
	solnS0.R =								(G.^2+K*G.*x);
	[solnS0.Rmax, solnS0.locRmax] =			max(solnS0.R);
	solnS0.K =								K;
	solnS0.barq =							barq;
	solnS0.baralpha =						sqrt(IntNr/IntDr);
	solnS0.baralphaapprox =					0.979197710209806*(1./sqrt(tildPe))*(K/sqrt(tildPe)).^(-1/3);
	solnS0.barbeta =						G(1+(nx-1)/2)^2;
	solnS0.barbetaapprox =					0.5454*(K/sqrt(tildPe))^(2/3);
	solnS0.barzeta_0 =						K/Pe;
	solnS0.barzeta_1 =						0;
	solnS0.barZ_0 =							3*solnS0.barzeta_0;
	solnS0.barZ_1 =							0;

end