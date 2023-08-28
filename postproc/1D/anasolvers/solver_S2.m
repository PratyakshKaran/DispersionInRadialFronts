% tests the front-attached frame solution for phi_A corresponding to the similarity solution 
% [cylindrical front -> dispersion dominated regime]

function solnS2 =		solver_S2(simul,fltr)

%	obtaining computational parameters
	nx =				simul.simil.nx;
	absmaxx =			simul.simil.absmaxx;
	power =				simul.simil.power;
	eta =				fltr.eta;
	gamm =				fltr.nondim.gamm;
	Pe =				fltr.nondim.Pe;
	derivwarn =			simul.derivwarn;
	tolrootfind =		simul.simil.tolrootfind;
	rootfindrelax =		simul.simil.rootfindrelax;
	leftBCtuned =		simul.simil.leftBCtuned;
	rootfindmeth =		simul.simil.rootfindmeth;
	maxrootfinditer =	simul.simil.maxrootfinditer;

%	obtaining parameter values
	tildQ =				1/3;
	barq =				((9*eta)^(1/3))*	...
						((gammaincinv(gamm/(1+gamm),tildQ,'upper'))^(1/3));
% 	K =					(1+gamm)*(3/((9*eta)^(2/3)))*(1/gamma(tildQ))*(((barq^3)/(9*eta))^(tildQ-1)*exp(-((barq^3)/(9*eta))));
	K =					(1+gamm)*((barq^2)/(3*eta))*(1/gamma(tildQ))*((((barq^3)/(9*eta))^(tildQ-1))*exp(-((barq^3)/(9*eta))));
	tildPe =			(Pe*barq)/(eta*Pe);

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
			acoeffx1 =		differentcoeff(x,ix,1,'cd',derivwarn);
			acoeffx2 =		differentcoeff(x,ix,2,'cd',derivwarn);
			for jx = -1:1
				row(icounter) =		ix;
				col(icounter) =		ix+jx;
				jaco(icounter) =	acoeffx2(3+jx)+((1-Pe)/(eta*Pe))*acoeffx1(3+jx)-	...
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
		disp(['solver_S2 Root finding error is ',num2str(errrootfind)]);
	end

%	obtaining front width integral
	IntNr =				trapz(x,x.^2.*(G.^2+K*G.*x));
	IntDr =				trapz(x,G.^2+K*G.*x);

%	containing output in output data structure
	solnS2.x =								x;
	solnS2.G =								G;
	solnS2.R =								G.^2+K*G.*x;
	[solnS2.Rmax, solnS2.locRmax] =			max(solnS2.R);
	solnS2.K =								K;
	solnS2.barq =							barq;
	solnS2.baralpha =						sqrt(IntNr/IntDr);
	solnS2.baralphaapprox =					0.979197710209806*(1./sqrt(tildPe))*(K/sqrt(tildPe)).^(-1/3);
	solnS2.barbeta =						G(1+(nx-1)/2)^2;
	solnS2.barbetaapprox =					0.5454*(K/sqrt(tildPe))^(2/3);
	solnS2.barzeta_0 =						(2*pi*barq*(1/tildPe)*K);
	solnS2.barZ_0 =							solnS2.barzeta_0/(2/3);
	
end