% tests the front-attached frame solution for phi_A corresponding to the similarity solution 
% [cylindrical front -> no dispersion scenario & dispersion negligible regime]

function solnS1 =		solver_S1(simul,fltr)

%	obtaining computational parameters
	nx =				simul.simil.nx;
	absmaxx =			simul.simil.absmaxx;
	power =				simul.simil.power;
	gamm =				fltr.nondim.gamm;
	Pe =				fltr.nondim.Pe;
	derivwarn =			simul.derivwarn;
	tolrootfind =		simul.simil.tolrootfind;
	rootfindrelax =		simul.simil.rootfindrelax;
	leftBCtuned =		simul.simil.leftBCtuned;
	rootfindmeth =		simul.simil.rootfindmeth;
	maxrootfinditer =	simul.simil.maxrootfinditer;
	barQthres =			simul.simil.barQthres;

	barQ =				Pe/2;
	if (barQ > barQthres)
		barQ =			1/2;
		z1 =			sqrt(barQ+erfinv((1-gamm)/(1+gamm)));
		barq =			2*z1;
		if (gamm ~= 1)
			K =			((1+gamm)/sqrt(pi))*(exp(-(erfinv((1-gamm-eps)/(1+gamm+eps)))^2)*	...
						(z1/erfinv((1-gamm-eps)/(1+gamm+eps)))*(1-(barQ/(z1^2))));
		else
			K =			1.0;
		end
	else
		barq =			2*sqrt(gammaincinv(gamm/(gamm+1),barQ,'upper'));
		K =				(1+gamm)*(barq/2)*(1/gamma(barQ))*((((barq^2)/4)^(barQ-1))*exp(-((barq^2)/4)));
	end

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
			acoeffx2 =		differentcoeff(x,ix,2,'cd',derivwarn);
			for jx = -1:1
				row(icounter) =		ix;
				col(icounter) =		ix+jx;
				jaco(icounter) =	acoeffx2(3+jx)-Pe*(K*x(ix)+(2-double(strcmp(rootfindmeth,'STL')==1))*G(ix))*double(jx==0);
				icounter =			icounter+1;
			end
			if (strcmp(rootfindmeth,'STL') == 0)
				B(ix) =				differentiation(G,x,ix,2,'cd',derivwarn)-Pe*G(ix)*(K*x(ix)+G(ix));
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
%		numerical solution of phiA
		if (strcmp(rootfindmeth,'STL') == 1)
			G =			(1-rootfindrelax)*G+(rootfindrelax*sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B)';	
			errrootfind=sqrt(sum((G-Gprev).^2));
		else
			G =			G-(rootfindrelax*sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B)';
			errrootfind=sqrt(sum((B).^2));
		end
		disp(['solver_S1 Root finding error is ',num2str(errrootfind)]);
	end
%	obtaining front width integral
	IntNr =				trapz(x,x.^2.*(G.^2+K*G.*x));
	IntDr =				trapz(x,G.^2+K*G.*x);

%	containing output in output data structure
	solnS1.x =								x;
	solnS1.G =								G;
	solnS1.R =								G.^2+K*G.*x;
	[solnS1.Rmax, solnS1.locRmax] =			max(solnS1.R);
	solnS1.K =								K;
	solnS1.barq =							barq;
	solnS1.baralpha =						sqrt(IntNr/IntDr);
	solnS1.baralphaapprox =					0.979197710209806*(1./sqrt(Pe))*(K/sqrt(Pe)).^(-1/3);
	solnS1.barbeta =						G(1+(nx-1)/2)^2;
	solnS1.barbetaapprox =					0.5454*(K/sqrt(Pe))^(2/3);
	solnS1.barzeta_0 =						(2*pi*barq*(1/Pe)*K);
	solnS1.barZ_0 =							solnS1.barzeta_0;
	
end