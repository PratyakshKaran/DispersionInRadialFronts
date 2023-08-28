% tests the front-attached frame solution for phi_A corresponding to the similarity solution 
% [cylindrical front -> dispersion dominated regime]

function solnS2 =		solver_S2_ET(simul,fltr)

%	obtaining computational parameters
	nx =				simul.simil.nx;
	absmaxx =			simul.simil.absmaxx;
	power =				simul.simil.power;
	eta =				fltr.eta;
	gamm =				fltr.nondim.gamm;
	Pe =				fltr.nondim.Pe;
	derivwarn =			simul.derivwarn;
	leftBCtuned =		simul.simil.leftBCtuned;

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
	if (leftBCtuned == 1)
		B(1) =			-K*x(1);
	else
		B(1) =			1.0;
	end
	for ix = 2:nx-1
		acoeffx1 =		differentcoeff(x,ix,1,'cd',derivwarn);
		acoeffx2 =		differentcoeff(x,ix,2,'cd',derivwarn);
		for jx = -1:1
			row(icounter) =		ix;
			col(icounter) =		ix+jx;
% 			jaco(icounter) =	eta*acoeffx2(3+jx)+((barq*x(ix)+x(ix)^2)/3)*acoeffx1(3+jx);
			jaco(icounter) =	eta*acoeffx2(3+jx)+((barq*x(ix))/3)*acoeffx1(3+jx);
			icounter =			icounter+1;
		end
	end
	icounter =					icounter-1;
	G =							(sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B)';
 	G =							G/max(G);
	disp('solver_S2 Early Time solved');

%	obtaining front width integral
	IntNr =				trapz(x,x.^2.*(G.^2+K*G.*x));
	IntDr =				trapz(x,G.^2+K*G.*x);
	IntNr_alt =			trapz(x,x.^2.*(barq+x).*(G.^2+K*G.*x));
	IntDr_alt =			trapz(x,(barq+x).*(G.^2+K*G.*x));
%	containing output in output data structure
	solnS2.x =								x;
	solnS2.G =								G;
	solnS2.R =								G.^2+K*G.*x;
	[solnS2.Rmax, solnS2.locRmax] =			max(solnS2.R);
	solnS2.K =								K;
	solnS2.barq =							barq;
	solnS2.baralpha =						sqrt(IntNr/IntDr);
	solnS2.baralpha_alt =					sqrt(IntNr_alt/IntDr_alt);
	solnS2.baralphaapprox =					0.979197710209806*(1./sqrt(tildPe))*(K/sqrt(tildPe)).^(-1/3);
	solnS2.barbeta =						G(1+(nx-1)/2)^2;
	solnS2.barbetaapprox =					0.5454*(K/sqrt(tildPe))^(2/3);
	solnS2.barzeta_0 =						(2*pi*IntDr_alt);
	solnS2.barZ_0 =							solnS2.barzeta_0/(5/3);

end