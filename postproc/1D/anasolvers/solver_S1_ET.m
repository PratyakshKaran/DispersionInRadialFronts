% tests the front-attached frame solution for phi_A corresponding to the similarity solution 
% [cylindrical front -> no dispersion scenario & dispersion negligible regime]

function solnS1 =		solver_S1_ET(simul,fltr)

%	obtaining computational parameters
	nx =				simul.simil.nx;
	absmaxx =			simul.simil.absmaxx;
	power =				simul.simil.power;
	gamm =				fltr.nondim.gamm;
	Pe =				fltr.nondim.Pe;
	derivwarn =			simul.derivwarn;
	leftBCtuned =		simul.simil.leftBCtuned;
	barQthres =			simul.simil.barQthres;

	barQ =				Pe/2;
	if (barQ > barQthres)
		z1 =			sqrt(2*barQ+erfinv((1-gamm)/(1+gamm)));
		barq =			2*z1;
		if (gamm ~= 1)
			K =			((1+gamm)/sqrt(pi))*(exp(-(erfinv((1-gamm)/(1+gamm)))^2)*	...
						(z1/erfinv((1-gamm)/(1+gamm)))*(1-(barQ/(z1^2))));
		else
			K =			1.0;
		end
	else
		barq =			2*sqrt(gammaincinv(gamm/(gamm+1),barQ,'upper'));
		K =				(1+gamm)*(barq/2)*(1/gamma(barQ))*(((barq^2)/4)^(barQ-1)*exp(-((barq^2)/4)));
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
	if (leftBCtuned == 1)
		B(1) =		-K*x(1);
	else
		B(1) =		1.0;
	end
	for ix = 2:nx-1
		acoeffx1 =		differentcoeff(x,ix,1,'cd',derivwarn);
		acoeffx2 =		differentcoeff(x,ix,2,'cd',derivwarn);
		for jx = -1:1
			row(icounter) =		ix;
			col(icounter) =		ix+jx;
			jaco(icounter) =	(barq+x(ix))*acoeffx2(3+jx)+	...
								((1-Pe)+((Pe*(barq*x(ix)+x(ix)^2))/2))*acoeffx1(3+jx);
			icounter =			icounter+1;
		end
	end
	icounter =			icounter-1;
	G =					(sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B)';	
	G =					G/max(G);
	disp('solver_S1 Early Time solved');
%	obtaining front width integral
	IntNr =				trapz(x,x.^2.*(G.^2+K*G.*x));
	IntDr =				trapz(x,G.^2+K*G.*x);
	IntNr_alt =			trapz(x,x.^2.*(barq+x).*(G.^2+K*G.*x));
	IntDr_alt =			trapz(x,(barq+x).*(G.^2+K*G.*x));
%	containing output in output data structure
	solnS1.x =								x;
	solnS1.G =								G;
	solnS1.R =								G.^2+K*G.*x;
	[solnS1.Rmax, solnS1.locRmax] =			max(solnS1.R);
	solnS1.K =								K;
	solnS1.barq =							barq;
	solnS1.baralpha =						sqrt(IntNr/IntDr);
	solnS1.baralpha_alt =					sqrt(IntNr_alt/IntDr_alt);
	solnS1.baralphaapprox =					0.979197710209806*(1./sqrt(Pe))*(K/sqrt(Pe)).^(-1/3);
	solnS1.barbeta =						G(1+(nx-1)/2)^2;
	solnS1.barbetaapprox =					0.5454*(K/sqrt(Pe))^(2/3);
	solnS1.barzeta_0 =						(2*pi*IntDr_alt);
	solnS1.barzeta_0 =						(2*pi*barq*(1/Pe)*K);
	solnS1.barZ_0 =							solnS1.barzeta_0/2;
	
end