% tests the front-attached frame solution for phi_A corresponding to the early-time similarity solution
% [spherical front -> dispersion dominated regime]

function solnS4 =		solver_S4(simul,fltr)

%	obtaining computational parameters
	nx =				simul.simil.nx;
	absmaxx =			simul.simil.absmaxx;
	power =				simul.simil.power;
	eta =				fltr.eta;
	gamm =				fltr.nondim.gamm;
	Pe =				fltr.nondim.Pe;
	derivwarn =			simul.derivwarn;
	
%	generating x-grid
	xl =				flip(-absmaxx*linspace(0,1,(nx-1)/2+1).^power);
	xr =				absmaxx*linspace(0,1,(nx-1)/2+1).^power;
	x =					[xl,xr(2:end)];
	
%	initiating equation solver matrices
	row =				zeros(3*nx-4,1);
	col =				zeros(3*nx-4,1);
	jaco =				zeros(3*nx-4,1);
	B =					zeros(nx,1);
	
%	obtaining front and mass of product properties
	Rf =				(gamm^2)/(1+gamm^2);
	Rmax =				gamm/4;
	baralpha_half =		sqrt(2)*(16*eta)^(1/4);
	baralpha =			sqrt(2/5);
	barq_adv =			(3)^(1/3);
	barq_minusadv =		2*((1-gamm)/(1+gamm))*eta^(1/4);

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
	B(1) =				1.0;
	icounterresum =		icounter;
	
%	populating loop-varying jacobian entries
	icounter =			icounterresum;
	for ix = 2:nx-1
		acoeffx1 =		differentcoeff(x,ix,1,'cd',derivwarn);
		acoeffx2 =		differentcoeff(x,ix,2,'cd',derivwarn);
		for jx = -1:1
			row(icounter) =	ix;
			col(icounter) =	ix+jx;
			if (eta == 0)
				jaco(icounter) =		...
				acoeffx2(3+jx)+((Pe*x(ix))/2)*acoeffx1(3+jx);
			else
				jaco(icounter) =		...
				acoeffx2(3+jx)+((x(ix)^3)/(4*eta))*acoeffx1(3+jx);
			end
			icounter =		icounter+1;
		end
	end
	icounter =			icounter-1;

%	numerical solution of G
	G =					sparse(row(1:icounter),col(1:icounter),jaco(1:icounter))\B;

%	containing output in output data structure
	solnS4.x =								x;
	solnS4.G =								G;
	solnS4.phiA =							G;
	solnS4.phiB =							(gamm/2)+gamm*((1/2)-solnS4.phiA);
	solnS4.R =								(solnS4.phiA.*solnS4.phiB);
	[solnS4.Rmax, solnS4.locRmax] =			max(solnS4.R);
	solnS4.theta =							solnS4.phiA-solnS4.phiB;
	solnS4.barq_adv =						barq_adv;
	solnS4.barq_minusadv =					barq_minusadv;
	solnS4.Rf =								Rf;
	solnS4.Rmax =							Rmax;
	solnS4.baralpha_half =					baralpha_half;
	solnS4.baralpha =						baralpha;
	if (eta ~= 0)
		solnS4.barzeta_0 =					3^(-1/3)*8*pi*eta^(1/4);
		solnS4.barzeta_1 =					((8*pi)/15)*eta^(1/4);
		solnS4.barZ_0 =						(3/5)*solnS4.barzeta_0;
		solnS4.barZ_1 =						(2/3)*solnS4.barzeta_1;
	else
		solnS4.barzeta_0 =					(3^(2/3)*4*sqrt(pi))/sqrt(Pe);
		solnS4.barzeta_1 =					(5/3)*sqrt(pi)*((4*sqrt(2))/(sqrt(Pe)^3));
		solnS4.barZ_0 =						(3/5)*solnS4.barzeta_0;
		solnS4.barZ_1 =						(1/2)*solnS4.barzeta_1;
	end

end