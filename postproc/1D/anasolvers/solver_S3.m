% solves the front-attached frame solution for phi_A corresponding to the stationary solution
% [spherical front]

function solnS3 =		solver_S3(simul,fltr)

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

%	computing reuqired parameters
	if (eta < 0.1*Pe)
		K =				((log(1+gamm))^2)/Pe;
		r_fS =			Pe/log(1+gamm);
		tildPe =		1;
 	elseif (eta > 10*Pe)
 		K =				2/(pi*sqrt(eta*Pe));
 		r_fS =			sqrt(eta*Pe);
 		tildPe =		2;
	else
		K =				(gamm+exp((pi/2)*sqrt(Pe/eta)))/	...
						(eta*(exp((pi/2)*sqrt(Pe/eta))-1)*	...
						(1+(tan(sqrt(eta/Pe)*log(1-(1/(1+gamm))*(1-exp((pi/2)*sqrt(Pe/eta))))))^2));
		r_fS =			sqrt(eta*Pe)*tan(sqrt(eta/Pe)*log(1-(1/(1+gamm))*(1-exp((pi/2)*sqrt(Pe/eta)))));
		tildPe =		1/((((eta*Pe)/(r_fS^2))+1));
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
% 	if (leftBCtuned == 1)
% 		G =				-x(1)+x(1)*((x-x(1))/(x(nx)-x(1)));
% 	else
% 		G =				1-((x-x(1))/(x(nx)-x(1)));
% 	end
	G =					0*x;

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
									tildPe*(x(ix)+(2-double(strcmp(rootfindmeth,'STL')==1))*G(ix))*double(jx==0);
				icounter =			icounter+1;
			end
			if (strcmp(rootfindmeth,'STL') == 0)
				B(ix) =				differentiation(G,x,ix,2,'cd',derivwarn)-tildPe*G(ix)*(x(ix)+G(ix));
			end
		end
		icounter =			icounter-1;
		if (strcmp(rootfindmeth,'STL') == 0)
			if (leftBCtuned == 1)
				B(1) =		G(1)+x(1);
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
		disp(['solver_S3 Root finding error is ',num2str(errrootfind)]);
	end

%	obtaining front width integral
	IntNr =				trapz(x,(K^(-1/3)*x+r_fS).^2.*x.^2.*(G.^2+K*G.*x));
	IntDr =				trapz(x,(K^(-1/3)*x+r_fS).^2.*(G.^2+K*G.*x));
	IntNr_alt =			trapz(x,x.^2.*(G.^2+K*G.*x));
	IntDr_alt =			trapz(x,G.^2+K*G.*x);

%	containing output in output data structure
	solnS3.x_orig =				x;
	solnS3.x =					K^(-1/3)*x;
	solnS3.r =					solnS3.x+r_fS;
	solnS3.G =					G/max(G);
	solnS3.phiAS =				(K^(2/3))*G;
	solnS3.thetaS =				-K*solnS3.x;
	solnS3.phiBS =				(K^(2/3))*(G+x);
	solnS3.RS =					(K^(4/3))*(G.^2-G.*x);
	solnS3.K =					K;
	solnS3.r_fS =				r_fS;
	solnS3.baralphaa =			sqrt(IntNr/IntDr);
	solnS3.baralphaa_alt =		sqrt(IntNr_alt/IntDr_alt);
	solnS3.baralphaa_approx =	0.979197710209806*(1./sqrt(tildPe))*(K/sqrt(tildPe)).^(-1/3);
	solnS3.w_fS =				solnS3.baralphaa*K^(-1/3);
	solnS3.w_fS_alt =			solnS3.baralphaa_alt*K^(-1/3);
	solnS3.w_fS_approx =		solnS3.baralphaa_approx*K^(-1/3);
	solnS3.RFS =				solnS3.RS(floor((nx-1)/2+1));
	solnS3.locRmax =			solnS3.r(floor((nx-1)/2+1));
	solnS3.barzeta_0 =			(((eta*Pe)/(r_fS^2))+1)*K^(4/3);
	solnS3.barZ_0 =				solnS3.barzeta_0;
	
end