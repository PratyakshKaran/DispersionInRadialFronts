%-----------------------------------------------------------------------------------------------------------------------------------
% 1D Transient Reactive Transport Solver for Planar/Cylindrical/Spherical Injection using Finite Difference Method (FDM)
% 
%	written by --
%	uddipta ghosh (uddipta.ghosh@iitgn.ac.in)
%	pratyaksh karan (pratyakshkaran@gmail.com)
%
% solves the temporal evolution of reactive properties (including concentrations of species A and B) for A+B->C reaction in the 
% pre-specified planar/cylindrical/spherical injection flow field
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
function solnreactive =			solver1D_reactive(geomdom,simul,fltr,grd,flow)

%	stripping structures to variables
	nx =						simul.nx;
	nt =						simul.nt;
	Pe =						fltr.nondim.Pe;
	v =							flow.vel.v;
	x =							grd.x;
	t =							grd.t;
	transpbc_left =				geomdom.transpbc.left;
	transpbc_right =			geomdom.transpbc.right;
 	transpbc_phileft0 =			fltr.transpbc.phileft0;
 	transpbc_phiright0 =		fltr.transpbc.phiright0;
	transpbc_inject_singular =	fltr.transpbc.inject.singular;
	STLtol =					simul.STLtol;
	errSTLmode =				simul.errSTLmode;
	STexplicit =				simul.STexplicit;
	STLrelax =					simul.STLrelax;
	derivwarn =					simul.derivwarn;
	isradial =					geomdom.isradial;
	solveGPU = 					simul.solveGPU;
	compMassC =					simul.compMassC;
	compfront =					simul.compfront;
	tstepprevs =				simul.tstepprevs;
	eta =						fltr.eta;
	geomdom_transpinit_rsm =	geomdom.transpinit.rsm;
	frcdispers =				simul.frcdispers;
	nodiffusion =				simul.nodiffusion;
	if ((frcdispers ~= 0) && (isradial == 0))
		velfrcdispers =			fltr.transpbc.frcdispers;
	else
		velfrcdispers =			0.0;
	end
	confmap_space =					simul.confmap.space;
	confmap_time =					simul.confmap.time;
	iterateC =						simul.iterateC;

%	Checking if Time Grid is Regular
	regt = linspace(t(1),t(end),length(t));
	if (sum(regt~=t)==0)
		istreg = 1;
	else
		istreg = 0;
	end
	dt =						t(2)-t(1);
	if (((istreg == 0) || (confmap_time == 1)) && (STexplicit == 1))
		disp('Explicit Source Term Solution is not useful for irregular/conformal mapped time grids, resetting to Implicit');
		STexplicit =			0;
	end
	if (tstepprevs > 1)
		acoefft =				differentcoeff(t,3,1,'bd',derivwarn);
	end

%	Setting Injection Point as zero if specified
	if (transpbc_inject_singular == 1)
		grd.x(1) =				0.0;
	end

%	Pre-Allocating Field Variables
	phi =						zeros(nx,3);
	phiprev =					zeros(nx,3);
	if (compMassC == 1)
		solnreactive.MassC =	zeros(nt,1);
		solnreactive.barR =		zeros(nt,1);
		solnreactive.MassC_4mR=	zeros(nt,1);
		solnreactive.MassC_4mR_twk =	...
								zeros(nt,1);
		solnreactive.dMassCdt =	zeros(nt,1);
	end
	if (compfront == 1)
		solnreactive.xf =		zeros(nt,1);
		solnreactive.Rf =		zeros(nt,1);
		solnreactive.wf =		zeros(nt,1);
		solnreactive.xforig =	zeros(nt,1);
		solnreactive.Rforig =	zeros(nt,1);
		solnreactive.wforig =	zeros(nt,1);
		solnreactive.R =		zeros(nx,1);
	end
%	Pre-Allocating Simulation Variables
	nsparse =					2*(nx-2)+6;
	row =						zeros(nsparse,1);
	col =						zeros(nsparse,1);
	jaco =						zeros(nsparse,1);
	resid =						zeros(nx,1);
	jacoC =						zeros(nx,5,3);
	residC =					zeros(nx,3);
	jacoT =						zeros(nx,5,3);
	residT =					zeros(nx,3);
	Sp =						zeros(nx,3);
	Sc =						zeros(nx,3);
	if (STexplicit == 1)
		jacospinv =				zeros(nx,nx,3);
	end
%	loading initial conditions
	for ix = 1:nx
		phi(ix,1) =				fltr.init.phiA(ix);
		phi(ix,2) =				fltr.init.phiB(ix);
		phi(ix,3) =				fltr.init.phiC(ix);
	end
	if (geomdom_transpinit_rsm == 1)
		for ix = 1:nx
			phiprev(ix,1) =		fltr.init.phiA(ix);
			phiprev(ix,2) =		fltr.init.phiB(ix);
			phiprev(ix,3) =		fltr.init.phiC(ix);
		end
	end

%	opening progress log file
	proglog =			fopen([simul.outfoldrname,'/progress.log'],'a+t');

%	Creating Output Folder
% 	mkdir('../outdata/1D');
% 	mkdir(['../outdata/1D/',simul.casename]);

	if (geomdom_transpinit_rsm == 0)
%		Opening Dependent Variables Output Files
		phiAout =				fopen([simul.outfoldrname,'/phiA_1.dat'],'wt');
		phiBout =				fopen([simul.outfoldrname,'/phiB_1.dat'],'wt');
		phiCout =				fopen([simul.outfoldrname,'/phiC_1.dat'],'wt');
%		Writing Intial Conditions in Output Files
		for ix = 1:nx
			fprintf(phiAout,'%d \t ',phi(ix,1));
			fprintf(phiBout,'%d \t ',phi(ix,2));
			fprintf(phiCout,'%d \t ',phi(ix,3));
		end
		fprintf(phiAout,'\n');
		fprintf(phiBout,'\n');
		fprintf(phiCout,'\n');
%		Closing Files
		fclose(phiAout);
		fclose(phiBout);
		fclose(phiCout);
	end

%	Constant Jacobian and Residular Matrices
	for ispecies = 1:3
%		Updating Inner Region Jacobian
		for ix = 2:nx-1
			if (ix >= 3)
				acoeffx1bd =		differentcoeff(x,ix,1,'bd',derivwarn);
			else
				acoeffx1bd =		zeros(1,5);
				acoeffx1bd(3+(+0))=	1.0/(x(ix)-x(ix-1));
				acoeffx1bd(3+(-1))=	-1.0/(x(ix)-x(ix-1));				
			end
			acoeffx1 =				differentcoeff(x,ix,1,'cd',derivwarn);
			acoeffx2 =				differentcoeff(x,ix,2,'cd',derivwarn);
			if ((confmap_space == 0) || (isradial == 0))
				C2i =				eta*Pe*(v(ix)+velfrcdispers)+double(1-nodiffusion);
			else
				C2i =				(eta*Pe*(v(ix)+velfrcdispers)+double(1-nodiffusion))*((1-x(ix))^4);
			end
			if (isradial == 0)
				C1i =				0.0;
			else
				if (confmap_space == 0)
					C1i =			eta*Pe*differentiation(v,x,ix,1,'cd',derivwarn)+	...
									(eta*Pe*(v(ix)+velfrcdispers)+double(1-nodiffusion))/(x(ix)/double(isradial));
				else
					C1i =			eta*Pe*((1-x(ix))^4)*differentiation(v,x,ix,1,'cd',derivwarn)+	...
									(eta*Pe*(v(ix)+velfrcdispers)+double(1-nodiffusion))*	...
									((double(isradial)-2*x(ix))/x(ix))*((1-x(ix))^3);
				end
			end
			for jx = -(1+double(ix>=3)):1
				if ((confmap_space == 0) || (isradial == 0))
					jacoC(ix,3+(jx),ispecies) =	v(ix)*acoeffx1bd(3+jx)-(1/Pe)*(C2i*acoeffx2(3+jx)+C1i*acoeffx1(3+jx));
				else
					jacoC(ix,3+(jx),ispecies) =	((1-x(ix))^2)*v(ix)*acoeffx1bd(3+jx)-(1/Pe)*(C2i*acoeffx2(3+jx)+C1i*acoeffx1(3+jx));
				end
			end
		end
%		Implementing Boundary Condition for the Injection Plane/Line/Point
		ix = 1;
		if (transpbc_left(ispecies) == 1)
			jacoC(ix,3+(+0),ispecies) =		1.0;
		else
			acoeffx1 =						differentcoeff(x,ix,1,'fd',derivwarn);
			for jx = 0:2
				if ((confmap_space == 0) || (isradial == 0))
					jacoC(ix,3+jx,ispecies)=acoeffx1(3+jx);
				else
					jacoC(ix,3+jx,ispecies)=((1-x(ix))^2)*acoeffx1(3+jx);
				end
			end
		end
		residC(ix,ispecies) =				transpbc_phileft0(ispecies);
%		Implementing Boundary Condition for the Right End
		ix = nx;
		if (transpbc_right(ispecies) == 1)
			jacoC(ix,3+(+0),ispecies) =		1.0;
		else
			acoeffx1 =						differentcoeff(x,ix,1,'bd',derivwarn);
			for jx = -2:0
				if ((confmap_space == 0) || (isradial == 0))
					jacoC(ix,3+jx,ispecies)=acoeffx1(3+jx);
				else
					jacoC(ix,3+jx,ispecies)=((1-x(ix))^2)*acoeffx1(3+jx);
				end
			end
		end
		residC(ix,ispecies) =				transpbc_phiright0(ispecies);
	end

%	Time-Marching
	if (geomdom_transpinit_rsm == 1)
		itstart =		geomdom.transpinit.rsmit;
	else
		itstart =		2;
	end
	for it = itstart:simul.nt
%		Obtaining Time Step Size
		if ((it == 2) || (tstepprevs == 1))
			if (istreg == 1)
				coeffcurrent =				1.0/dt;
				coeffprev =					-1.0/dt;
				coeffprevprev =				0.0;
			else
				coeffcurrent =				1.0/(t(it)-t(it-1));
				coeffprev =					-1.0/(t(it)-t(it-1));
				coeffprevprev =				0.0;
			end
		else
			if (istreg == 0)
				acoefft =				differentcoeff(t,it,1,'bd',derivwarn);
			end
			coeffcurrent =				acoefft(3+0);
			coeffprev =					acoefft(3-1);
			coeffprevprev =				acoefft(3-2);
		end
%		Saving Previous Time Step Solution
		phiprevprev =					phiprev;
		phiprev =						phi;
		if ((STexplicit == 0) || ((STexplicit == 1) && (it <= tstepprevs+1)))
%			Resetting Jacobian and Residual
			jacoT =						jacoT*0.0;
			residT =					residT*0.0;
%			Time-Step Dependent Jacobian and Residular Matrices
			for ispecies = 1:3
%				Updating Inner Region Jacobian
				for ix = 2:nx-1
					if (confmap_time == 0)
						jacoT(ix,3+(+0),ispecies) =	...
										coeffcurrent;
					else
						jacoT(ix,3+(+0),ispecies) =	...
										((1-t(it))^2)*coeffcurrent;
					end
				end
				if (confmap_time == 0)
					for ix = 2:nx-1
						residT(ix,ispecies) =	...
										-coeffprev*phiprev(ix,ispecies)-coeffprevprev*phiprevprev(ix,ispecies);
					end
				else
					for ix = 2:nx-1
						residT(ix,ispecies) =	...
										((1-t(it))^2)*(-coeffprev*phiprev(ix,ispecies)-coeffprevprev*phiprevprev(ix,ispecies));
					end
				end
			end
			jacoT =						jacoT+jacoC;
			residT =					residT+residC;
		end
%		Source Term Linearized Iterations
		iter =							0;
		for istage = 1:2
			STLtol_in =					STLtol;
			errSTL =					STLtol_in*10;
			while (errSTL > STLtol_in)
%				setting error to zero when iterations are not needed
				if ((istage==2) || (STexplicit == 1))
					STLtol_in =				1.0E100;
				end
%				Iteration Counter
				iter =						iter+1;
%				Saving Iteration Solution
				philast =					phi;
				if (STexplicit == 0)
%					Source Term Expressions
					for ix = 1:nx
						Sp(ix,1) =			-phi(ix,2);
						Sp(ix,2) =			-phi(ix,1);
						Sc(ix,3) =			phi(ix,1)*phi(ix,2);
					end
				else
%					Source Term Expressions
					for ix = 1:nx
						Sc(ix,1) =			-phi(ix,1)*phi(ix,2);
						Sc(ix,2) =			-phi(ix,1)*phi(ix,2);
						Sc(ix,3) =			phi(ix,1)*phi(ix,2);
					end
				end
%				Iterating over Reactants and Product
				if (iterateC == 1)
					if (istage == 1)
						leftispecies =		1;
						rightispecies =		3;
					else
						leftispecies =		1;
						rightispecies =		0;
					end
				else
					leftispecies =			2*istage-1;
					rightispecies =			istage+1;
				end
				for ispecies = leftispecies:rightispecies
%					Setting Jacobian
					if ((STexplicit == 0) || ((STexplicit == 1) && (it <= tstepprevs+1)))
%						Resetting Jacobian
						jaco =					jaco*0.0;
						icounter =				1;
%						Updating Inner Region
						for ix = 2:nx-1
							for jx = -(1+double(ix>=3)):1
								row(icounter) =	ix;
								col(icounter) =	ix+jx;
								jaco(icounter)=	jacoT(ix,3+jx,ispecies)-double(jx==0)*Sp(ix,ispecies);
								icounter =		icounter+1;
							end
						end
%						Implementing Boundary Conditions
						for ix = [1,nx]
							for jx = 0-2*double(ix==nx):2-2*double(ix==nx)
								row(icounter) =	ix;
								col(icounter) =	ix+jx;
								jaco(icounter)=	jacoT(ix,3+jx,ispecies);
								icounter =		icounter+1;
							end
						end
%						Trimming icounter to get Filled Sparse Matrix
						icounter =				icounter-1;
%						Updating Solution
						jacosp =				sparse(row(1:icounter),col(1:icounter),jaco(1:icounter));
						if (STexplicit == 1)
							jacospinv_in =		inv(jacosp);
							for ix = 1:nx
								for jx = 1:nx
									jacospinv(ix,jx,ispecies) =		...
												jacospinv_in(ix,jx);
								end
							end
							if ((STexplicit == 1) && (it > tstepprevs+1))
								clear jaco;
								clear jacosp;
								clear jacospinv_in;
							end
						end
					end
%					Setting Residual
%					Resetting Residual
					resid =						resid*0.0;
%					Updating Inner Region
					for ix = 2:nx-1
						resid(ix) =				residT(ix,ispecies)+Sc(ix,ispecies);
					end
%					Implementing Boundary Conditions
					for ix = [1,nx]
						resid(ix) =				residT(ix,ispecies);
					end
%					Updating Solution
					if (STexplicit == 0)
						if (solveGPU == 0)
							phi(1:nx,ispecies)=	(1.0-STLrelax*STLrelax^(1-istage))*phi(1:nx,ispecies)+	...
												STLrelax*STLrelax^(1-istage)*(jacosp\resid);
						else
							gjacosp =			gpuArray(jacosp);
							gresid =			gpuArray(resid);
							phi(1:nx,ispecies)=	(1.0-STLrelax*STLrelax^(1-istage))*phi(1:nx,ispecies)+	...
												STLrelax*STLrelax^(1-istage)*(gjacosp\gresid);
						end
					else
						phi(1:nx,ispecies) =	jacospinv*resid;
					end
				end
%				Computing STL Iteration Error
				if (errSTLmode == 1)
					errSTL =				sqrt((1.0/nx)*sum(sum((phi(:,2*istage-1:istage+1)-philast(:,2*istage-1:istage+1)).^2)));
				elseif (errSTLmode == 2)
					errSTL =				sqrt(sum(sum((phi(:,2*istage-1:istage+1)-philast(:,2*istage-1:istage+1)).^2)));
				elseif (errSTLmode == 3)
					errSTL =				max(max(abs(phi(:,2*istage-1:istage+1)-philast(:,2*istage-1:istage+1))));
				end
				prog = [	'time step ', num2str(it), ' of ', num2str(nt),					...
							' iteration number ', num2str(iter), ' with STL error ', num2str(errSTL)	];
				fprintf(proglog,prog);
				fprintf(proglog,'\n');
				disp(prog);
			end
		end
		if (sum(it == simul.tstepsave) ~= 0)
%			Opening Dependent Variables Output Files
			phiAout =					fopen([simul.outfoldrname,'/phiA_',num2str(it),'.dat'],'wt');
			phiBout =					fopen([simul.outfoldrname,'/phiB_',num2str(it),'.dat'],'wt');
			phiCout =					fopen([simul.outfoldrname,'/phiC_',num2str(it),'.dat'],'wt');
%			Writing Concentration Output to File
			for ix = 1:nx
				fprintf(phiAout,'%d \t ',phi(ix,1));
				fprintf(phiBout,'%d \t ',phi(ix,2));
				fprintf(phiCout,'%d \t ',phi(ix,3));
			end
			fprintf(phiAout,'\n');
			fprintf(phiBout,'\n');
			fprintf(phiCout,'\n');
%			Closing Files
			fclose(phiAout);
			fclose(phiBout);
			fclose(phiCout);	
%			displaying and logging progress of time marching
			prog = ['Written out time step ',num2str(it),' out of ',num2str(nt)];
			disp(prog); fprintf(proglog,[prog,'\n']);
		end
%		displaying and logging progress of time marching
		prog = ['Solved time step ',num2str(it),' out of ',num2str(nt)];
		disp(prog); fprintf(proglog,[prog,'\n']);
%		obtaining mass of C
		if (compMassC == 1)
			solnreactive.MassC(it) = ...
								((double(isradial)*2*pi)^double(isradial~=0))*trapz(x',phi(:,3).*(x'.^isradial));
			solnreactive.barR(it) =	...
								((double(isradial)*2*pi)^double(isradial~=0))*trapz(x',phi(:,1).*phi(:,2).*(x'.^isradial));	
			if (it ~= 1)
				solnreactive.MassC_4mR(it) = ...
									0.5*(solnreactive.barR(it)+solnreactive.barR(it-1))*(t(it)-t(it-1));
				solnreactive.MassC_4mR_twk(it) = ...
									solnreactive.MassC_4mR(it)+2*pi*phi(1,3)*(t(it)-t(it-1));
			end
		end
%		obtaining front properties
		if (compfront == 1)
			for ix = 1:simul.nx
				solnreactive.R(ix) =		phi(ix,1)*phi(ix,2);
			end
			[~,inodfrorig] =				min((phi(:,1)-phi(:,2)).^2);
			solnreactive.xforig(it) =		x(inodfrorig);
			solnreactive.Rforig(it)	=		solnreactive.R(inodfrorig);
			IntNr =							0.0;
			IntDr =							0.0;
			for ix = 1:simul.nx-1
				IntNr =						IntNr+0.5*	...
											((x(ix+1)-x(inodfrorig))^2*solnreactive.R(ix+1)+	...
											(x(ix)-x(inodfrorig))^2*solnreactive.R(ix))*	...
											(x(ix+1)-x(ix));
				IntDr =						IntDr+0.5*(solnreactive.R(ix+1)+solnreactive.R(ix))*(x(ix+1)-x(ix));
			end
			solnreactive.wforig(it) =		sqrt(IntNr/IntDr);
			[~,maxRloc] =					max(solnreactive.R);
			solnreactive.xf(it) =			grd.x(maxRloc);
			solnreactive.Rf(it) =			solnreactive.R(maxRloc);
			leftswitch =					1;
			rightswitch =					1;
			ixleft =						1;
			ixright = nx;
			for ix = 1:nx
				if ((solnreactive.R(ix) >= 0.5*solnreactive.Rf(it)) && (leftswitch == 1))
					ixleft =				ix;
					leftswitch =			0;
				end
				if ((solnreactive.R(ix) < 0.5*solnreactive.Rf(it)) && (leftswitch == 0) && (rightswitch == 1))
					ixright =				ix-1;
					rightswitch =			0;
				end
			end
			solnreactive.wf(it) =			x(ixright)-x(ixleft);
		end
	end
%	obtaining rate of production of C
	if (compMassC == 1)
		solnreactive.dMassCdt(2) =			(solnreactive.MassC(2)-solnreactive.MassC(1))/(t(2)-t(1));
		for it = 3:nt
			solnreactive.dMassCdt(it) =		differentiation(solnreactive.MassC,t,it,1,'bd',derivwarn);
		end
	end
	solnreactive.phi =						phi;
	fclose(proglog);
end
%-----------------------------------------------------------------------------------------------------------------------------------