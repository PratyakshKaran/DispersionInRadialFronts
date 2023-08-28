%-----------------------------------------------------------------------------------------------------------------------------------
% 1D Stationary/Steady Reactive Transport Solver for Spherical Injection using Finite Difference Method (FDM)
% 
%	written by --
%	uddipta ghosh (uddipta.ghosh@iitgn.ac.in)
%	pratyaksh karan (pratyakshkaran@gmail.com)
%
% solves the steady state reactive properties (including concentrations of species A and B) for the A+B->C reaction in the 
% pre-specified spherical injection flow field
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
function phi =					stationary1D_reactive(geomdom,simul,fltr,grd,flow)

%	stripping structures to variables
	nx =						simul.nx;
	Pe =						fltr.nondim.Pe;
	v =							flow.vel.v;
	x =							grd.x;
	transpbc_left =				geomdom.transpbc.left;
	transpbc_right =			geomdom.transpbc.right;
 	transpbc_phileft0 =			fltr.transpbc.phileft0;
 	transpbc_phiright0 =		fltr.transpbc.phiright0;
	STLtol =					simul.STLtol;
	errSTLmode =				simul.errSTLmode;
	STexplicit =				simul.STexplicit;
	STLrelax =					simul.STLrelax;
	derivwarn =					simul.derivwarn;
	isradial =					geomdom.isradial;
	solveGPU = 					simul.solveGPU;
	eta =						fltr.eta;
	nodiffusion =				simul.nodiffusion;

%	Pre-Allocating Field Variables
	phi =						zeros(nx,3);
	phiprev =					zeros(nx,3);
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
		phi(ix,3) =				0*fltr.init.phiA(ix);
	end

%	Creating Output Folder
%	mkdir('../outdata/1D');
%	mkdir(['../outdata/1D/',simul.casename]);

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
			C2i =					eta*Pe*v(ix)+double(1-nodiffusion);
			if (isradial == 0)
				C1i =				0.0;
			else
				C1i =				eta*Pe*differentiation(v,x,ix,1,'cd',derivwarn)+	...
									(eta*Pe*v(ix)+double(1-nodiffusion))/(x(ix)/double(isradial));
			end
			for jx = -(1+double(ix>=3)):1
				jacoC(ix,3+(jx),ispecies) =	v(ix)*acoeffx1bd(3+jx)-(1/Pe)*(C2i*acoeffx2(3+jx)+C1i*acoeffx1(3+jx));
			end
		end
%		Implementing Boundary Condition for the Injection Plane/Line/Point
		ix = 1;
		if (transpbc_left(ispecies) == 1)
			jacoC(ix,3+(+0),ispecies) =		1.0;
		else
			acoeffx1 =						differentcoeff(x,ix,1,'fd',derivwarn);
			for jx = 0:2
				jacoC(ix,3+jx,ispecies) =	acoeffx1(3+jx);
			end
		end
		residC(ix,ispecies) =				double(ispecies~=3)*transpbc_phileft0(ispecies)+	...
											double(ispecies==3)*(transpbc_phileft0(1)-transpbc_phileft0(2));
%		Implementing Boundary Condition for the Right End
		ix = nx;
		if (transpbc_right(ispecies) == 1)
			jacoC(ix,3+(+0),ispecies) =		1.0;
		else
			acoeffx1 =						differentcoeff(x,ix,1,'bd',derivwarn);
			for jx = -2:0
				jacoC(ix,3+jx,ispecies) =	acoeffx1(3+jx);
			end
		end
		residC(ix,ispecies) =				double(ispecies~=3)*transpbc_phiright0(ispecies)+	...
											double(ispecies==3)*(transpbc_phiright0(1)-transpbc_phiright0(2));
	end

%	Time-Marching
	for it = 1:1
%		Obtaining Time Step Size
		coeffcurrent =				0.0;
		coeffprev =					0.0;
		coeffprevprev =				0.0;
%		Saving Previous Time Step Solution
		phiprevprev =				phiprev;
		phiprev =					phi;
		if ((STexplicit == 0) || ((STexplicit == 1) && (it <= tstepprevs+1)))
%			Resetting Jacobian and Residual
			jacoT =						jacoT*0.0;
			residT =					residT*0.0;
%			Time-Step Dependent Jacobian and Residular Matrices
			for ispecies = 1:3
%				Updating Inner Region Jacobian
				for ix = 2:nx-1
					jacoT(ix,3+(+0),ispecies) =	...
										coeffcurrent;
				end
				for ix = 2:nx-1
					residT(ix,ispecies)=-coeffprev*phiprev(ix,ispecies)-coeffprevprev*phiprevprev(ix,ispecies);
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
					end
				else
%					Source Term Expressions
					for ix = 1:nx
						Sc(ix,1) =			-phi(ix,1)*phi(ix,2);
						Sc(ix,2) =			-phi(ix,1)*phi(ix,2);
					end
				end
%				Iterating over Reactants and Product
				for ispecies = 2*istage-1:istage+1
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
					errSTL =					sqrt((1.0/nx)*	...
												sum(sum((phi(:,2*istage-1:istage+1)-philast(:,2*istage-1:istage+1)).^2)));
				elseif (errSTLmode == 2)
					errSTL =					sqrt(sum(sum((phi(:,2*istage-1:istage+1)-philast(:,2*istage-1:istage+1)).^2)));
				elseif (errSTLmode == 3)
					errSTL =					max(max(abs(phi(:,2*istage-1:istage+1)-philast(:,2*istage-1:istage+1))));
				end
				prog = [	' iteration number ', num2str(iter), ' with STL error ', num2str(errSTL)	];
				disp(prog);
			end
		end
	end
end
%-----------------------------------------------------------------------------------------------------------------------------------