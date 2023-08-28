%-----------------------------------------------------------------------------------------------------------------------------------
% 1D Steady State Darcy Flow Solver corresponding for Planar/Cylindrical/Spherical Injection using Finite Difference Method (FDM)
% 
%	written by --
%	uddipta ghosh (uddipta.ghosh@iitgn.ac.in)
%	pratyaksh karan (pratyakshkaran@gmail.com)
%
% solves the piezometric head (H) and discharge rate (vel) for Darcy flow in a domain with a specified UNIFORM permeability field
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
function [H, vel, Kperm] = solver1D_darcy(geomdom,simul,fltr,grd)

%	stripping structures to variables
	x =								grd.x;
	nx =							simul.nx;
	Kperm =							fltr.Kperm;
	flowbc_left =					geomdom.flowbc.left;
	flowbc_right =					geomdom.flowbc.right;
	flowbc_Hleft0 =					fltr.flowbc.Hleft0;
	flowbc_Hright0 =				fltr.flowbc.Hright0;
	isradial =						geomdom.isradial;
	solveGPU = 						simul.solveGPU;
	derivwarn =						simul.derivwarn;
	confmap_space =					simul.confmap.space;
	
%	obtaining gridsizes
	nsparse =						nx*3;

%	initializing matrices
	resid =							zeros(nx,1);
	jaco =							zeros(nsparse,1);
	row =							zeros(nsparse,1);
	col =							zeros(nsparse,1);
	H =								zeros(nx,1);

%	solving flow field if required
	if (simul.givenvel == 1)
		disp 'CAUTION: Velocity Field to be Specified, overriding to Zero'
		H =							zeros(nx,1);
		vel.v =						zeros(nx,1);
	else
%		populating jacobian and residual
%		initializing sparse matrix element entry counter
		icounter =					1;
%		populating inner grid points
		for ix = 2:nx-1
%			obtaining derivative coefficients
			acoeff1 =				differentcoeff(x,ix,1,'cd',derivwarn);
			acoeff2 =				differentcoeff(x,ix,2,'cd',derivwarn);
			for jx = -1:1
				row(icounter) =		ix;
				col(icounter) =		ix+jx;
				if (isradial == 0)
					jaco(icounter)=	acoeff2(3+jx)+((1.0/Kperm(ix))*differentiation(Kperm,x,ix,1,'cd',derivwarn))*acoeff1(3+jx);
				else
					if ((confmap_space == 0) || (isradial == 0))
						jaco(icounter)=	...
									acoeff2(3+jx)+((1.0/Kperm(ix))*differentiation(Kperm,x,ix,1,'cd',derivwarn)+	...
									(double(isradial)/x(ix)))*acoeff1(3+jx);
					else
						jaco(icounter)=	...
									(1-x(ix))*x(ix)*acoeff2(3+jx)+	...
									((((1-x(ix))*x(ix))/Kperm(ix))*differentiation(Kperm,x,ix,1,'cd',derivwarn)+	...
									(double(isradial)-2*x(ix)))*acoeff1(3+jx);
					end
				end
				icounter =			icounter+1;
			end
		end
%		populating injection plane/line/point grid point
		ix = 1;
		if (flowbc_left == 1)
			row(icounter) =			ix;
			col(icounter) =			ix;
			jaco(icounter) =		1.0;
			icounter =				icounter+1;
		else
			acoeff1 =				differentcoeff(x,ix,1,'fd',derivwarn);
			for jx = 0:2
				row(icounter) =		ix;
				col(icounter) =		ix+jx;
				if ((confmap_space == 0) || (isradial == 0))
					jaco(icounter)=	acoeff1(3+jx);
				else
					jaco(icounter)=	((1-x(ix))^2)*acoeff1(3+jx);
				end
				icounter =			icounter+1;
			end
		end
%		populating residual
		resid(ix) =					flowbc_Hleft0;
%		populating far-end grid point
		ix = nx;
		if (flowbc_right == 1)
			row(icounter) =			ix;
			col(icounter) =			ix;
			jaco(icounter) =		1.0;
			icounter =				icounter+1;
		else
			acoeff1 =				differentcoeff(x,ix,1,'bd',derivwarn);
			for jx = -2:0
				row(icounter) =		ix;
				col(icounter) =		ix+jx;
				if ((confmap_space == 0) || (isradial == 0))
					jaco(icounter)=	acoeff1(3+jx);
				else
					jaco(icounter)=	((1-x(ix))^2)*acoeff1(3+jx);
				end
				icounter =			icounter+1;
			end
		end
%		populating residual
		resid(ix) =					flowbc_Hright0;
%		trimming icounter to get filled sparse matrix
		icounter =					icounter-1;

%		obtaining numerical solution
		jacosp =					sparse(row(1:icounter),col(1:icounter),jaco(1:icounter));
		if (solveGPU == 0)
			Hsoln =					jacosp\resid;
		else
			gjacosp =				gpuArray(jacosp);
			gresid =				gpuArray(resid);
			Hsoln =					gjacosp\gresid;
		end
	
%		updating pressure head matrix
		H(1:nx) =					Hsoln(1:nx);
	
%		obtaining velocities
		v =							zeros(nx,1);
		if ((confmap_space == 0) || (isradial == 0))
			for ix = 2:nx-1
				v(ix) =				-Kperm(ix)*differentiation(H,x,ix,1,'cd',derivwarn);
			end
			ix = 1;
			v(ix) =					-Kperm(ix)*differentiation(H,x,ix,1,'fd',derivwarn);
			ix = nx;
			v(ix) =					-Kperm(ix)*differentiation(H,x,ix,1,'bd',derivwarn);
		else
			for ix = 2:nx-1
				v(ix) =				-((1-x(ix))^2)*Kperm(ix)*differentiation(H,x,ix,1,'cd',derivwarn);
			end
			ix = 1;
			v(ix) =					-((1-x(ix))^2)*Kperm(ix)*differentiation(H,x,ix,1,'fd',derivwarn);
			ix = nx;
			v(ix) =					-((1-x(ix))^2)*Kperm(ix)*differentiation(H,x,ix,1,'bd',derivwarn);
		end
	
%		updating output velocity structure
		vel.v =		v;
	end
end
%-----------------------------------------------------------------------------------------------------------------------------------