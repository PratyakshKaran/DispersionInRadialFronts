%-----------------------------------------------------------------------------------------------------------------------------------
% 1D Reactive Transport Warped Time Solver using Finite Difference Method (FDM)
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
function solnWT =			solver_WT(geomdom,simul,fltr,grd)

%	stripping structures to variables
%	WT specific parameters
	indtgrid =					simul.WT.indtgrid;
	nt =						simul.WT.nt;
	nx =						simul.WT.nx;
	absmaxx =					simul.WT.absmaxx;
	power =						simul.WT.power;
	m_WT =						simul.WT.m_WT;
	r0 =						simul.WT.r0;
	delt0 =						simul.WT.delt0;
	k_WT =						geomdom.isradial;
	num_rf =					simul.WT.num_rf;
	num_wf =					simul.WT.num_wf;
	if (num_wf == 1)
		wf_in =					simul.WT.wf;
	end
	if (num_rf == 1)
		rf_in =					simul.WT.rf;
	end
	if ((num_wf == 1) || (num_rf == 1))
		tf_in =					simul.WT.tf;
	end
%	common parameters
	if (indtgrid == 0)
		nt =					simul.nt;
	end
	Da =						fltr.nondim.Da;
	Pe =						fltr.nondim.Pe;
	Q =							fltr.flowbc.Q;
% 	xact =						grd.x;
	t =							grd.t;
	transpbc_left =				geomdom.transpbc.left;
	transpbc_right =			geomdom.transpbc.right;
 	transpbc_phileft0 =			fltr.transpbc.phileft0;
 	transpbc_phiright0 =		fltr.transpbc.phiright0;
	STLtol =					simul.STLtol;
	errSTLmode =				simul.errSTLmode;
	STexplicit =				simul.STexplicit;
	STLrelax =					simul.STLrelax;
	derivwarn =					simul.derivwarn;
	solveGPU = 					simul.solveGPU;
	tstepsave =					unique(simul.tstepsave);
	tstepprevs =				simul.tstepprevs;
	eta =						fltr.eta;

	if (m_WT == 1)
		C_WT =					(eta*Pe*Q)/(delt0*(r0^k_WT));
	else
		C_WT =					1.0;
	end
	A =							C_WT/(Pe*Q*(delt0^(m_WT))*(r0^(m_WT*k_WT))*(m_WT*k_WT+k_WT+1));

%	obtaining time grid, space grid, advective front and strip thickness
	rf =						(r0^(k_WT+1)+(k_WT+1)*Q*t).^(1/(k_WT+1));
	xl =						flip(-absmaxx*linspace(0,1,(nx-1)/2+1).^power);
	xr =						absmaxx*linspace(0,1,(nx-1)/2+1).^power;
	x =							[xl,xr(2:end)];
	if (indtgrid == 1)
 		tstart =				A*(rf(1).^(m_WT*k_WT+k_WT+1)-r0^(m_WT*k_WT+k_WT+1));
 		tend =					A*(rf(end).^(m_WT*k_WT+k_WT+1)-r0^(m_WT*k_WT+k_WT+1));
		t =						linspace(tstart,tend,nt);
		tact =					(((t/A)+r0^(m_WT*k_WT+k_WT+1)).^((k_WT+1)/(m_WT*k_WT+k_WT+1))-r0^(k_WT+1))/((k_WT+1)*Q);
	else
		tact =					t;
 		t =						A*(rf.^(m_WT*k_WT+k_WT+1)-r0^(m_WT*k_WT+k_WT+1));
	end
	rf =						(r0^(k_WT+1)+(k_WT+1)*Q*tact).^(1/(k_WT+1));

%	Checking if Time Grid is Regular
	regt = linspace(t(1),t(end),length(t));
	if (sum(regt~=t)==0)
		istreg = 1;
	else
		istreg = 0;
	end
	dt =						t(2)-t(1);
	if ((istreg == 0) && (STexplicit == 1))
		disp('Explicit Source Term Solution is not useful for irregular time grid, resetting to Implicit Source Term Solution');
		STexplicit =			0;
	end
	if (tstepprevs > 1)
		acoefft =				differentcoeff(t,3,1,'bd',derivwarn);
	end

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
		phi(ix,3) =				fltr.init.phiC(ix);
	end
	for ix = 1:nx
		phiprev(ix,1) =			fltr.init.phiA(ix);
		phiprev(ix,2) =			fltr.init.phiB(ix);
		phiprev(ix,3) =			fltr.init.phiC(ix);
	end

%	opening progress log file
	proglog =			fopen([simul.outfoldrname,'/progress.log'],'a+t');

%	Creating Output Folder
	simul.outfoldrnameorig =		simul.outfoldrname;
	simul.outfoldrname =			[simul.outfoldrname,'/WTSoln'];
	mkdir(simul.outfoldrname);

%	Opening Dependent Variables Output Files
	r_WT =					rf(1)+x*((delt0*(r0^k_WT))/(rf(1)^k_WT));
	rout =					fopen([simul.outfoldrname,'/r_1.dat'],'wt');
	phiAout =				fopen([simul.outfoldrname,'/phiA_1.dat'],'wt');
	phiBout =				fopen([simul.outfoldrname,'/phiB_1.dat'],'wt');
	phiCout =				fopen([simul.outfoldrname,'/phiC_1.dat'],'wt');
%	Writing Intial Conditions in Output Files
	for ix = 1:nx
		fprintf(rout,'%d \t ',r_WT(ix));
		fprintf(phiAout,'%d \t ',phi(ix,1));
		fprintf(phiBout,'%d \t ',phi(ix,2));
		fprintf(phiCout,'%d \t ',phi(ix,3));
	end
	fprintf(rout,'\n');
	fprintf(phiAout,'\n');
	fprintf(phiBout,'\n');
	fprintf(phiCout,'\n');
%	Closing Files
	fclose(rout);
	fclose(phiAout);
	fclose(phiBout);
	fclose(phiCout);

%	Constant Jacobian and Residular Matrices
	if (num_wf == 0)
		for ispecies = 1:3
%			Implementing Boundary Condition for the Injection Plane/Line/Point
			ix = 1;
			if (transpbc_left(ispecies) == 1)
				jacoC(ix,3+(+0),ispecies) =	1.0;
			else
				acoeffx1 =					differentcoeff(x,ix,1,'fd',derivwarn);
				for jx = 0:2
					jacoC(ix,3+jx,ispecies)=acoeffx1(3+jx);
				end
			end
			residC(ix,ispecies) =			transpbc_phileft0(ispecies);
%			Implementing Boundary Condition for the Right End
			ix = nx;
			if (transpbc_right(ispecies) == 1)
				jacoC(ix,3+(+0),ispecies) =	1.0;
			else
				acoeffx1 =					differentcoeff(x,ix,1,'bd',derivwarn);
				for jx = -2:0
					jacoC(ix,3+jx,ispecies)=acoeffx1(3+jx);
				end
			end
			residC(ix,ispecies) =			transpbc_phiright0(ispecies);
		end
	end
	for ispecies = 1:3
%		Implementing Boundary Condition for the Injection Plane/Line/Point
		ix = 1;
		residC(ix,ispecies) =				transpbc_phileft0(ispecies);
%		Implementing Boundary Condition for the Right End
		ix = nx;
		residC(ix,ispecies) =				transpbc_phiright0(ispecies);
	end

%	Obtaining current time front width and front location if enforced
	if (num_rf == 1)
		rf =		pchip(tf_in,rf_in,tact);
	end
	if (num_wf == 1)
		wf =		pchip(tf_in,wf_in,tact);
	end

%	Time-Marching
	tstepsaveind = 1;
	for it = 2:length(t)
%		resetting absmax and x grid if width front is to be taken from numerical solution
		if (num_wf == 1)
			absmaxx =					10*(0.5*wf(it))*((rf(it)^k_WT)/(delt0*(r0^k_WT)));
			xl =						flip(-absmaxx*linspace(0,1,(nx-1)/2+1).^power);
			xr =						absmaxx*linspace(0,1,(nx-1)/2+1).^power;
			x =							[xl,xr(2:end)];
		end
%		Obtaining Time Step Size
		if ((it == 2) || (tstepprevs == 1))
			if (istreg == 1)
				coeffcurrent =			1.0/dt;
				coeffprev =				-1.0/dt;
				coeffprevprev =			0.0;
			else
				coeffcurrent =			1.0/(t(it)-t(it-1));
				coeffprev =				-1.0/(t(it)-t(it-1));
				coeffprevprev =			0.0;
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
			if (num_wf == 1)
				for ispecies = 1:3
%					Implementing Boundary Condition for the Injection Plane/Line/Point
					ix = 1;
					if (transpbc_left(ispecies) == 1)
						jacoT(ix,3+(+0),ispecies) =		1.0;
					else
						acoeffx1 =						differentcoeff(x,ix,1,'fd',derivwarn);
						for jx = 0:2
							jacoT(ix,3+jx,ispecies) =	acoeffx1(3+jx);
						end
					end
%					Implementing Boundary Condition for the Right End
					ix = nx;
					if (transpbc_right(ispecies) == 1)
						jacoT(ix,3+(+0),ispecies) =		1.0;
					else
						acoeffx1 =						differentcoeff(x,ix,1,'bd',derivwarn);
						for jx = -2:0
							jacoT(ix,3+jx,ispecies) =	acoeffx1(3+jx);
						end
					end
				end
			end
			for ispecies = 1:3
%				Updating Inner Region Jacobian
				for ix = 2:nx-1
%					coefficients for derivatives
					acoeffx2 =			differentcoeff(x,ix,2,'cd',derivwarn);
					acoeffx1 =			differentcoeff(x,ix,1,'cd',derivwarn);
					for jx = -1:1
						jacoT(ix,3+jx,ispecies) =	...
										coeffcurrent*double(jx==0)-	...
										(1/((m_WT*k_WT+k_WT+1)*Pe*Q*A*(rf(it)^(m_WT*k_WT))))*	...
										((((eta*Pe*Q)/(rf(it)^k_WT))+1)*(((rf(it)^k_WT)/(delt0*(r0^k_WT)))^2)*acoeffx2(3+jx)+	...
										(k_WT/rf(it))*((rf(it)^k_WT)/(delt0*(r0^k_WT)))*acoeffx1(3+jx));
					end
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
				philast =				phi;
				if (STexplicit == 0)
%					Source Term Expressions
					for ix = 1:nx
						Sp(ix,1) =			-(Da/((m_WT*k_WT+k_WT+1)*Pe*Q*A*rf(it)^(m_WT*k_WT)))*phi(ix,2);
						Sp(ix,2) =			-(Da/((m_WT*k_WT+k_WT+1)*Pe*Q*A*rf(it)^(m_WT*k_WT)))*phi(ix,1);
						Sc(ix,3) =			(Da/((m_WT*k_WT+k_WT+1)*Pe*Q*A*rf(it)^(m_WT*k_WT)))*phi(ix,1)*phi(ix,2);
					end
				else
%					Source Term Expressions
					for ix = 1:nx
						Sc(ix,1) =			-(Da/((m_WT*k_WT+k_WT+1)*Pe*Q*A*rf(it)^(m_WT*k_WT)))*phi(ix,1)*phi(ix,2);
						Sc(ix,2) =			-(Da/((m_WT*k_WT+k_WT+1)*Pe*Q*A*rf(it)^(m_WT*k_WT)))*phi(ix,1)*phi(ix,2);
						Sc(ix,3) =			(Da/((m_WT*k_WT+k_WT+1)*Pe*Q*A*rf(it)^(m_WT*k_WT)))*phi(ix,1)*phi(ix,2);
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
%		Writing Output to File
		if (((tact(it) <= grd.t(tstepsave(min(tstepsaveind+1,length(tstepsave))))) && ...
				(tact(min(it+1,length(tact))) > grd.t(tstepsave(min(tstepsaveind+1,length(tstepsave))))))	||	...
				(it == length(tact)))
%			obtaining real x-grid
			r_WT =						rf(it)+x*((delt0*(r0^k_WT))/(rf(it)^k_WT));
			tstepsaveind =				tstepsaveind+1;
%			Opening Dependent Variables Output Files
			rout =						fopen([simul.outfoldrname,'/r_',num2str(tstepsave(tstepsaveind)),'.dat'],'wt');
			phiAout =					fopen([simul.outfoldrname,'/phiA_',num2str(tstepsave(tstepsaveind)),'.dat'],'wt');
			phiBout =					fopen([simul.outfoldrname,'/phiB_',num2str(tstepsave(tstepsaveind)),'.dat'],'wt');
			phiCout =					fopen([simul.outfoldrname,'/phiC_',num2str(tstepsave(tstepsaveind)	),'.dat'],'wt');
%			Writing Concentration Output to File
			for ix = 1:nx
				fprintf(rout,'%d \t ',r_WT(ix));
				fprintf(phiAout,'%d \t ',phi(ix,1));
				fprintf(phiBout,'%d \t ',phi(ix,2));
				fprintf(phiCout,'%d \t ',phi(ix,3));
			end
			fprintf(rout,'\n');
			fprintf(phiAout,'\n');
			fprintf(phiBout,'\n');
			fprintf(phiCout,'\n');
%			Closing Files
			fclose(rout);
			fclose(phiAout);
			fclose(phiBout);
			fclose(phiCout);	
%			displaying and logging progress of time marching
			prog = ['Written out time step ',num2str(it),' out of ',num2str(length(t))];
			disp(prog); fprintf(proglog,[prog,'\n']);
		end
%		displaying and logging progress of time marching
		prog = ['Solved time step ',num2str(it),' out of ',num2str(length(t))];
		disp(prog); fprintf(proglog,[prog,'\n']);
	end

	solnWT.phi =				phi;
	fclose(proglog);

	solnWT.m_WT =				m_WT;
	solnWT.k_WT =				k_WT;
	solnWT.C_WT =				C_WT;
	solnWT.A =					A;
	solnWT.rf =					rf;
	solnWT.wf =					rf;
 	solnWT.t =					t;
 	solnWT.tact =				tact;
 	solnWT.x =					x;
	solnWT.outfoldrname =		simul.outfoldrname;
	solnWT.nx =					nx;
	solnWT.absmaxx =			absmaxx;
	solnWT.power =				power;
	solnWT.m_WT =				m_WT;
	solnWT.k_WT =				k_WT;
	solnWT.r0 =					r0;
	solnWT.delt0 =				delt0;
	
end