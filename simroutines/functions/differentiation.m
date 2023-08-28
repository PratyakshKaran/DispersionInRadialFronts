%-----------------------------------------------------------------------------------------------------------------------------------
% Generic Function to compute numerical derivative approximation using Finite Difference Method (FDM)
% 
%	written by --
%	pratyaksh karan (pratyakshkaran@gmail.com)
%
% computes non-mixed derivative of a dependent variable (u) that is a function of an independent variables (x), where the discetized
% grid for 'x' can be irregular
% inputs:
% u						array containing values of dependent variable whose derivative is to be computed
% x						array containing values of independent variable w.r.t. which the variable is to be computed
% ix					index representing location on the discretized x-axis at which the derivative is to be computed
% orddiff				order of differentiation -- 'x' represents first order and 'xx' represents second order
% dirdiff				direction of discretization -- 'fd' is 'forward-differenced', 'bd' is 'backward-differenced', 'cd' is
%						'centrally-differenced'; direction overridden at first grid point to 'fd' and at last grid point to 'bd'
% derivwarn				whether to derivative manual warnings occuring when computing derivative
% outputs:
% differentiation		numerical values of the derivative
%-----------------------------------------------------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------------------
function differentiation = differentiation(u,x,ix,orddiff,dirdiff,derivwarn)

	nx =							length(x);
	differentcoeffin =				differentcoeff(x,ix,orddiff,dirdiff,derivwarn);
	
	differentiation =				0.0;
	for jx = -2:2
		iix = ix+jx;
		if ((iix >= 1) && (iix <= nx))
			differentiation =		differentiation+differentcoeffin(3+jx)*u(iix);
		end
	end
end
%-----------------------------------------------------------------------------------------------------------------------------------