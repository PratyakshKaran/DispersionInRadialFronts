%-----------------------------------------------------------------------------------------------------------------------------------
% Generic Function to generate neigbour coefficients for numerical derivative approximation using Finite Difference Method (FDM)
% 
%	written by --
%	pratyaksh karan (pratyakshkaran@gmail.com)
%
% computes neighbour coefficients corresponding to FFM approxcimations of the derivative of an independent variables (x), 
% where the discetized grid for 'x' can be irregular
% inputs:
% x						array containing values of independent variable w.r.t. which the variable is to be computed
% ix					index representing location on the discretized x-axis at which the derivative is to be computed
% orddiff				order of differentiation -- 'x' represents first order and 'xx' represents second order
% dirdiff				direction of discretization -- 'fd' is 'forward-differenced', 'bd' is 'backward-differenced', 'cd' is
%						'centrally-differenced'; direction overridden at first grid point to 'fd' and at last grid point to 'bd'
% derivwarn				whether to derivative manual warnings occuring when computing derivative
% outputs:
% differentcoeff		five-element array containing neigbour coefficients, where the indices 1,2,3,4,5 respectively represents the
%						nodes ix-2,ix-1,ix,ix+1,ix+2
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
function differentcoeff = differentcoeff(x,ix,orddiff,dirdiff,derivwarn)

	nx = length(x);
	
	differentcoeff = zeros(1,5);

	if (orddiff == 1)
		if (strcmp(dirdiff,'bd') && ((ix == 1) || (ix == 2)))
			if derivwarn
				disp([	'Direction of differentiation is ', dirdiff, ' with the node index as ', num2str(ix),	...
						', thus overriding direction of differentiation to cd'				]);
			pause;
			end
			dirdiff = 'cd';
		end
		if (strcmp(dirdiff,'fd') && ((ix == nx) || (ix == nx-1)))
			if derivwarn
				disp([	'Direction of differentiation is ', dirdiff, ' with the node index as ', num2str(ix),	...
						' with grid length ', num2str(nx),														...
						', thus overriding direction of differentiation to cd'				]);
			end
			dirdiff = 'cd';
			pause;
		end
		if (strcmp(dirdiff,'cd') && (ix == 1))
			if derivwarn
				disp([	'Direction of differentiation is ', dirdiff, ' with the node index as ', num2str(ix),	...
						', thus overriding direction of differentiation to fd'				]);
			end
			dirdiff = 'fd';
			pause;
		end
		if (strcmp(dirdiff,'fd') && (ix == nx))
			if derivwarn
				disp([	'Direction of differentiation is ', dirdiff, ' with the node index as ', num2str(ix),	...
						' with grid length ', num2str(nx),														...
						', thus overriding direction of differentiation to bd'				]);
			end
			dirdiff = 'bd';
			pause;
		end
		if (strcmp(dirdiff,'cd'))
			dx =					x(ix+1)-x(ix);
			idx =					x(ix)-x(ix-1);
			differentcoeff(3+1) =		idx/(2.0*dx*idx);
			differentcoeff(3+0) =		(dx-idx)/(2.0*dx*idx);
			differentcoeff(3+-1) =		-dx/(2.0*dx*idx);
		elseif (strcmp(dirdiff,'fd'))
			dx =					x(ix+1)-x(ix);
			dxi =					x(ix+2)-x(ix+1);
			differentcoeff(3+2) =		-(dx^2.0)/(dxi*dx*(dx+dxi));
			differentcoeff(3+1) =		((dx+dxi)^2.0)/(dxi*dx*(dx+dxi));
			differentcoeff(3+0) =		-((dx+dxi)^2.0-dx^2.0)/(dxi*dx*(dx+dxi));
		elseif (strcmp(dirdiff,'bd'))
			idx =					x(ix)-x(ix-1);
			iidx =					x(ix-1)-x(ix-2);
			differentcoeff(3+0) =		((idx+iidx)^2.0-idx^2.0)/(iidx*idx*(idx+iidx));
			differentcoeff(3-1) =		-((idx+iidx)^2.0)/(iidx*idx*(idx+iidx));
			differentcoeff(3-2) =		(idx^2.0)/(iidx*idx*(idx+iidx));
		else
			disp([	'Direction of derivative given as ', dirdiff, ' which is invalid'	]);
			pause;
		end
	elseif (orddiff == 2)
		if (strcmp(dirdiff,'bd') || strcmp(dirdiff,'fd'))
			if derivwarn
				disp([	'Order of derivative is ', num2str(orddiff),	...
						', overriding direction of differentiation from ',dirdiff,' to cd'	]);
			end
			dirdiff = 'cd';
			pause;
		end
		if (strcmp(dirdiff,'cd'))
			if ((ix > 1) && (ix < nx))
				differentcoeff(3-1) =		2.0/((x(ix)-x(ix-1))*(x(ix+1)-x(ix-1)));
				differentcoeff(3+1) =		2.0/((x(ix+1)-x(ix))*(x(ix+1)-x(ix-1)));
				differentcoeff(3+0) =		-2.0/((x(ix+1)-x(ix))*(x(ix)-x(ix-1)));
			else
				disp([	'Second derivative cannot be computed at end node, and node inputted is ',num2str(ix)]);
				pause;
			end
		else
			disp([	'Direction of derivative given as ', dirdiff, ' which is invalid'	]);
			pause;
		end
	else
		disp([		'Order of derivative is ', num2str(orddiff), ', derivatives higher than order 2 not coded in'	]);
		pause;
	end

end
%-----------------------------------------------------------------------------------------------------------------------------------