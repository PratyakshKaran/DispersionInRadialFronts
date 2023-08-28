%-----------------------------------------------------------------------------------------------------------------------------------
% Wrapper [Variant - Testing] for
% differentiation
% [Non-Function Program]
%
%	written by --
%	pratyaksh karan (pratyakshkaran@gmail.com)
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	initiating cleaned matlab sessions
clc;
clear;
close all;
fclose all;
addpath('../functions');
%-----------------------------------------------------------------------------------------------------------------------------------

A =				rand(1,5);
B =				rand(1,5);
C =				rand(1,5);
nx =			101;
x =				linspace(0.1,1.0,nx).^(1/2);
y =				0.0*x;
dy_ana =		0.0*x;
d2y_ana =		0.0*x;
dy_num_CD =		0.0*x;
dy_num_FD =		0.0*x;
dy_num_BD =		0.0*x;
d2y_num =		0.0*x;
for i = 1:5
	y =			y+A(i)*x.^i;
	dy_ana =	dy_ana+i*A(i)*x.^(i-1);
	d2y_ana =	d2y_ana+i*(i-1)*A(i)*x.^(i-2);
	y =			y+B(i)*exp(-i*x);
	dy_ana =	dy_ana-i*B(i)*exp(-i*x);
	d2y_ana =	d2y_ana+i^2*B(i)*exp(-i*x);
	y =			y+C(i)*cos(2*pi*(i/5)*x);
	dy_ana =	dy_ana-2*pi*(i/5)*C(i)*sin(2*pi*(i/5)*x);
	d2y_ana =	d2y_ana-4*pi^2*(i/5)^2*C(i)*cos(2*pi*(i/5)*x);
end
for j = 3:length(x)-2
	dy_num_CD(j) =		differentiation(y,x,j,1,'cd',true);
	dy_num_FD(j) =		differentiation(y,x,j,1,'fd',true);
	dy_num_BD(j) =		differentiation(y,x,j,1,'bd',true);
	d2y_num(j) =		differentiation(y,x,j,2,'cd',true);
end

figure();
axes();
plot(x,y,'r-','linewidth',2);
hold on;
plot(x,dy_ana,'g-','linewidth',2);
plot(x(3:end-2),dy_num_CD(3:end-2),'gs','linewidth',0.5);
plot(x(3:end-2),dy_num_FD(3:end-2),'go','linewidth',0.5);
plot(x(3:end-2),dy_num_BD(3:end-2),'gp','linewidth',0.5);
plot(x,d2y_ana,'b-','linewidth',2);
plot(x(3:end-2),d2y_num(3:end-2),'bs','linewidth',0.5);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
legend({'$y$',	...
		'$\displaystyle \frac{{\rm d} y}{{\rm d} x}$ (ana)',		...
		'$\displaystyle \frac{{\rm d} y}{{\rm d} x}$ (num-CD)',	...
		'$\displaystyle \frac{{\rm d} y}{{\rm d} x}$ (num-FD)',	...
		'$\displaystyle \frac{{\rm d} y}{{\rm d} x}$ (num-BD)',	...
		'$\displaystyle \frac{{\rm d}^2 y}{{\rm d} x^2}$ (ana)',		...
		'$\displaystyle \frac{{\rm d}^2 y}{{\rm d} x^2}$ (num-CD)'},		...
		'interpreter','latex');
%-----------------------------------------------------------------------------------------------------------------------------------