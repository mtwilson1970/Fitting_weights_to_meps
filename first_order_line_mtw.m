% first_order_line_mtw.m

% Marcus Wilson
% 13 September 2021

%Fits a linear regression line and gives uncertainties

function [a, b, unca, uncb, R2]=first_order_line_mtw(x,y,sigmas)

%[a, b, unca, uncb]=first_order_line_mtw(x,y,sigmas)

%x - contains x values
%y - contains y values
%sigmas - contains sigmas
%
%output
%a = linear coeff
%b=constant
%unca = uncertainty in a
%uncb=uncertainty in b
% R2 is the R-squared residual value (no account of sigmas)

S=sum( sigmas.^(-2) ); 
Sx=sum( x.*sigmas.^(-2) );
Sy=sum( y.*sigmas.^(-2) );
Sxx=sum( x.*x.*sigmas.^(-2) );
Sxy=sum( x.*y.*sigmas.^(-2) );
Delta=S*Sxx - Sx^2;   

a=(S*Sxy-Sx*Sy)/Delta;
b=(Sxx*Sy - Sx*Sxy)/Delta;

unca=sqrt(   ( (1/(length(x)-2))*sum( ( (y-(a*x + b) ).^2 ).*sigmas.^(-2) ) / ( sum( ((x-mean(x)).^2).*sigmas.^(-2))) ) );
uncb=unca*sqrt(   (1/length(x))*(Sxx/S));

R2= 1 - sum((y - (a*x + b) ).^2)/sum(( y-mean(y)).^2);


end
