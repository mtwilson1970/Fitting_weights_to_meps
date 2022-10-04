function [m,c,R2]=fit_and_correlate(x,y,include_list,xlab,ylab,titlab,plotnum)
% Marcus Wilson 29 September 2021
% Do a correlation plot of variables x and y, find R2 value and extract
% line of best fit
%
% x - an array of x-values
% y - an array of y-values
% include_list is a list of which people to include in the analysis
% xlabs, ylabs are character strings for the x and y axis labels
% titlab gives the title label
% plotnum is the subplot number we use - up to (4)
%
% m is the gradient of the bestfit line
% c is the constant of the bestfit line
% R2 is the R-squared value

figure(51)


%which ones do we include
xtemp=x'; ytemp=y;  %copy the x and y, and transpose the x array at the same time
x=xtemp(include_list);
y=ytemp(include_list);


sigmas=x./x;   %an array of ones to use as the sigma value for our bestline fit

subplot(1,2,plotnum)
grid on; hold on;
plot(x,y,'kx');   %plot data
xlabel(xlab)
ylabel(ylab,'interpreter','tex')
box on;


%Find the bestfit line
[m, c, uncm, uncc, R2]=first_order_line_mtw(x,y,sigmas)

xdash=[min(x) max(x)];
ydash=xdash*m+c;
plot(xdash,ydash,'k--');   %put on the dashed line

title([titlab ' R^2=' num2str(R2)],'interpreter','tex');     %put on the title


end