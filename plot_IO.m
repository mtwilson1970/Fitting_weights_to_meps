%
% Script for plotting curves. Several people, do mean and s.e.m. plot
% Marcus Wilson 16 August 2021
% rmts are the resting motor thresholds as percent
% meps are the MEP values as mV
% deltax is the horizontal displacement for s.e.m. bars so they don't
% overlap
% type_of_curve is the character array for type of curve, e.g. k--
% width_of_curve is the width
%
function [status]=plot_IO(rmts,meps,deltax,type_of_curve,width_of_curve);
%First the curve
mean_curve=mean(meps,1);
plot(rmts,mean_curve,type_of_curve,'linewidth',width_of_curve);
grid on; hold on;
xlabel('TMS amplitude (percent RMT)')
ylabel('MEP (mV)')
%
%Extract the colour for the error bars
%type_colour=type_of_curve(1)
type_colour=type_of_curve;   %the full one
%
%The the sem
sem=std(meps,1)/sqrt(size(meps,1));   %find the standard error. std/sqrt(n)
for i=1:size(meps,2); %go over the rmt values
    plot([rmts(i)+deltax rmts(i)+deltax],[mean_curve(i)-sem(i) mean_curve(i)+sem(i)],type_colour,'linewidth',width_of_curve);
end
set(gca,'xlim',[min(rmts)-5   max(rmts)+5])

status=1;

end

