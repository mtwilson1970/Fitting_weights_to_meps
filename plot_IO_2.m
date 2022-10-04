%
% Script for plotting curves. Several people, do mean and s.e.m. plot
% Marcus Wilson 16 August 2021.
% Edited 11 January 2022 to make more versatile
% rmts are the resting motor thresholds as percent
% meps are the MEP values as mV
% deltax is the horizontal displacement for s.e.m. bars so they don't
% overlap
% type_of_curve is the character array for type of curve, e.g. k--
% width_of_curve is the width of the curve
% type_of_bars is the character array for the type of the sem bars, e.g. kx--
% width_of_bars is their width
% ytext_stuff is the text for the yaxis
function [status]=plot_IO_2(rmts,meps,deltax,type_of_curve,width_of_curve,type_of_bars,width_of_bars,ytext_stuff);
%First the curve

%Only plot it if width is bigger than zero. Some intro calculations now
 grid on; hold on;
 xlabel('TMS amplitude (percent RMT)')
 ylabel(ytext_stuff)
 mean_curve=mean(meps,1);

 if (width_of_curve==0)
    %create some axes, but don't plot
else
    plot(rmts,mean_curve,type_of_curve,'linewidth',width_of_curve);   
end
%
%
%The the sem, if we plot it
if (width_of_bars==0)    %so if there's zero width we don't plot the bars
    %Do nothing
else
    sem=std(meps,1)/sqrt(size(meps,1));   %find the standard error. std/sqrt(n)
    for i=1:size(meps,2); %go over the rmt values
        plot([rmts(i)+deltax rmts(i)+deltax],[mean_curve(i)-sem(i) mean_curve(i)+sem(i)],type_of_bars,'linewidth',width_of_bars);
    end
    set(gca,'xlim',[min(rmts)-5   max(rmts)+5])
end


status=1;

end

