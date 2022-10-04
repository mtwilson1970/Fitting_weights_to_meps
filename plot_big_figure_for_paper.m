%This goes with the fit_weights_to_goldsworthy_vallence.m and plots out the
%final figures


%11 January 2022
%plot out a final composite panel for the paper. Might still need a bit of
%editing after.
figure(200)
subplot(4,2,1)
%first plot is the IO curves.   We want the model fits as pink lines with no error bars
%to start
plot_IO_2(amps_in_rmt,output_curve_TBS_in_array,+1,'m-',1,'m-',0,'MEP (mV)');  grid on; hold on %The post-TBS model, as a pink line
plot_IO_2(amps_in_rmt,output_curve_base_in_array,-1,'m--',1,'m--',0,'MEP (mV)'); %The baseline model, as a pink dashed line

%Now the measured results
plot_IO_2(amps_in_rmt,iTBS_s_ave,+1,'k-',0,'k-',1,'MEP (mV)');   %post TBS measured
plot_IO_2(amps_in_rmt,base_s_ave,-1,'k--',0,'kx--',1,'MEP (mV)');   %baseline measured


%4 October 2022
%Find R2 values for part (a)

%baseline
R2base=1 - sum( (output_curve_base_in_array - base_s_ave).^2) / sum( ( base_s_ave - mean(base_s_ave) ).^2 );
R2post=1 - sum( (output_curve_TBS_in_array - iTBS_s_ave).^2) / sum( ( iTBS_s_ave - mean(iTBS_s_ave) ).^2 );

title(['(a), R^2_{pre}=' num2str(R2base) ',R^2_{post}=' num2str(R2post)],'interpreter','TeX')
box on



subplot(4,2,2)
%This one is PC1 expt vs PC1 model
plot(pc1_expt,pc1_model,'kx')
xlabel('PC1 experiment')
ylabel('PC1 model')
hold on; grid on;
%fit trendline
p=polyfit(pc1_expt,pc1_model,1);
xfits=[min(pc1_expt) max(pc1_expt)];   %the x-values for the line
yfits=xfits*p(1)+p(2);   %the yfits
plot(xfits,yfits,'k--')
yfitted_all=pc1_expt*p(1)+p(2);   %fit to all points
R2=1 - sum((pc1_model-yfitted_all).^2)/sum((pc1_model-mean(pc1_model)).^2);   %the R2 value
title(['(b), R-squared=' num2str(R2)])

subplot(4,2,3)
%This one the relative change post-TBS
to_plot=[1:18];  %the people to put on the plot. Sometimes we have problems here
%Plot the post/pre for each of the expt and model conditions
%After Nigel's email 17 August 2021
expt_post_over_pre = iTBS_s_ave(to_plot,:)./base_s_ave(to_plot,:);
model_post_over_pre = output_curve_TBS_in_array(to_plot,:)./output_curve_base_in_array(to_plot,:);
plot_IO_2(amps_in_rmt,expt_post_over_pre,-1,'kx-',1,'k-',1,'Relative change post-TBS')
plot_IO_2(amps_in_rmt,model_post_over_pre,1,'mo-',1,'m-',1,'Relative change post-TBS')
title('(c)')
box on


subplot(4,2,4)
%This one the change in wee-wie versus PC1 expt
delta=(wee_iTBS_array-wie_iTBS_array)' - (wee_base_array-wie_base_array)';  %after - before. Change in net excitation
%now plot against PC1
plot(pc1_expt(1:18),delta(1:18),'kx'); grid on; hold on
xlabel('PC1 experiment')
ylabel('change in w_{ee}')%-wie')
p=polyfit(pc1_expt(1:18),delta(1:18),1);
xfits=[min(pc1_expt(1:18)) max(pc1_expt(1:18))];   %the x-values for the line
yfits=xfits*p(1)+p(2);   %the yfits
plot(xfits,yfits,'k--')
yfitted_all=pc1_expt(1:18)*p(1)+p(2);   %fit to all points
R2=1 - sum((delta(1:18)-yfitted_all).^2)/sum((delta(1:18)-mean(delta(1:18))).^2);   %the R2 value
title(['(d), R-squared=' num2str(R2)])

subplot(4,2,5)
% %This the changes in wee and wie as trajectories
% for i=1:length(ee_move_list)  %length of people list
%     plot(ee_move_list(:,i),ie_move_list(:,i),'kx-')  %draw a connecting line for the time periods
%     xlabel('change in wee');
%     ylabel('change in wie');
%     grid on; hold on;
%     title('(e)')
%     pbaspect([1 1 1])
%     set(gca,'xlim',[-0.2 0.2])
%     set(gca,'ylim',[-0.2 0.2])
% end
pre=wee_base_array'; post=wee_iTBS_array';
tbl=table(pre,post)
p=parallelplot(tbl,'DataNormalization','none','Color','k');
grid on
title('(e) w_{ee}');

subplot(4,2,6);
histogram(wee_iTBS_array-wee_base_array,binedges)
xlabel('Movement in w_{ee}')
ylabel('Number of subjects')
title('(f)'); grid on

subplot(4,2,7)
 pre=(amp_for_this_person').^(-1); post=(amp_after_TBS').^(-1);
 tbl=table(pre,post)
 p=parallelplot(tbl,'DataNormalization','none','Color','k');
 grid on
 title('(g)  \mu');
 
subplot(4,2,8)
(amp_after_TBS - amp_for_this_person)'
[h,p]=ttest(amp_after_TBS,amp_for_this_person)
disp('')
histogram(amp_after_TBS.^(-1)-amp_for_this_person.^(-1),binedges)
xlabel('Movement in \mu')
ylabel('Number of subjects')
grid on; 
title('(h)')
%

figure(201)
%plot residuals as histograms
subplot(2,2,1)
resids_bin_edges=[-1.05:0.1:1.01]

resids_base_110=(base_s_ave(:,3) - output_curve_base_in_array(:,3));   
histogram(resids_base_110,resids_bin_edges);  grid on; hold on;
ylabel('Number of subjects')
xlabel('Experimental - Model')
[ha,pa]=adtest(resids_base_110) %anderson darling
[ht,pt]=ttest(resids_base_110) %ttest
title(['(a) baseline 110 p=' num2str(pa) ' n, p=' num2str(pt) ' m'])
set(gca,'ylim',[0 10])

subplot(2,2,2)
resids_base_180=(base_s_ave(:,10) - output_curve_base_in_array(:,10));   
histogram(resids_base_180,resids_bin_edges);  grid on; hold on;
ylabel('Number of subjects')
xlabel('Experimental - Model')
[ha,pa]=adtest(resids_base_180) %anderson darling
[ht,pt]=ttest(resids_base_180) %ttest
title(['(b) baseline 180 p=' num2str(pa) ' n, p=' num2str(pt) ' m'])
set(gca,'ylim',[0 10])



subplot(2,2,3)
resids_TBS_110=(iTBS_s_ave(:,3) - output_curve_TBS_in_array(:,3));   
histogram(resids_TBS_110,resids_bin_edges);  grid on; hold on;
ylabel('Number of subjects')
xlabel('Experimental - Model')
[ha,pa]=adtest(resids_TBS_110) %anderson darling
[ht,pt]=ttest(resids_TBS_110) %ttest
title(['(c) post-TBS 110, p=' num2str(pa) ' n, p=' num2str(pt) ' m'])
set(gca,'ylim',[0 10])


subplot(2,2,4)
resids_TBS_180=(iTBS_s_ave(:,10) - output_curve_TBS_in_array(:,10));   
histogram(resids_TBS_180,resids_bin_edges);  grid on; hold on;
ylabel('Number of subjects')
xlabel('Experimental - Model')
[ha,pa]=adtest(resids_TBS_180) %anderson darling
[ht,pt]=ttest(resids_TBS_180) %ttest
set(gca,'ylim',[0 10])

title(['(d) post-TBS 180 p=' num2str(pa) ' n, p=' num2str(pt) ' m'])

