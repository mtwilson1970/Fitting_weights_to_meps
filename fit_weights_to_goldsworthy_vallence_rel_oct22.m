% 
%
% Marcus Wilson 17 July 2020.  Fit weights to Goldsworthy and Vallence data
%
%
%  10 Feb 2022. Fit the post-TBS curve to the relative values. See p 7734.
%  Hence the '_rel' in the name of the code. Calls
%  fit_weights_to_meps_rel.m
% 
%   4 October 2022. Tidy up initial bits to select parameters. Now calls
%   fit_weights_to_meps_rel_oct22.m
% 
close all; clear

%%%%%MAIN INPUT OPTIONS HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pick_goldsworthy=true;  
    %if true we pick goldsworthy iTBSdata, if false we pick the vallence cTBS data

select_a_session=false;   session_number=1;
    %Applies to Vallence cTBS data only. If true we only use data from a
    %single session (defined by session_number); if false we average over all sessions.
    %Standard is false.

    
S_select_para = 0;  %Standard = 0.     (Need to ensure do_S_approach=true or it defaults to option 1.)
    %This specifies the 'S' approach. Options are:
          %0 - We fit 'S' to the pre- and post- data independently (i.e. pre and post have different values of \mu)
          %1 - We don't fit S to either pre- or post- data - it is just it standard value (i.e. \mu=1 for both pre and post) 
          %2 - We fit 'S' to the pre- and use the same value for the post (so \mu_post = \mu_pre).
          
          
fix_ee=0;  fix_ie=2; 
   %fix_ee and fix ie options. 
       %If 0 we fix both pre- and -post independently
       %if 1 we don't fit either to pre- or post - just use their standard values of 1.92e-4 V s
       %if 2 we fit the value to the baseline and use the same value for the post-TBS
   %standard is fix_ee=0 (fit to both pre- and post-) and fix_ie=2 (fit to pre- and then use same value for -post)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%OTHER INPUT OPTIONS (generally leave as they are)%%%%%%%

do_we_scale=false; 
   %do we scale the data to make it all have the same average 180% baseline response (true), 
   %or do we choose to use the model to find the best fit amplitude or RMT for
   %each person (20 July 2020, false)? Standard = false

select_which_time_period=false;
   %If true we choose which of the three time period after the stimulation to
   %analyze.
   %If false we analyze all of them
time_period_min=3; time_period_max=3; %If we select a time period above, then this is the one that we use (1 to 3)

lsq_limit=9999999; %what is the cutoff for counting a person?   (%10 might be okay).  If really big (9999999) we accept all of them

do_S_approach = true; %Needs to be true to activate the 'S' approach
%if true (and do_we_scale is false, the standard) we stretch or squash the input axis - An 'S' model. 
%Standard = true  See p7107

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Revert back to older variable for selecting the approach
amp_or_squash=~do_S_approach;     %take the reverse of the S selection parameter for this one

%option check
if (  (S_select_para==0) || (S_select_para==1) || (S_select_para==2))
    disp(['S_select_para is ' num2str(S_select_para)])
else
    disp('S_select_para must be 0, 1 or 2')
    stop
end
if (  (fix_ee==0) || (fix_ee==1) || (fix_ee==2))
    disp(['fix_ee is ' num2str(fix_ee)])
else
    disp('fix_ee must be 0, 1 or 2')
    stop
end
if (  (fix_ie==0) || (fix_ie==1) || (fix_ie==2))
    disp(['fix_ie is ' num2str(fix_ie)])
else
    disp('fix_ie must be 0, 1 or 2')
    stop
end




if (pick_goldsworthy)
    name_of_data='Goldsworthy iTBS';
else
    name_of_data='Vallence cTBS sess 1';
end


%20 July 2020. We need to get back a fitted amplitude. Use amp_to_use=0,
%unless we choose not to, in which case it is 1;
if (do_we_scale)
    amp_to_use=1;
else
    amp_to_use=0;  %this means the fitting will find an amplitude, which it will report
end


for time_period=time_period_min:time_period_max

%read in data

if (~pick_goldsworthy)
 M=csvread('vallence_cTBS_3_sessions_10aug20.csv');   %all three sessions. 10 August. Correct by using the monophasic RMT not biphasic
 %create just one session if we need it
 if (select_a_session)
     %pick one of the vallence sessions
     N=M(:,(session_number-1)*50+2:(session_number-1)*50+51); 
     nn=[M(:,1)  N];   
     M=nn;   %create a new N file
 end
 
 if (select_a_session)   %just one session
     goldsworthy=2+0.1*session_number; %2.1 = session 1, 2.2 = session 2, 2.3 = session 3
 else
 goldsworthy=3;  %it's the vallence stuff with all 3 sessions
 end
% 
 %stdunc_base=[0.0033 0.0162 0.0307 0.0750 0.1142 0.1407 0.1524 0.1643 0.1626 0.1672];  %st unc in the baseline mean for Vallence, with mean over all people, trials and sessions and 
 stunc_base=[0.10 0.10 0.07 0.05 0.05 0.04 0.03 0.04 0.02 0.02]; %the approx relative expt uncerts
 %stunc_base=[1 1 1 1 1 1 1 1 1 1];   %a test case
 stunc_per_person = stunc_base*sqrt(size(M,1))   %See p7050.

else
    
% Goldsworthy here
 M = csvread('goldsworthy_iTBS_10aug20.csv');
 %stunc_base=[0.0031 0.0267 0.1012 0.1787 0.2468 0.2798 0.3029 0.2935 0.2928 0.3031];  %st unc in baseline mean averaged over all people and trials (2 of them). Just one session.
 stunc_base=[0.25 0.15 0.35 0.1 0.05 0.05 0.05 0.05 0.05 0.05];
  %stunc_base=[1 1 0.01 1  1 1 1 1 1 1 1 1];   %a test
 stunc_per_person = stunc_base*sqrt(size(M,1))   %See p7050.
 goldsworthy=1;   %it's the goldsworthy stuff
% %Also need to adjust rmt in the _meps code

end    %end the if for picking vallence or goldsworhty


%partition data
%structure for each row is RMT, baseline1 (90-180%), baseline 2, after iTBS
%1, 2 and then 3.
%For the vallence_cTBS_3_sessions data there are three sessions recorded.
if (length(M) > 51)
    
    %its the vallence 3 session data. Average over the three sessions.
   rmt_list(:) = M(:,1);
baseline1(:,:) = ( M(:,2:11) + M(:,52:61) + M(:,102:111) )/3;
baseline2(:,:) = ( M(:,12:21) + M(:,62:71) + M(:,112:121) )/3;

iTBS1(:,:) = ( M(:,22:31) + M(:,72:81) + M(:,122:131) )/3;
iTBS2(:,:) = ( M(:,32:41) + M(:,82:91) + M(:,132:141) )/3; 
iTBS3(:,:) = ( M(:,42:51) + M(:,92:101) + M(:,142:151) )/3;

    

else

    %Just read in the one-session case
    
    rmt_list(:) = M(:,1);

baseline1(:,:) = M(:,2:11);
baseline2(:,:) = M(:,12:21);

iTBS1(:,:) = M(:,22:31);
iTBS2(:,:) = M(:,32:41);
iTBS3(:,:) = M(:,42:51);



end


if (do_we_scale)  %do we try scaling data to start, or fit an amplitude (20 July 2020)
    
    %Find the normalization factor (p7041)
ave_people_base1_180=mean(baseline1(:,10));  %the mean of the 180% MEPs for baseline1. Scalar
ave_people_base2_180=mean(baseline2(:,10));  %the mean of the 180% MEPs for baseline1
ave_people_and_bases_180 = 0.5*(ave_people_base1_180 + ave_people_base2_180);  %the average 180% MEP responses over all people and baselines

ave_over_two_baselines_180=0.5*(baseline1(:,10) + baseline2(:,10));  %the average 180% MEP response for the baselines for all people. (Vector).

scaling(:)=ave_people_and_bases_180 * ave_over_two_baselines_180(:).^(-1);  %a vector. A different scaling for each person, but just one scaling per person

    
 %Now scale the data
for i=1:size(M,1)  %go over all rows (people)
    base1_s(i,:) = baseline1(i,:)*scaling(i);
    base2_s(i,:) = baseline2(i,:)*scaling(i);
    iTBS1_s(i,:) = iTBS1(i,:)*scaling(i);
    iTBS2_s(i,:) = iTBS2(i,:)*scaling(i);
    iTBS3_s(i,:) = iTBS3(i,:)*scaling(i);
end

else
    for i=1:size(M,1);  %go over all rows (people)
    base1_s(i,:) = baseline1(i,:);
    base2_s(i,:) = baseline2(i,:);
    iTBS1_s(i,:) = iTBS1(i,:);
    iTBS2_s(i,:) = iTBS2(i,:);
    iTBS3_s(i,:) = iTBS3(i,:);
    end
end




%Now find the average before and after points
base_s_ave = 0.5*(base1_s + base2_s); 

if (select_which_time_period)
    if (time_period==1)
        iTBS_s_ave=iTBS1;  %early result
    else if (time_period==2)
            iTBS_s_ave=iTBS2;   %mid result
        else
            iTBS_s_ave=iTBS3;   %later result
        end
    end
    
else
    %average over all
    iTBS_s_ave = (iTBS1 + iTBS2 + iTBS3)/3;
end

%Now have just a small set of data.  Find weights for each.

wee_b=[]; wie_b=[]; wee_i=[]; wie_i=[];

output_curve_base_in_array=[];
output_curve_TBS_in_array=[];

amp_for_this_person=NaN*ones(1,size(M,1)); %preallocate size
amp_after_TBS=NaN*ones(1,size(M,1));  

% %%%%for plots of expt results.  1 October 2021
% rmts=[90:10:180];
% meps=base_s_ave;
% figure
% hold on; box on;
% [status]=plot_IO(rmts,base_s_ave,+1,'k--',1);
% [status]=plot_IO(rmts,iTBS_s_ave,-1,'k-',1);
% stop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


base_null = ones(1,10);   % an array of ones that we will divide by first time through.
%Second call to fit_weights_to_meps_rel_oct22 will take the actual modeled baseline data.
mep_null=base_null;  %second call it will take the actual values of expt baseline data
%
%for i=18:18
for i=1:size(M,1);   %loop over all people
    person=i;   %which person
%for i=0:-1;  %a test line to skip this bit and go to average



disp(' ')
disp(['Person ' num2str(i)])

disp('baseline')
[wee_base, wie_base, lsq_base, fittedamp, RMT_actually_used,output_curve_base_pp]=fit_weights_to_meps_rel_oct22(base_s_ave(i,:),rmt_list(i),amp_to_use,stunc_per_person,amp_or_squash,goldsworthy,person,0.9,0.9,base_null,mep_null,fix_ee,fix_ie,S_select_para);   %the weights at end will be overwritten
output_curve_base_pp=max(output_curve_base_pp,0.001*ones(1,10)) ; %make the minimum value 0.01 for a baseline mep



wee_base_array(i) = wee_base;  wie_base_array(i)=wie_base;  lsq_base_array(i)=lsq_base;
amp_for_this_person(i) = fittedamp;   %load up the amplitude to use for this person.

output_curve_base_in_array=[output_curve_base_in_array; output_curve_base_pp]; %the baseline curve from the fitiing
%The fitting for the TBS

disp('TBS')
[wee_iTBS, wie_iTBS, lsq_iTBS, fittedamp, RMT_actually_used,output_curve_TBS_pp]=fit_weights_to_meps_rel_oct22(iTBS_s_ave(i,:),rmt_list(i),amp_for_this_person(i),stunc_per_person,amp_or_squash,goldsworthy,person,wee_base,wie_base,output_curve_base_pp,base_s_ave(i,:),fix_ee,fix_ie,S_select_para);  %20 July 2020. Fit with this amplitude.   10 Jan 2022. Bring in baseline wee and wie into the fitting code. 10feb22 bring in the baseline curve too


%This line is for fitting the scale factor to the post TBS, rather than copying from the baseline [wee_iTBS, wie_iTBS, lsq_iTBS, fittedamp, RMT_actually_used]=fit_weights_to_meps(iTBS_s_ave(i,:),rmt_list(i),amp_to_use,stunc_per_person,amp_or_squash,goldsworthy,person);  %20 July 2020. Fit with this amplitude
amp_after_TBS(i)=fittedamp;   %record the amplitude after the iTBS curve fitting.

wee_iTBS_array(i) = wee_iTBS;  wie_iTBS_array(i)=wie_iTBS;  lsq_iTBS_array(i)=lsq_iTBS;

output_curve_TBS_in_array=[output_curve_TBS_in_array; output_curve_TBS_pp];  %the output curve from the fitting

%check if lsq are less than chosen limit
%if ( (lsq_base<lsq_limit) & (lsq_iTBS<lsq_limit) )   %Check both
if (lsq_base<lsq_limit)  %just the baseline. So same number accepted each time
    wee_b=[wee_b  wee_base]; %add the next one
    wie_b=[wie_b wie_base];
    wee_i=[wee_i  wee_iTBS];
    wie_i=[wie_i   wie_iTBS];
end


end

%22 July 2020.
% Now compare with a fit to the average curve.
disp(' ')
disp('Fitting the average curve, not per person')
person=person+1; %increment again
[wee_base, wie_base, lsq_base, fittedamp, RMT_actually_used, output_curve_base]=fit_weights_to_meps_rel_oct22(mean(base_s_ave,1),mean(rmt_list),amp_to_use,stunc_per_person,amp_or_squash,goldsworthy,person,0.9,0.9,base_null,mep_null,fix_ee,fix_ie,S_select_para);
before_single_wee=wee_base; before_single_wie = wie_base; before_single_lsq = lsq_base;

amp_for_this = fittedamp;   %load up the amplitude to use for this person.




%The fitting for the TBS

disp('TBS')
[wee_iTBS, wie_iTBS, lsq_iTBS, fittedamp, RMT_actually_used, output_curve_TBS]=fit_weights_to_meps_rel_oct22(mean(iTBS_s_ave,1),mean(rmt_list),amp_for_this,stunc_per_person,amp_or_squash,goldsworthy,person,wee_base,wie_base,output_curve_base,mean(base_s_ave,1),fix_ee,fix_ie,S_select_para);  %20 July 2020. Fit with this amplitude. 10jan21 bring in the fitted baseline weights too
%THIS WAS A DEBUG CALL ONLY [wee_iTBS, wie_iTBS, lsq_iTBS, fittedamp]=fit_weights_to_meps(mean(base_s_ave,1),mean(rmt_list),amp_for_this,stunc_per_person,amp_or_squash,goldsworthy);  %20 July 2020. Fit with this amplitude

after_single_wee=wee_iTBS; after_single_wie = wie_iTBS; after_single_lsq=lsq_iTBS;





%Now show means
disp('  ')
disp('Only the accepted people included  ')
disp(['number of accepted people is ' num2str(length(wee_b))])
disp(['average wee before is ' num2str(mean(wee_b)) ' and ave wie before is ' num2str(mean(wie_b))])
disp(['average wee after is ' num2str(mean(wee_i)) ' and ave wie after is ' num2str(mean(wie_i))])



%plot graph for accepted people only
figure(1)
plot(wee_b,wie_b,'bo'); %initial values
grid on; hold on; pbaspect([1 1 1]); set(gca,'xlim',[0.75 1.05]); set(gca,'ylim',[0.75 1.05]);
plot(wee_i,wie_i,'rx'); %final values
for i=1:length(wee_b)
    plot( [wee_b(i) wee_i(i)],[wie_b(i) wie_i(i)],'k-')
end
xlabel('w_{ee} multiplier')
ylabel('w_{ie} multiplier')
title('blue-before, red-after, accepted people')

%plot graph with lines scaled according to the lsq weights. Use area of
%circile proportional to weight.
figure(2)
grid on; hold on; set(gca,'xlim',[0.75 1.05]); set(gca,'ylim',[0.75 1.05]); pbaspect([1 1 1]);

for i=1:length(wee_base_array)
   plot(wee_base_array(i),wie_base_array(i),'bo','MarkerSize',(60/sqrt(lsq_base_array(i) + lsq_iTBS_array(i)))  ); %initial values
plot(wee_iTBS_array(i),wie_iTBS_array(i),'rx','MarkerSize',(60/sqrt(lsq_base_array(i) + lsq_iTBS_array(i)))  ); %final values
    plot( [wee_base_array(i) wee_iTBS_array(i)],[wie_base_array(i) wie_iTBS_array(i)],'k-','linewidth', (20/(lsq_base_array(i) + lsq_iTBS_array(i)) )   )
end
xlabel('w_{ee} multiplier')
ylabel('w_{ie} multiplier')
title('blue-before, red-after, all people lsq weighted')

%Weight each point according to lsq (sum of both base and itBS)
wee_before = sum( wee_base_array.*(lsq_base_array + lsq_iTBS_array).^(-1) )/sum( (lsq_base_array + lsq_iTBS_array).^(-1)) ;
wie_before = sum( wie_base_array.*(lsq_base_array + lsq_iTBS_array).^(-1) )/sum( (lsq_base_array + lsq_iTBS_array).^(-1));

wee_after = sum( wee_iTBS_array.*(lsq_base_array + lsq_iTBS_array).^(-1) )/sum( (lsq_base_array + lsq_iTBS_array).^(-1));
wie_after = sum( wie_iTBS_array.*(lsq_base_array + lsq_iTBS_array).^(-1) )/sum( (lsq_base_array + lsq_iTBS_array).^(-1));
%Weighted according to the lsq sum output
disp('  ')
disp('All people using the inverse sum of squares as a weighting factor ')

disp(['average wee before is ' num2str(wee_before) ' and ave wie before is ' num2str(wie_before)])
disp(['average wee after is ' num2str(wee_after) ' and ave wie after is ' num2str(wie_after)])


disp(' ')
disp('Weights for fit to average MEP data across all people and sessions')
disp(['average wee before is ' num2str(before_single_wee) ' and ave wie before is ' num2str(before_single_wie) ' with lsq ' num2str(before_single_lsq)])
disp(['average wee after is ' num2str(after_single_wee) ' and ave wie after is ' num2str(after_single_wie) ' with lsq ' num2str(after_single_lsq)])


%13 August 2020
%Present in terms of amount moved. Put all befores at 1.0,1.0
%wee_b and wie_b are 'before' values,  wee_i and wie_i are the after ones
ee_movement=wee_i-wee_b;
ie_movement=wie_i-wie_b;
figure(3); 
plot(ee_movement,ie_movement,'kx'); hold on
xlabel('w_{ee} movement')
ylabel('w_{ie} movement')
grid on;  pbaspect([1 1 1]);
title('Net movement in w_{ee} and w_{ie}')


%if it's first time through, create the ee_move_list
if (time_period==time_period_min)
    ee_move_list=[0*ee_movement; ee_movement];
    ie_move_list=[0*ie_movement; ie_movement];
    traj_of_mean_ee=[0; after_single_wee-before_single_wee];
    traj_of_mean_ie=[0; after_single_wie-before_single_wie];
    diff_in_time=0*ee_movement;
else
    %if not first time, add them to list
    ee_move_list=[ee_move_list; ee_movement];  %update array
    ie_move_list=[ie_move_list; ie_movement];
    traj_of_mean_ee=[traj_of_mean_ee; after_single_wee-before_single_wee];
    traj_of_mean_ie=[traj_of_mean_ie; after_single_wie-before_single_wie];
end



%Now look at the ee ie difference movement
diff_movement=ee_movement-ie_movement;
figure(4); 
binedges=[-0.21:0.02:0.21];
subplot(4,1,time_period)
numbers=histogram(diff_movement,binedges);  %put the differential movement into bins and plot it
title(['time period ' num2str(time_period)])

diff_in_time=[diff_in_time; diff_movement];
xlabel('movement in w_{ee}-w_{ie} from baseline')
ylabel('frequency')

end   %end time period

time_slices=[0 [time_period_min:time_period_max] ]';
for i=1:length(ee_move_list)  %length of people list
    figure(3)
    plot(ee_move_list(:,i),ie_move_list(:,i),'k-')  %draw a connecting line for the time periods
    figure(5)
    plot3(time_slices,ee_move_list(:,i),ie_move_list(:,i)); hold on; grid on
    xlabel('time period')
    ylabel('w_{ee} movement')
    zlabel('w_{ie} movement')
end
%plot mean trajectory
mean_wee=mean(ee_move_list,2);  %average over columns
mean_wie=mean(ie_move_list,2);
plot3(time_slices,mean_wee,mean_wie,'k-','linewidth',3) %mean trajectory
plot3(time_slices,0*time_slices,0*time_slices,'k--','linewidth',2)
title([name_of_data ' changes in w_{ee} and w_{ie} for all people'])

%plot trajectory of mean
plot3(time_slices,traj_of_mean_ee,traj_of_mean_ie,'bx-','linewidth',3) %trajectory of mean


figure(4)
subplot(4,1,4)
title('movement in time')
plot(time_slices',diff_in_time); grid on; hold on
%Now plot the trajectory of the average, and the average trajectory
%Plotting the mean trajectory
plot(time_slices',mean_wee-mean_wie,'k-','linewidth',2); 
plot(time_slices',0*time_slices','k--','linewidth',2) %Plot out a zero line
%Plot the trajectory of mean
plot(time_slices',traj_of_mean_ee-traj_of_mean_ie,'bx-','linewidth',2)
xlabel('time period (0=baseline)')
ylabel('change in w_{ee}-w_{ie}')
set(gca,'ylim',[-0.25 0.25])



remove_some=false
if (remove_some)
%Remove people 3, 5 from iTBS. Replace with ave of others
keepsies=[1:10 12:18];  %people to keep
base_s_ave=base_s_ave(keepsies,:);
iTBS_s_ave=iTBS_s_ave(keepsies,:);
output_curve_base_in_array=output_curve_base_in_array(keepsies,:);
output_curve_TBS_in_array=output_curve_TBS_in_array(keepsies,:);
end


%16 August 2021

%Now take modelled values and pull out a curve for them

%We have the model curves in the form [90:10:180]

%Now just plot the modelled meps

figure
amps_in_rmt=[90:10:180];
plot_IO(amps_in_rmt,output_curve_TBS_in_array,+1,'b-',1);



%16 August 2021
%Now take experimental values and pull out a curve for them

plot_IO(amps_in_rmt,output_curve_base_in_array,-1,'k-',1);


% figure
% %fitted relative to expt plot
% % 
% %first the post TBS one
% relative_post_TBS=output_curve_TBS_in_array./iTBS_s_ave;   %divide the best fit model curve by actual measurement
% plot_IO(amps_in_rmt,relative_post_TBS,1,'b-',1);
% %
% %Then the pre TBS (baseline) one
% relative_baseline=output_curve_base_in_array./base_s_ave;
% plot_IO(amps_in_rmt,relative_baseline,-1,'k-',1);
% ylabel('model relative to baseline')

figure
to_plot=[1:18];  %the people to put on the plot. Sometimes we have problems here
%Plot the post/pre for each of the expt and model conditions
%After Nigel's email 17 August 2021
expt_post_over_pre = iTBS_s_ave(to_plot,:)./base_s_ave(to_plot,:);
model_post_over_pre = output_curve_TBS_in_array(to_plot,:)./output_curve_base_in_array(to_plot,:);
plot_IO(amps_in_rmt,expt_post_over_pre,-1,'kx-',1)
plot_IO(amps_in_rmt,model_post_over_pre,1,'bo-',1)


figure
%expt results
plot_IO(amps_in_rmt,base_s_ave,-1,'k-',1);
plot_IO(amps_in_rmt,iTBS_s_ave,+1,'b-',1);
title('experimental results; bars are stdevs')



%17 August 2021. Look at PCA for the post/pre ratio for expt, then describe
%modelling fits with it. 
%z=zscore(iTBS_s_ave./base_s_ave);    %ratio of expt postTBS to expt baseline 
z=zscore(expt_post_over_pre)
[coeff score latent]=pca(z);   
pc1_expt=score(:,1);   %the first pcas. Describes the general deviation from the trend for each individual. 
disp('relative percent contributions are ')
disp(100*latent./sum(latent))

%Now look at the model results
%zm=zscore(output_curve_TBS_in_array./output_curve_base_in_array);   %the ratio of the modeled fit postTBS to the model fit baseline
zm=zscore(model_post_over_pre);   %use the calculated value from earlier

%page 7521.  The scores for the model in the expt PCs basis are 
% score_model = zmodel*(coeff')^-1 = Vmodel*coeff;
score_model=zm*coeff;     %just find those scores. 18 x 10 matrix.
pc1_model=score_model(:,1);  %first column

figure
plot(pc1_expt,pc1_model,'kx')
xlabel('PC1 experiment')
ylabel('PC1 model')
title('person-by-person variation in relative MEP post/pre')
hold on; grid on;
%fit trendline
p=polyfit(pc1_expt,pc1_model,1);
xfits=[min(pc1_expt) max(pc1_expt)];   %the x-values for the line
yfits=xfits*p(1)+p(2);   %the yfits
plot(xfits,yfits,'k--')
yfitted_all=pc1_expt*p(1)+p(2);   %fit to all points
R2=1 - sum((pc1_model-yfitted_all).^2)/sum((pc1_model-mean(pc1_model)).^2);   %the R2 value
title(['R-squared=' num2str(R2)])

figure
%19 August 2021. Show post/pre for some different people
PC_vals=pc1_expt+pc1_model;
[M,I3]=max(PC_vals); %pick index of maximum
[M,I1]=min(PC_vals); %pick index list of minimum
[M,I2]=min(abs(pc1_expt)+abs(pc1_model)); %pick a near zero case

peoples=[I1 I2 I3]
colour_list='krb'; %order of colours
for j=1:length(peoples)
    linekey=[colour_list(j) '--']; %put a dashed line indicator on the end 
plot(amps_in_rmt,expt_post_over_pre(peoples(j),:),linekey); %measured
grid on; hold on;
xlabel('amplitude (pc RMT)')
ylabel('MEP (mV)')
set(gca,'xlim',[85 185])
    linekey=[colour_list(j) '-']; %put a solid line indicator
plot(amps_in_rmt,model_post_over_pre(peoples(j),:),linekey); %modelled
end
title(['Individuals, low (black, person ' num2str(I1) '), medium, (red, person ' num2str(I2) '), high, (blue, person ' num2str(I3) ')'])


%3 September 2021
%Record the wee, wie and S values that we get out. 

figure
subplot(3,1,1)
disp('Change in wee:')
(wee_iTBS_array - wee_base_array)'
[h,p]=ttest(wee_iTBS_array,wee_base_array)
disp('')
histogram(wee_iTBS_array-wee_base_array,binedges)
xlabel('Movement in wee')
ylabel('Number of subjects')

disp('Change in wie:')
subplot(3,1,2)
(wie_iTBS_array - wie_base_array)'
[h,p]=ttest(wie_iTBS_array,wie_base_array)
disp('')
histogram(wie_iTBS_array-wie_base_array,binedges)
xlabel('Movement in wie')
ylabel('Number of subjects')


disp('Change in amp factor')
subplot(3,1,3)
(amp_after_TBS - amp_for_this_person)'
[h,p]=ttest(amp_after_TBS,amp_for_this_person)
disp('')
histogram(amp_after_TBS.^(-1)-amp_for_this_person.^(-1),binedges)
xlabel('Movement in \mu')
ylabel('Number of subjects')


%Now look for changes in principal components. 3 September 2021
before_list=[wee_base_array' wie_base_array'  amp_for_this_person' ];
after_list=[wee_iTBS_array' wie_iTBS_array' amp_after_TBS'];
z=zscore([before_list; after_list]);  %do the lot as zscore
[coeff,score,latent]=pca(z)   %do a pca on these
figure; 
plot(score(1:18,1),score(19:36,1),'kx'); grid on; hold on;
plot([min(score(1:18)) max(score(1:18))],[min(score(1:18)) max(score(1:18))] , 'k--')
xlabel('PC1 baseline')
ylabel('PC1 after TBS')
title('Changes in S, wee and wie as PCs')

figure(99)
%here we plot the change in wee-wie against pc1expt

delta=(wee_iTBS_array-wie_iTBS_array)' - (wee_base_array-wie_base_array)';  %after - before. Change in net excitation
%now plot against PC1
plot(pc1_expt(1:18),delta(1:18),'kx'); grid on; hold on
xlabel('PC1 expt')
ylabel('change in wee-wie')
p=polyfit(pc1_expt(1:18),delta(1:18),1);
xfits=[min(pc1_expt(1:18)) max(pc1_expt(1:18))];   %the x-values for the line
yfits=xfits*p(1)+p(2);   %the yfits
plot(xfits,yfits,'k--')
yfitted_all=pc1_expt(1:18)*p(1)+p(2);   %fit to all points
R2=1 - sum((delta(1:18)-yfitted_all).^2)/sum((delta(1:18)-mean(delta(1:18))).^2);   %the R2 value
title(['R-squared=' num2str(R2)])


%29 October 2021
%Do some parallel plots
figure(60)
subplot(2,2,2)
pre=wee_base_array'; post=wee_iTBS_array';
tbl=table(pre,post)
p=parallelplot(tbl,'DataNormalization','none');
grid on
title('(b) cTBS, wee');

% subplot(2,2,3)
% pre=wie_base_array'; post=wie_iTBS_array';
% tbl=table(pre,post)
% p=parallelplot(tbl,'DataNormalization','none');
% grid on
% title('(c) iTBS, wie');

% subplot(4,2,6)
% pre=wee_base_array'-wie_base_array'; post=wee_iTBS_array'-wie_iTBS_array';
% tbl=table(pre,post)
% p=parallelplot(tbl,'DataNormalization','none');
% grid on
% title('(f) cTBS, wee-wie');
% 
 subplot(2,2,4)
 pre=amp_for_this_person'; post=amp_after_TBS';
 tbl=table(pre,post)
 p=parallelplot(tbl,'DataNormalization','none');
 grid on
 title('(d) cTBS, 1/\mu');


%29 September 2021
%Make up a correlation plot of various thing
[m,c,R2]=fit_and_correlate(pc1_expt,wee_iTBS_array-wee_base_array,[1:18],'PC1 experiment','\Delta w_{ee}','(a) iTBS S+ee',1);
%[m,c,R2]=fit_and_correlate(pc1_expt,wie_iTBS_array-wie_base_array,[1:18],'PC1 experiment','\Delta w_{ie}','(b)',2);
[m,c,R2]=fit_and_correlate(pc1_expt,amp_after_TBS-amp_for_this_person,[1:18],'PC1 experiment','\Delta (1/\mu)','(b)',2);
%[m,c,R2]=fit_and_correlate(pc1_expt,(wee_iTBS_array-wie_iTBS_array)-(wee_base_array-wie_base_array),[1:18],'PC1 experiment','\Delta (w_{ee}-w_{ie})','(d)',4);


plot_big_figure_for_paper;  %11 January 2022. Plot the nice composite figure at the end