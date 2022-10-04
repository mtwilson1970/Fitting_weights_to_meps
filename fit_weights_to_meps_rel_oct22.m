function [wee_out, wie_out, lsq_out, amp_factor, RMT, my_sampled]=fit_weights_to_meps_rel_oct22(mep_value_list,rmt_for_this_person,chosen_amp,stunc_per_person,amp_or_squash,goldsworthy,person,wee_in,wie_in,base_normalize,expt_normalize,fix_ee,fix_ie,S_select_para);



%10 feb 2022. When we look at evaluating least squares we divide by the
%base_normalize coming in. First time it's all ones, second time it's the baseline
%curve.
%expt_normalize is the experimental curve second time around, or all ones
%the first time

%16 August 2021. Spit out the RMT used now.

%11 August 2020. Chosen_amp could be a squash on the x-axis to. That's what
%the amp_or_squash flag does. True is the mep (y) stretch, false is
%the TMS amplitude (x) squash.

%17 July 2020, p7041. Weight array as a function

% 20 July. Allow it to fit an amplitude as well if we want. Output paras
% are wee, wie and the fitted amplitude.

%chosen_amp in the input is the multiplicative scale we need to use - if
%zero then we find it and report as amp_factor.

%stunc_per_person gives the standard uncertainty for one person in the mep.
%Helps weight the fitting. If we have it 0 then we fit a specified stunc
%from the Goldsworthy iTBS plot.


% weight_array.m
% See p 6952. Calculates an IO curve for an array of w_ee and w_ie values
%
% Marcus Wilson
% 29 May 2020
%
  

%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%

just_do_plots = true;   
    %if true (standard) we skip evaluating the IO curves. 
    %If false we rather than evaluate IO curves using neurofield

actually_show_plots = false;  
    %If true we show lots of subisdary plots - useful for debugging.
    %If false we don't.
    %Plotting takes a long time so standard is false. 

load_up_matrix = false;   %collect together a single file with all results in it
    %If true we compile and save the calculated results into a look-up table
    %If false (standard) we don't save them.
file_to_save_lookup_table_in='big_matrix_of_results.mat'
    %This is the filename where the look-up table will be saved (if load_up_matrix=true)
    
use_matrix = true;  
    %If true (standard) we will use a previously calculated look-up table. 
    %If false we won't
file_with_results_in='enormous_interpolated_results_p7122.mat'
    %This is the filename of the look-up table (used if use_matrix=true)

use_rmt_value = true;  %use a selected value for the RMT amplitude
    %If true (standard) use a previously-calculated RMT value
    %If false calculate the RMT for each IO curve. (Prone to error - not
    %recommended). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







  
  %What do we fix in terms of ee and ie?
  if (fix_ee==0)   %We will fit it both pre- and post-. Release it.
      lockee=0;
  end
  if (fix_ee==1)  %We fit it to neither pre- nor post-. Lock it. 
      lockee=1;
  end
  if (fix_ee==2)  %We fit it to pre- but not post-
      if (chosen_amp==0)  %first time through
          lockee=0; 
      else  %second time through
          lockee=1;   %lock it
      end
  end
  
  if (fix_ie==0)   %We will fit it both pre- and post-. Release it.
      lockie=0;
  end
  if (fix_ie==1)  %We fit it to neither pre- nor post-. Lock it. 
      lockie=1;
  end
  if (fix_ie==2)  %We fit it to pre- but not post-
      if (chosen_amp==0)  %first time through
          lockie=0; 
      else  %second time through
          lockie=1;   %lock it
      end
  end
  
  
  
  %do a calc of scale factor each time.    Choose the right 'chosen_amp' value. See page 7535 to do a S fit or (S) fit to the model. 
   if (S_select_para==0);
      chosen_amp=0;   %Here we fit S both to baseline and post-iTBS independently
   end
   if (S_select_para==1);
       chosen_amp=1;  %Here we don't fit S at all, to either case
   end
   if (S_select_para==2);
       %Here we don't respecify chosen_amp. It will be whatever value it
       %comes in to the routine as. This means we fit S the first time
       %through (post TBS) but not the second. 
   end
  %0 = do an S fit both to baseline and post-TBS independently
%1 = don't fit S at all, either to baseline or post-TBS
%unspecified (i.e. comment the line out) is fit the baseline and use the same value for the post-TBS

  

%standard values of ee, ie; might be overwritten
   ee_mult=[0.75:0.01:1.05];  standee=16; 
   ie_mult=[0.75:0.01:1.05];  standie=16;

  
  
  if (lockee==1)
      ee_mult=ones(1,31)*wee_in;     %make it all the same 
      locked_i = round((wee_in - 0.90)/0.01) + 16    %find the corresponding index. This is the one to look up in matr
  end
      
  if (lockie==1)
      ie_mult= ones(1,31)*wie_in;    %here are our ie_mult values to come out
      locked_j = round((wie_in - 0.90)/0.01) + 16    %find the corresponding index that's been locked
  end
      

  
  
  










amp_factor=NaN;  %null to start with.





if (goldsworthy==1)  %do we use the goldsworthy (true) or vallence (false) values for rmt in terms of s-1 rate.
rmt_set_keep =  rmt_for_this_person*620.0/42.28; %Goldsworthy iTBS data. in this case it's this value here in /s. Otherwise this line is ignored. 609 is average
else
    if (goldsworthy==3)
     % rmt_set_keep =  rmt_for_this_person*653.4/48.1;  %Vallence all 3 sessions cTBS data, see p7055. Assume average is 50
                   rmt_set_keep = rmt_for_this_person*623.6/39.06; %Vallence for just session 1.
    % disp('****WARNING. Using wrong RMT values. Need to correct***'); %Vallence for just session 2.    
             
    else
        if (goldsworthy==2.1)
            rmt_set_keep = rmt_for_this_person*623.6/39.06; %Vallence for just session 1.
        else
            if (goldsworthy==2.2)
             % rmt_set_keep = rmt_for_this_person*657.9/48.1111;
                rmt_set_keep = rmt_for_this_person*623.6/39.06; %Vallence for just session 1.

       %      disp('****WARNING. Using wrong RMT values. Need to correct***'); %Vallence for just session 2.    
              %stop
            else
              %  rmt_set_keep=rmt_for_this_person*654.8/48.1111;  %vallence for just session 3
               rmt_set_keep = rmt_for_this_person*623.6/39.06; %Vallence for just session 1.

         %     disp('****WARNING. Using wrong RMT values. Need to correct***');
              %stop
            end
        end
    end
end

master_name = 'IO_master_conf_file.conf';  %a master file with standard parameters. DON'T EDIT THIS ONE.
new_name_base = 'person_1';  %the copied file basename which can be edited.
new_name = [new_name_base '.conf'];  %complete file name
new_name_temp_c = [new_name_base '_after_ee_change.conf'];
new_name_finished = [new_name_base '_after_ie_change.conf'];
end_name_base = [new_name_base '_for_saving_grand_matrix'];

weight_ee_standard = 1.92e-4;
weight_ie_standard = 1.92e-4;  %standard couplings
amp_values=[400:40:1200]; %a row vector of amplitude values; Must be integers




flagfirst=chosen_amp;   %9sept20. Copy it before it gets overwritten

%ee_mult=[0.91 0.92 0.93]; standee = 2;   %multiplier for the ee weights; hich index is 'standard'
%ie_mult=[0.87 0.88 0.89]; standie = 2; %multiplier for ie; which index is 'standard'

%ee_mult=[0.86 0.89 0.92 0.95 0.98] ; standee = 3;   %multiplier for the ee weights; hich index is 'standard'
%ie_mult=[0.83 0.86 0.89 0.92 0.95]; standie = 3; %multiplier for ie; which index is 'standard'

% ee_mult=[0.82:0.02:0.98];  standee=5;   %p7028 Grand matrix of curves
% ie_mult=[0.82:0.02:0.98];  standie=5;   
% 
%  %p7062  enormouse matrix
%  ee_mult=[0.82:0.01:0.98]; standee=9;
%  ie_mult=[0.82:0.01:0.98]; standie=9;


if (~just_do_plots)    %if we're not just plotting, we need to calculate

%create a conf file from the master
[status_of_cp]=create_person(master_name, new_name);  %creates a copy of the conf file with the new name.


for i=1:length(ee_mult); %go over all ee values
    weight_ee = weight_ee_standard * ee_mult(i); 
    % 
    % Now write a new conf file
    
    [status] = change_weights_ee(new_name,new_name_temp_c,weight_ee);  %changes the weights for the ee
    
    for j=1:length(ie_mult); %go over all ie values
        weight_ie = weight_ie_standard * ie_mult(j);  %pick the weight we need
        
        %write the conf file
       
        [status] = change_weights_ie(new_name_temp_c,new_name_finished,weight_ie);
        
        %Now we can run an IO curve

        [IOcurve_baseline,dummycouples_ee,dummy_couples_ie]=run_IO_curve(amp_values,new_name_finished);  %constructs an IO curve with amps given by ampvalues and puts it in IO_curve_result
        
        %Now save results
        name_to_save_under = [end_name_base num2str(i) num2str(j) '.mat'];   %name to save the data under
        save(name_to_save_under);  %save the values
    
    end   %next ie value
    
    
end    %next ee value

end   %end the just doing plots


%Now we show results
%close all;

if (load_up_matrix)
    matr=zeros( length(ee_mult), length(ie_mult), length(amp_values) ); % a blank matrix
end

if (use_matrix)
%   load 'big_matrix_of_results.mat'
   %  load 'enormous_interpolated_results.mat'
     load(file_with_results_in);   %5 oct 22 - specify the name in the optiosn
     matr=bigger;   %put in the p7062 29july2020 interpolated matrix
end


k=0;  %subplot number
%find the standard curve
if (use_matrix)
    standard_out=matr(standee,standie,:);  %if we are using a preloaded matrix
else
load([end_name_base num2str(standee) num2str(standie) '.mat']);
standard_out = IOcurve_baseline;   %pick out the standard response
end



%Goldsworthy cTBS data 2016
amps_in_rmt=[90:10:180];  amps_in_rate=amps_in_rmt*625/100;

if (stunc_per_person(end)==0);   %22 July 2020, page 7050. Specify the uncertainty
    
     sem_mep = [0 0.03 0.05 0.20 0.29 0.38 0.42 0.45 0.46 0.46];  %from goldsworthy. Use these
else
    sem_mep=stunc_per_person;
end

if (chosen_amp == 0)
    %if the amp selected is zero, we need to calculate it
    if (amp_or_squash)
        %if on the yaxis
        amp_min = 0.05; amp_max = 10.0;   %we'll search these limits
    else
        %on the x-axis
        amp_min=0.5; amp_max = 2;
        %amp_min=1.0;  amp_max=1.0;
    end
    
else
    %otherwise use selected one
    amp_min = chosen_amp; amp_max = amp_min;   %just use the chosen one

end




for i=1:length(ee_mult)
    for j=1:length(ie_mult)
        k=k+1; %increment subplot
        
        
             if (use_matrix)  %if we use the big matrix
                 
                 %11 january 2022. Need to use the correct values of i and
                 %j. They might be locked.            
                 if (lockee)
                     i_to_use=locked_i;
                 else
                     i_to_use=i;   %otherwise use the value i
                 end
                 
                 %11 january 2022. Need to use the correct values of i and
                 %j. They might be locked.            
                 if (lockie)
                     j_to_use=locked_j;
                 else
                     j_to_use=j;   %otherwise use the value i
                 end

                IOcurve_baseline(1:length(amp_values)) = matr(i_to_use,j_to_use,:);
             else
                load([end_name_base num2str(i) num2str(j) '.mat']);  %load the file
             end
        
        
        if (load_up_matrix)
           matr(i,j,:) = IOcurve_baseline(:); % load up matrix
        end
        
  

        
        
        
 if (actually_show_plots)       
        figure(1)
         subplot(length(ee_mult),length(ie_mult),k);  %select subplot
         plot(amp_values,IOcurve_baseline,'k-'); hold on; grid on;
         plot(amp_values,standard_out,'k--');  %standard response
         set(gca,'xlim',[400 1200]);
         set(gca,'ylim',[0 5]);
 end
 
 
       
         %find the %RMT value for the xaxis. Need to interpolate. 
       
        if (use_rmt_value)
            RMT=rmt_set_keep;  %set it to user defined.
        else
            min_mep = 0.20;  %threshold for calling it a MEP
      
            [RMT]=find_rmt(min_mep,amp_values,IOcurve_baseline); %finds the RMT for this plot
        end
   
 if (actually_show_plots)       
         figure(2)
         subplot(length(ee_mult),length(ie_mult),k);  %select subplot
     %plot w.r.t. rmt
       plot(100*amp_values/RMT,IOcurve_baseline,'k-');grid on; hold on
       plot(amps_in_rmt,mep,'rx-');   %goldsworthy
        plot(amps_in_rmt,mep+sem_mep,'r--');
        plot(amps_in_rmt,mep-sem_mep,'r--');  %with standard deviations
        set(gca,'xlim',[80 200])
        set(gca,'ylim',[0 5])
 end
 
 %interpolate my values to the Goldsworthy values;
      
%What we do here depends on whether we stretch y-axis or x-axis
 if (amp_or_squash)
     %if true we multipuly the amplitudes (yaxis). Need to sample first
my_sampled = interp1(100*amp_values/RMT,IOcurve_baseline,amps_in_rmt,'linear','extrap');  %extrapolate if necessary 27 July 2020



 end
       
        
        %How good is it?  20 July 2020. Run over amp values
       leastsq=9999999999;  %initial value to be overridden
       chosen_amp=[];
       for chosen_amp = amp_min:0.01:amp_max;   %go from min to max in 0.01 steps
        
           if (amp_or_squash)
           
             mep=mep_value_list*chosen_amp;  %load up the mep and multiply by chosen amp
           else
           %if we are squashing x-axis, need to do a sample based on chosen
           %amp
   
           %This one for scaling the axis 27aug20         
      my_sampled = interp1(100*amp_values/RMT,IOcurve_baseline,amps_in_rmt*chosen_amp,'linear','extrap');  %extrapolate if necessary 27 July 2020
      
      
      
   %this one for shifting the position of the xaxis %27aug20
   %my_sampled = interp1(100*amp_values/RMT,IOcurve_baseline,amps_in_rmt*(chosen_amp-1.25)*10,'linear','extrap');  %extrapolate if necessary 27 July 2020

      mep=mep_value_list;
             
      
          
             
           end
           
     %if false we scale the x-axis
             
          leastsq_forthisamp=sum( ((  my_sampled(1:10)./base_normalize -mep(1:10)./expt_normalize).^2)./[ sem_mep(1:10)].^2   );  %normalize the assessed fit by the baseline curve
  %      leastsq_forthisamp=sum( ((my_sampled(1:5)-mep(1:5)).^2)./[sem_mep(2) sem_mep(2:5)].^2   );    %the first five points
       
        
        if (leastsq_forthisamp < leastsq);  
                 %if its smaller than current smallest, update
                 leastsq=leastsq_forthisamp;
                 best_amp_so_far = chosen_amp;   
        end
          
      
        
        
        
        end   %end the loop for amp. Now have a best_amp_so_far and a leastsq value
        leastsq_matrix(i,j)=leastsq;   %load up result in a matrix
        best_amp_for_this_ij(i,j)  = best_amp_so_far;
      

        
        
        if (actually_show_plots)
        
            if (i==length(ee_mult))  %last row
               figure(1); xlabel('amplitude (/s)')
                figure(2); xlabel('amplitude (\%RMT)','interpreter','latex')
            end
            if (j==1);  %on first of each row
               figure(1); ylabel('MEP (mV)')
                figure(2); ylabel('MEP (mV)')
            end
            figure(1);
            title( ['$w_{ee} \times$' num2str(ee_mult(i)) ' $w_{ie} \times$' num2str(ie_mult(j))   ], 'interpreter','latex')
            figure(2);
            title( ['$w_{ee} \times$' num2str(ee_mult(i)) ' $w_{ie} \times$' num2str(ie_mult(j)) ' lsq' num2str(leastsq)  ], 'interpreter','latex')

        end %end actually_show_plots
        
    end
    
    
    
    
end
  
 
        
if (load_up_matrix)
    %save results
    save(file_to_save_lookup_table_in)
end

%leastsq_matrix

%find minimum of leastsq_matrix and its position. 16 July 2020
pointee=[];
pointie=[];
min_of_lsq=999999;  %nulls to start with
for i=1:length(ee_mult);
    for j=1:length(ie_mult);
        if  (leastsq_matrix(i,j) < min_of_lsq)  %we have a new minimum found
            pointee=i; pointie=j;  min_of_lsq=leastsq_matrix(i,j);   %set the new minimum
       %pointee 
       %pointie
        end
    end
end
disp(['Best fit position is ee=' num2str(ee_mult(pointee)) 'x and ie=' num2str(ie_mult(pointie)) 'x with lsq=' num2str(min_of_lsq)])
disp(['Multiplicative amplitude is ' num2str(best_amp_for_this_ij(pointee,pointie))])

wee_out=ee_mult(pointee);
wie_out=ie_mult(pointie);
lsq_out=min_of_lsq;
amp_factor=best_amp_for_this_ij(pointee,pointie);  %19 August 2020. Pick the best amp for the best position


%figure(7)
%contourf(ee_mult,ie_mult,leastsq_matrix'/sqrt(10),[0:2:20]); colorbar
%stop


%%%%%9 Sept 2020%%%% Now plot the graph we want for the best fit only

i=pointee;
j=pointie;   %select best point from what we've evaluated. Some of the points might be locked

if (lockee==1)
    i_to_use=locked_i;
else
    i_to_use=i;
end
if (lockie==1)
    j_to_use=locked_j;
else
    j_to_use=j;
end
       
  amp_min=best_amp_for_this_ij(i,j); amp_max=amp_min ;  %select amps
        
             if (use_matrix)  %if we use the big matrix

                IOcurve_baseline(1:length(amp_values)) = matr(i_to_use,j_to_use,:);   %use the locked values
             else
                load([end_name_base num2str(i) num2str(j) '.mat']);  %load the file
             end
        
        
        if (load_up_matrix)
           matr(i,j,:) = IOcurve_baseline(:); % load up matrix
        end
            
       
       
         %find the %RMT value for the xaxis. Need to interpolate. 
       
        if (use_rmt_value)
            RMT=rmt_set_keep;  %set it to user defined.
        else
            min_mep = 0.20;  %threshold for calling it a MEP
      
            [RMT]=find_rmt(min_mep,amp_values,IOcurve_baseline); %finds the RMT for this plot;
            %this needs recording
        end
   

 
 %interpolate my values to the Goldsworthy values;
      
%What we do here depends on whether we stretch y-axis or x-axis
 if (amp_or_squash)
     %if true we multipuly the amplitudes (yaxis). Need to sample first
my_sampled = interp1(100*amp_values/RMT,IOcurve_baseline,amps_in_rmt,'linear','extrap');  %extrapolate if necessary 27 July 2020
 end
       
        
        %How good is it?  20 July 2020. Run over amp values
       leastsq=999999;  %initial value to be overridden
       chosen_amp=[];
       for chosen_amp = amp_min:0.01:amp_max;   %go from min to max in 0.01 steps
        
           if (amp_or_squash)
           
             mep=mep_value_list*chosen_amp;  %load up the mep and multiply by chosen amp
           else
           %if we are squashing x-axis, need to do a sample based on chosen
           %amp
   
           %This one for scaling the axis 27aug20         
      my_sampled = interp1(100*amp_values/RMT,IOcurve_baseline,amps_in_rmt*chosen_amp,'linear','extrap');  %extrapolate if necessary 27 July 2020
   
   %this one for shifting the position of the xaxis %27aug20
   %my_sampled = interp1(100*amp_values/RMT,IOcurve_baseline,amps_in_rmt*(chosen_amp-1.25)*10,'linear','extrap');  %extrapolate if necessary 27 July 2020
   
      %16 August 2021;
   
   
   
      mep=mep_value_list;                   
             
           end
           
     %if false we scale the x-axis
             figure(6);  %now plot this out
             subplot(4,5,person)
      
             if (flagfirst==0)   %if we are first time through
        plot(amps_in_rmt,my_sampled,'m--');grid on; hold on
        plot(amps_in_rmt,mep,'k--')  %expt line as well
       
        for z=1:length(amps_in_rmt)
        %displace in x slightly so no overlap with the next colour
            plot([amps_in_rmt(z)+1 amps_in_rmt(z)+1 amps_in_rmt(z)+1],[mep(z)-sem_mep(z) mep(z) mep(z)+sem_mep(z)], 'k-','linewidth',0.5);   %goldsworthy
        end
    
        
         else   %(else plot in blue)
           % plot(amps_in_rmt,my_sampled,'b-');grid on; hold on
       
           
        for z=1:length(amps_in_rmt)
      %  plot([amps_in_rmt(z)-1 amps_in_rmt(z)-1 amps_in_rmt(z)-1],[mep(z)-sem_mep(z) mep(z) mep(z)+sem_mep(z)], 'b-','linewidth',0.5);   %goldsworthy
        end 
     
        %Now plot the lsq landscape
        figure(7)
        contourf(ee_mult,ie_mult,leastsq_matrix',[0:2:20]); colorbar
        xlabel('w_{ee} multiplier')
        ylabel('w_{ie} multiplier')
        title('fit landscape for average person')
        
        
         end
             
        figure(6)
        set(gca,'xlim',[80 190])
        set(gca,'ylim',[0 6])
        if (flagfirst==0)
        if (person<19)  %i.e. 18 or below - all people
            title(['pn ' num2str(person) ', S=' num2str(leastsq_matrix(i,j)) ])
        else
            title(['average pn, S=' num2str(leastsq_matrix(i,j))])
        end
        end
                    
                
                xlabel('percent RMT')          
                 ylabel('MEP (mV)')    
                 
     end   %end the loop for amp. Now have a best_amp_so_far and a leastsq value
         
           
     
   


end







