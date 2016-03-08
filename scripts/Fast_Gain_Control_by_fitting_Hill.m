% 
% 
% created by Srinivas Gorur-Shandilya at 1:45 , 04 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

% load the data

load('/local-data/DA-paper/natural-flickering/without-lfp/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);
tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

% remove the baseline from the PID, and remember the error
PID_baseline = mean(mean(PID(1:5e3,:)));
PID = PID - PID_baseline;

%% Fitting Excursions with Hill functions
% In this document, we attempt to quantify the degree and timescale of dynamic gain control by fitting Hill functions to each excursion in the space of ORN response vs. projected stimulus. The rationale is that this is the simplest description of the ORN response that shows changing gain, and that fitting Hill functions should prevent our estimation of dynamic gain timescale being contaminated by effects of static nonlinearities. 

%% Negative Control: LN model
% In this section we generate some synthetic data using a LN model, and back out filters and then fit nonlinearities in each excursion and attempt to find the timescale of dynamic gain control using these nonlinearities. In the following figure, we plot the ORN response. vs the projected stimulus for each excursion (defined as crossing 10% of peak response), and then attempt to fit a Hill function to each excursion (shown in the second plot). The third plot shows the best fit K_D vs. the r2. The last two plots show how K_D varies with the mean stimulus in a history length (500ms) and how the Spearman correlation of this varies with history length. 

%%
% In this particular plot, the K_D barely changes, as we would expect from a LN model. We keep the cooperativity (n) of the Hill function fixed here, and just fit the K_D. 

clear p
p.   tau1 = 25;
p.   tau2 = 60;
p.      A = 0.4000;
p.      n = 2;
p. offset = 0.05;
p. Hill_A = 100;
p.Hill_Kd = 0.3;
p. Hill_n = 4;

[R,K] = pLNModel(nanmean(PID,2),p);
filtertime = (1:length(K))*1e-3;

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

plotHillExcursions(nanmean(PID,2),fp,R,false);
suptitle('-ve control: LN model')
prettyFig('fs=20;','fixLogX=true;');

if being_published
	snapnow
	delete(gcf)
end

%% Positive Control: DA model 
% In this section we generate synthetic data using a DA model. The same figure as described before is generated. Here, we see that the best-fit Hill functions to each excursion change. 

% fit DA model
clear p
p.   s0 = 7.8242e-04;
p.  n_z = 2;
p.tau_z = 151.1249;
p.  n_y = 2;
p.tau_y = 26.7002;
p.    C = 0.5457;
p.    A = 163.2252;
p.    B = 2.4703;
DA = DAModelv2(nanmean(PID,2),p);

% make a linear filter
R = mean(DA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

plotHillExcursions(nanmean(PID,2),fp,R,true);
suptitle('+ve control: DA model')
prettyFig('fs=20;','fixLogX=true;');

if being_published
	snapnow
	delete(gcf)
end

%%
% While the K_D does change with the mean stimulus in the preceding 500ms, the Spearman correlation does not go to zero as history lengths go to zero, which is weird. The peak of the Spearman correlation vs. history length plot is also at the wrong location (the nominal gain timescale in the model is ~300ms)

%%
% What if we allow n to also vary as we fit these excursions? Will that make the Spearman correlation go to zero at very small history lengths? 

plotHillExcursions(nanmean(PID,2),fp,R,false);
suptitle('+ve control: DA model. K_D and N are allowed to vary')
prettyFig('fs=20;','fixLogX=true;');

if being_published
	snapnow
	delete(gcf)
end

%% Actual Data
% Now we perform the analysis on real data, fitting Hill functions keeping the n the same. 


% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);


% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);


plotHillExcursions(nanmean(PID,2),fp,R,true);
suptitle('ORN Data')
prettyFig('fs=20;','fixLogX=true;');


if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;


