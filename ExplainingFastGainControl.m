% ExplainingFastGainControl.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Explaining Fast Gain Control
% We have shown that there is a lot of evidence that ORNs modulate gain as a function of recently experienced stimuli. We show that the LN model, despite explaining >95% of the data (for Gaussian/binary inputs), fails to account for this fast gain control. Can we develop a model to explain this fast gain control?

%%
% In the following sections, we fit a variety of models (see da-pdfs/list-of-models.pdf) to the data, to see which model explains the data best, and then analyse how well these models account for observed fast gain changes. 


% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end


% load data form the large variance flickering experiment
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);
B_spikes = vertcat(B_spikes,spikes(4).B);

% A spikes --> firing rate
hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

% assemble the data for fitting 
clear d
d.stimulus = mean2(PID);
d.response = mean2(fA);


%           ##       ##    ##    ##     ##  #######  ########  ######## ##       
%           ##       ###   ##    ###   ### ##     ## ##     ## ##       ##       
%           ##       ####  ##    #### #### ##     ## ##     ## ##       ##       
%           ##       ## ## ##    ## ### ## ##     ## ##     ## ######   ##       
%           ##       ##  ####    ##     ## ##     ## ##     ## ##       ##       
%           ##       ##   ###    ##     ## ##     ## ##     ## ##       ##       
%           ######## ##    ##    ##     ##  #######  ########  ######## ######## 

%% Fitting a LN Model
% Here, we fit a non-parametric LN model to the data. 


if ~exist('K')
	[K, ~, filtertime] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=10;','regmin=.1;','filter_length=999;');
	filtertime = filtertime*mean(diff(tA));
	K = K/max(K);
end
fp  =convolve(tA,mean2(PID),K,filtertime);
fp = fp + 18.6131;
fp = fp*1.2246;
clear p
p.A= 56.9845;
p.k= 23.5616;
p.n= 2.9577;

fp_LN = hill(p,fp);
fp_K = fp; clear fp

figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on

subplot(2,4,1:3), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp_K,'r');
r2 = rsquare(fp_K,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
title('Linear Prediction')
ylabel('Firing Rate (Hz)')

subplot(2,4,4), hold on
plot(filtertime,K,'r')
xlabel('Lag (s)')
ylabel('Filter Amplitude')

subplot(2,4,5:7), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp_LN,'r');
r2 = rsquare(fp_LN,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
title('LN Prediction')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

subplot(2,4,8), hold on
plot(sort(fp_K),hill(p,sort(fp_K)),'r')
xlabel('Linear Prediction (Hz)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%           ########     ###       ##     ##  #######  ########  ######## ##       
%           ##     ##   ## ##      ###   ### ##     ## ##     ## ##       ##       
%           ##     ##  ##   ##     #### #### ##     ## ##     ## ##       ##       
%           ##     ## ##     ##    ## ### ## ##     ## ##     ## ######   ##       
%           ##     ## #########    ##     ## ##     ## ##     ## ##       ##       
%           ##     ## ##     ##    ##     ## ##     ## ##     ## ##       ##       
%           ########  ##     ##    ##     ##  #######  ########  ######## ######## 

%% Fitting a DA Model
% In the following section, we fit a DA Model (DAModelv2), where we add a stimulus offset term (because we don't know the offset of the PID), and constrain $n_z$ and $n_y$ to 2. 

clear p
p.   s0= -0.1663;
p.tau_y= 22.9687;
p.  n_y= 2;
p.    C= 0.5701;
p.tau_z= 150.1172;
p.  n_z= 2;
p.    A= 556.6875;
p.    B= 10.9296;


[fp_DA,~,~,Ky,Kz] = DAModelv2(mean2(PID),p);

figure('outerposition',[0 0 1300 500],'PaperUnits','points','PaperSize',[1300 500]); hold on
subplot(1,4,1:3), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp_DA,'r');
r2 = rsquare(fp_LN,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
title('DA Prediction')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

subplot(1,4,4), hold on
plot((0:300)*1e-3,Ky,'r')
plot((0:300)*1e-3,Kz,'b')
legend('K_y','K_z')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% 


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))