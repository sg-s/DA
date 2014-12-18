% 
% 
% created by Srinivas Gorur-Shandilya at 8:38 , 12 November 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
% redo = 1; deliberately unset

data_root = '/local-data/DA-paper/fast-flicker/orn/';
allfiles  = dir(strcat(data_root,'*.mat'));

% combine all data
% combined_data = ReduceORNData(data_root,allfiles);
% paradigm_names = unique(combined_data.paradigm);
% save('MeanShiftedGaussians.mat','combined_data','paradigm_names');




% load cached data
load('MeanShiftedGaussians.mat')


%% Mean Shifted Gaussians
% How to ORNs respond to mean shifted gaussians? Specifcally, how do they vary their input-output curve? Is the adaptation to this mean optimal (a la Laughlin etc)? Can we find evidence for fast gain adaptation in the curves themselves? 

%         ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%        ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%        ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%         ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%              ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%        ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%         ######     ##    #### ##     ##  #######  ########  #######   ######  

%% Stimulus Distributions 
% Three gaussians with different means are chosen. The following figure shows the distributions for every trial of the data analysed here, colour coded by the experimental paradigm. 

% some global parameters
nbins = 50;
histx = [];
histy = [];
paradigm = [];
dt = 3e-3;
all_pid = [];

a = floor(30/3e-3);
z = floor(55/3e-3);

% assemble all histograms
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	for j = 1:length(plot_these)
		this_pid = combined_data.PID(plot_these(j),:);
		this_pid = this_pid(a:z);
		[y,x] = hist(this_pid,nbins);
		if mean(x) > .5
			histx = [histx; x];
			histy = [histy; y];
			paradigm = [paradigm i];
			all_pid = [all_pid; this_pid];
		end
	end
end

					



c = parula(length(paradigm_names));
paradigm = paradigm -  min(paradigm);
paradigm = paradigm  + 1;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
opacity = .2;
for i = 1:length(unique(paradigm))
	x = histx(paradigm==i,:);
	y = histy(paradigm==i,:);

	for j = 1:width(x)
		plot(x(j,:),y(j,:),'Color',[c(i,:) opacity])
	end
	plot(mean(x),mean(y),'Color',[c(i,:)])
end

xlabel('PID (V)')
ylabel('Count')
titlestr = strcat(mat2str(length(paradigm)),' trials');
title(titlestr)


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%% Stimulus Reproducibility 
% In this section, we look at the reproducibility of the stimulus. The following figure shows the stimulus for all the data we look at here, plotted on top of each other, colour-coded by experimental paradigm. 

c = parula(length(paradigm_names));
time = 1:length(all_pid);
time = time*dt;

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on

opacity = .2;
for i = 1:length(unique(paradigm))
	x = all_pid(paradigm==i,:);
	y = histy(paradigm==i,:);

	for j = 1:width(x)
		plot(time,x(j,:),'Color',[c(i,:) opacity]) % using undocumented opacity value, may break in future
	end
	plot(time,mean(x),'Color',[c(i,:)])
end

ylabel('PID (V)')
xlabel('Time (s)')
titlestr = strcat(mat2str(length(paradigm)),' trials');
title(titlestr)

set(gca,'XLim',[15 20])

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%          ########  ########  ######  ########   #######  ##    ##  ######  ######## 
%          ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
%          ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
%          ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
%          ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
%          ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
%          ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 

%% Neuron Responses: Overview
% The following figure shows the responses of the ORNs to this stimuli, and their distributions. 

c = parula(length(paradigm_names));
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,5,1:4), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_orn));
	plot(time,this_orn,'Color',c(i,:))
end
ylim = get(gca,'YLim');
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

subplot(1,5,5), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_orn = this_orn(a:z);
	[x,y] = hist(this_orn,50);
	plot(x,y,'Color',c(i,:))

end
set(gca,'YLim',ylim);
xlabel('Count')
title('Distribution')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%      ######   #######  ########  ########     ######## #### ##     ## ########  ######  
%     ##    ## ##     ## ##     ## ##     ##       ##     ##  ###   ### ##       ##    ## 
%     ##       ##     ## ##     ## ##     ##       ##     ##  #### #### ##       ##       
%     ##       ##     ## ########  ########        ##     ##  ## ### ## ######    ######  
%     ##       ##     ## ##   ##   ##   ##         ##     ##  ##     ## ##             ## 
%     ##    ## ##     ## ##    ##  ##    ##        ##     ##  ##     ## ##       ##    ## 
%      ######   #######  ##     ## ##     ##       ##    #### ##     ## ########  ######  

%% Correlation Times 
% On what timescales are the stimulus and the response correlated? In the following figure, we plot the autocorrelation function of the stimulus and the responses, for each stimulus presented:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

c = parula(length(paradigm_names));

subplot(1,2,1), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	this_pid= this_pid(a:z);
	[y,x]=autocorr(this_pid,500);
	x = x*3e-3;
  	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','XMinorTick','on')
xlabel('Time (s)')
ylabel('Autocorrelation')
title('Stimulus')


subplot(1,2,2), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_orn = this_orn(a:z);
	[y,x]=autocorr(this_orn,500);
	x = x*3e-3;
  	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','XMinorTick','on')
xlabel('Time (s)')
ylabel('Autocorrelation')
title('ORN Responses')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%          ######## ########  ######## ##    ## ########   ######  
%             ##    ##     ## ##       ###   ## ##     ## ##    ## 
%             ##    ##     ## ##       ####  ## ##     ## ##       
%             ##    ########  ######   ## ## ## ##     ##  ######  
%             ##    ##   ##   ##       ##  #### ##     ##       ## 
%             ##    ##    ##  ##       ##   ### ##     ## ##    ## 
%             ##    ##     ## ######## ##    ## ########   ######  


%% Trends in Data
% Are there any trends in the data? In the following figure, we coarse-grain the data by binning everything along 5-second bins to look at long-term trends in the data. The various colors correspond to various stimulus means, and correspond to other figures in this document. 

% plot_data is indexed by where we start
all_start = [15:5:50];
all_end = all_start+5;


for i = 1:length(paradigm_names)
	plot_data(i).stim_slope = [];
	plot_data(i).stim_slope_err = [];
	plot_data(i).stim_mean = [];
	plot_data(i).stim_mean_err = [];
	plot_data(i).resp_slope = [];
	plot_data(i).resp_slope_err = [];
	plot_data(i).resp_mean = [];
	plot_data(i).resp_mean_err = [];

	for j = 1:length(all_start)
		a = floor(all_start(j)/3e-3);
		z = floor(all_end(j)/3e-3);
		n = sqrt(z-a);

		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		these_pid=mean2(combined_data.PID(plot_these,:));
		these_resp=mean2(combined_data.fA(:,plot_these));

		cropped_pid = these_pid(a:z);
		cropped_resp = these_resp(a:z);

		plot_data(i).stim_mean = 		[plot_data(i).stim_mean mean(cropped_pid)];
		plot_data(i).stim_mean_err = 	[plot_data(i).stim_mean_err std(cropped_pid)/n];

		plot_data(i).resp_mean = 		[plot_data(i).resp_mean mean(cropped_resp)];
		plot_data(i).resp_mean_err = 	[plot_data(i).resp_mean_err std(cropped_resp)/n];


	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(plot_data)
	errorbar(all_start+2.5,plot_data(i).stim_mean,plot_data(i).stim_mean_err,'Color',c(i,:))

end
xlabel('Time (s)')
ylabel('PID (V)')

subplot(1,2,2), hold on
c = parula(length(paradigm_names));
for i = 1:length(plot_data)
	errorbar(all_start+2.5,plot_data(i).resp_mean,plot_data(i).resp_mean_err,'Color',c(i,:))

end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% OK, there are clearly trends in some of the data, especially in the responses. Now, we will attempt to remove all the trends by fitting a second-degree polynomial to the stimulus and response from 35 to 55 seconds. 


% remove trend
b = floor(5/3e-3);
a = floor(35/3e-3);
z = floor(55/3e-3);
clear detrended_data
detrended_data.time = [];
detrended_data.stim = [];
detrended_data.resp = [];
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_resp));
	baseline = mean(this_pid(1:b));
	time = time(a:z);
	this_pid = this_pid(a:z);
	this_resp = this_resp(a:z);

	detrended_data(i).time = time;
	ff = fit(time(:),this_pid(:),'poly2');
	detrended_data(i).stim = this_pid - ff(time)' + mean(ff(time)) - baseline;

	ff = fit(time(:),this_resp(:),'poly2');
	detrended_data(i).resp = this_resp - ff(time) + mean(ff(time));
end

%     ########   #######  ##      ## ######## ########      ######  ########  ########  ######  
%     ##     ## ##     ## ##  ##  ## ##       ##     ##    ##    ## ##     ## ##       ##    ## 
%     ##     ## ##     ## ##  ##  ## ##       ##     ##    ##       ##     ## ##       ##       
%     ########  ##     ## ##  ##  ## ######   ########      ######  ########  ######   ##       
%     ##        ##     ## ##  ##  ## ##       ##   ##            ## ##        ##       ##       
%     ##        ##     ## ##  ##  ## ##       ##    ##     ##    ## ##        ##       ##    ## 
%     ##         #######   ###  ###  ######## ##     ##     ######  ##        ########  ######  


%% Stimulus and Response: Power Spectral Density
% In this section, we compute the power spectral density of the stimulus and the response for each trace, during the times when the stimulus is flickering. The blue traces are the stimulus, and the red is the response. 


all_start = [1:floor(length(detrended_data(i).time)/20):(length(detrended_data(i).time))];
all_start(end) = [];
all_end = all_start + floor(length(detrended_data(i).time)/20);

figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
for i = 1:8
	subplot(2,4,i), hold on
	x = zeros(1,257);
	y = zeros(1,257);
	for j = 1:length(all_start)
		a = floor(all_start(j));
		z = floor(all_end(j));
		n = sqrt(z-a);
		[x,y2]=powerspec(1/3e-3,detrended_data(i).stim(a:z)/max(detrended_data(i).stim(a:z)),0);
		y = y+y2;
	end
	y = y/length(all_start);
	plot(x,y,'b')
	x = zeros(1,257);
	y = zeros(1,257);
	for j = 1:length(all_start)
		a = floor(all_start(j));
		z = floor(all_end(j));
		n = sqrt(z-a);
		[x,y2]=powerspec(1/3e-3,detrended_data(i).resp(a:z)/max(detrended_data(i).resp(a:z)),0);
		y = y+y2';
	end
	y = y/length(all_start);
	plot(x,y,'r')
	title(paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end))
	xlabel('Frequency (Hz)')
	ylabel('Power (norm)')
	set(gca,'XScale','log','YScale','log')
end 
PrettyFig;


if being_published
	snapnow
	delete(gcf)
end

%         ######      ###    #### ##    ##       ##  ######  ########  ######## ######## ########  
%        ##    ##    ## ##    ##  ###   ##      ##  ##    ## ##     ## ##       ##       ##     ## 
%        ##         ##   ##   ##  ####  ##     ##   ##       ##     ## ##       ##       ##     ## 
%        ##   #### ##     ##  ##  ## ## ##    ##     ######  ########  ######   ######   ##     ## 
%        ##    ##  #########  ##  ##  ####   ##           ## ##        ##       ##       ##     ## 
%        ##    ##  ##     ##  ##  ##   ###  ##      ##    ## ##        ##       ##       ##     ## 
%         ######   ##     ## #### ##    ## ##        ######  ##        ######## ######## ########  

%% Neuron Responses: Gain-Speed Tradeoff
% It is believed that when the stimulus is low, the gain needs to be high, which means the ORN integrates the signal over a longer time, which in turn means that the response kinetics are slower. To check if this is the case in this dataset, we compute the cross-correlation functions between the stimulus and the response in each case and see if the width depends in an obvious way on the mean stimulus. 

c = parula(length(detrended_data));
figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
subplot(1,3,1:2), hold on
peak_loc = zeros(length(detrended_data),1);
mean_stim = zeros(length(detrended_data),1);
for i = 1:length(detrended_data)
	mean_stim(i) = mean(detrended_data(i).stim);
	a = detrended_data(i).resp - mean(detrended_data(i).resp);
	a = a/std(a);
	b = detrended_data(i).stim - mean(detrended_data(i).stim);
	b = b/std(b);
	x = xcorr(a,b); % positive peak means a lags b
	t = 3e-3*(1:length(x));
	t = t-mean(t);
	x = x/max(x);
	plot(t,x,'Color',c(i,:))
	[~,loc] = max(x);
	peak_loc(i) = t(loc);
end
set(gca,'XLim',[-.2 .5])
xlabel('Lag (s)')
ylabel('Cross Correlation (norm)')
L = paradigm_names;
for i = 1:length(L)
	L{i} = L{i}(strfind(L{i},'-')+1:end);
end
legend(L)
subplot(1,3,3), hold on
plot(mean_stim,peak_loc*1000,'+-k')
xlabel('Mean Stimulus (V)')
ylabel('Peak of cross correlation (ms)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%      ####       ##  #######      ######  ##     ## ########  ##     ## ########  ######  
%       ##       ##  ##     ##    ##    ## ##     ## ##     ## ##     ## ##       ##    ## 
%       ##      ##   ##     ##    ##       ##     ## ##     ## ##     ## ##       ##       
%       ##     ##    ##     ##    ##       ##     ## ########  ##     ## ######    ######  
%       ##    ##     ##     ##    ##       ##     ## ##   ##    ##   ##  ##             ## 
%       ##   ##      ##     ##    ##    ## ##     ## ##    ##    ## ##   ##       ##    ## 
%      #### ##        #######      ######   #######  ##     ##    ###    ########  ######  

%% Neuron Responses: Input-Output Curve Changes
% How do ORNs change their input-output curves when presented with these different stimuli? If the ORN is a perfect sensor, that adapts perfectly, it would move and stretch its I/O curve so that it is equal to cumulative density function of the stimulus. However, we know that ORNs don't adapt perfectly, so they must do something else. What is it? If the change dominated by a lateral movement of a stretch?

% this stores the LN model parameters for all the data
if redo
	clear LNModel
	LNModel.p = []; % parameterised LN Model
	LNModel.LinearFit = [];
	LNModel.LNFit = [];


	for i = 1:length(detrended_data)
		this_orn = detrended_data(i).resp;
		this_pid = detrended_data(i).stim;
		time     = detrended_data(i).time;

		clear p
  		p.   tau1= 4.5591;
  		p.    K_n= 4.9181;
  		p.   tau2= 19.6875;
  		p.      A= max(this_orn);
  		p.      n= 2.4297;
  		p.     Kd= mean(this_orn);
  		p. offset= mean(this_orn);
  		p.    K_A= 0.3354;

  		% fit the model
  		clear data
  		data.stimulus = this_pid;
  		data.response = this_orn;
  		data.time = time;
  		LNModel(i).p=FitModel2Data(@pLNModel,data,p);

  		% solve for best fit parameters
  		[LNModel(i).LNFit,LNModel(i).K,LNModel(i).LinearFit] = pLNModel(this_pid,LNModel(i).p);


		% save all of this for later
		LNModel(i).this_orn = this_orn;
		LNModel(i).LNFit_r2 = rsquare(LNModel(i).LNFit,this_orn);
	end

	% calculate filters as in Baccus and Meister (zero mean, unit variance i/o)
	clear BMModel
	BMModel.K = [];
	BMModel.LinearFit = [];


	for i = 1:length(detrended_data)
		this_orn = detrended_data(i).resp;
		this_pid = detrended_data(i).stim;
		time     = detrended_data(i).time;

		this_orn = this_orn - mean(this_orn);
		this_orn = this_orn/std(this_orn);

		this_pid = this_pid - mean(this_pid);
		this_pid = this_pid/std(this_pid);

		[K,~,filtertime] = FindBestFilter(this_pid,this_orn,[],'filter_length=299;');
		filtertime = filtertime*mean(diff(time));
		BMModel(i).LinearFit = convolve(time,this_pid,K,filtertime);
		BMModel(i).K = K;


	end

end

%%
% The following figure shows each of the ORN responses to each stimulus together with the best LN Model fits:


figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
for i = 1:8
	subplot(2,4,i), hold on
	time = detrended_data(i).time;
	plot(time,LNModel(i).this_orn,'k');
	aa = min(time);
	time = 3e-3*(1:length(LNModel(i).LNFit)) + aa;
	lh=plot(time,LNModel(i).LNFit,'r');
	set(gca,'XLim',[40 50])
	legend(lh,oval(LNModel(i).LNFit_r2,2))
	title(paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end))
	xlabel('Time (s)')
	ylabel('Firing Rate (Hz)')
end 
PrettyFig;


if being_published
	snapnow
	delete(gcf)
end


%%
% The following figure shows the filters backed out of each data set (corresponding to each mean), and the corresponding non-linearities.  

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(LNModel));
for i = 1:length(LNModel)
	plot(3e-3*(1:length(LNModel(i).K)),LNModel(i).K,'Color',c(i,:))
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (norm)')

subplot(1,2,2), hold on
for i = 1:length(LNModel)
	x = sort(LNModel(i).LinearFit);
	plot(x,hill([LNModel(i).p.A LNModel(i).p.Kd LNModel(i).p.n],x),'Color',c(i,:))
end
xlabel('Filter Output (a.u.)')
ylabel('ORN Firing Rate (Hz)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% Why does the LN model do so poorly in higher stimulus conditions? Two possibilities are: 1) there is something interesting going on in higher doses, making the LN model fit less informative or 2) because the deviations in the signal are smaller in higher doses, the poor fit arises from counting noise. To figure out what's going on, we take the filter from the highest concentration case and use it to predict the response to the lowest concentration:
 
clear p ph L
p=LNModel(8).p;
p.A = 41.3682;
p.offset = 20.0278;
p.n = 4.4;
p.Kd = 20;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(detrended_data(1).time,detrended_data(1).resp,'k');
ph(1)=plot(detrended_data(1).time,LNModel(1).LNFit,'r');
ph(2)=plot(detrended_data(1).time,pLNModel(detrended_data(1).stim,p),'g');
r2 = rsquare(pLNModel(detrended_data(1).stim,p),detrended_data(1).resp);
L = {};
L{1} = oval(LNModel(1).LNFit_r2,2);
L{2} = oval(r2,2);
legend(ph,L)

set(gca,'XLim',[40 50])

xlabel('Time (s)')
ylabel('Firing Rate (Hz)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% The fit is not totally horrible, and is much better than the fit to the responses to higher concentrations. 

%%
% There seems to be a systematic variation of the filter shape with odor dose. In the following plot, we plot the location of the maximum and the minimum of the filters as a function of stimulus. 

[~,max_loc]=max(reshape([LNModel.K],300,8));
[~,min_loc]=min(reshape([LNModel.K],300,8));
ms = mean(reshape([detrended_data.stim],length(detrended_data(1).stim),8));


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(ms,max_loc*3e-3,'r-+')
plot(ms,min_loc*3e-3,'k-+')
legend({'Maximum','Minimum'})

set(gca,'XScale','log')
ylabel('Filter Peaks')
xlabel('Mean Stimulus')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% We see a similar plot to the auto-correlation analysis. Except for the first point, everything seems to be getting faster with increasing stimulus. 


 % ######      ###    #### ##    ##     ######  ##     ##    ###    ##    ##  ######   ########  ######  
% ##    ##    ## ##    ##  ###   ##    ##    ## ##     ##   ## ##   ###   ## ##    ##  ##       ##    ## 
% ##         ##   ##   ##  ####  ##    ##       ##     ##  ##   ##  ####  ## ##        ##       ##       
% ##   #### ##     ##  ##  ## ## ##    ##       ######### ##     ## ## ## ## ##   #### ######    ######  
% ##    ##  #########  ##  ##  ####    ##       ##     ## ######### ##  #### ##    ##  ##             ## 
% ##    ##  ##     ##  ##  ##   ###    ##    ## ##     ## ##     ## ##   ### ##    ##  ##       ##    ## 
 % ######   ##     ## #### ##    ##     ######  ##     ## ##     ## ##    ##  ######   ########  ######  

%%
% Now, we study the gain changes in this system by plotting the filtered stimulus vs. the actual response. The stimulus has been means subtracted and are divided through by their standard deviation. 

nbins = 10;
plot_range = 5;
figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:8
	subplot(2,4,i), hold on
	this_orn = detrended_data(i).resp;
	this_pid = detrended_data(i).stim;
	this_pid = this_pid - mean(this_pid);
	this_orn = this_orn - mean(this_orn);
	time = 3e-3*(1:length(this_pid));

	% recompute linear fits
	BMModel(i).LinearFit = convolve(time,this_pid,K,filtertime);

	%scatter(BMModel(i).LinearFit,this_orn,32,[0.8 0.8 0.8],'filled')
	plot(BMModel(i).LinearFit,this_orn,'.','Color',[.8 .8 .8]) 

	% bin the data
	minx = min(BMModel(i).LinearFit);
	maxx = max(BMModel(i).LinearFit);
	r = (maxx -minx)/nbins;
	a = minx;	

	BMModel(i).y = [];
	BMModel(i).ye = [];
	BMModel(i).x = linspace(minx,maxx,10);

	for j = 1:nbins
		z = a+r;
		BMModel(i).y = [BMModel(i).y mean(this_orn(BMModel(i).LinearFit < z & BMModel(i).LinearFit > a))];
		BMModel(i).ye = [BMModel(i).ye std(this_orn(BMModel(i).LinearFit < z & BMModel(i).LinearFit > a))];
		a = a+r;
	end

	% also fit a line to the middle 1/3 of the data
	a = minx + r/3;
	z = maxx - r/3;
	x = BMModel(i).LinearFit(BMModel(i).LinearFit < z & BMModel(i).LinearFit > a);
	y = this_orn(BMModel(i).LinearFit < z & BMModel(i).LinearFit > a);
	ff = fit(x,y,'poly1');
	temp=confint(ff);
	BMModel(i).gain_err = temp(2,1)-temp(1,1);
	BMModel(i).gain = ff.p1;

	errorbar(BMModel(i).x,BMModel(i).y,BMModel(i).ye,'k')
	xlabel('K\otimess')
	ylabel('f-mean(f) (Hz)')
	title(paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end))
	set(gca,'XLim',[-.7 .7],'YLim',[-20 20])
end

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% The following plot shows the gain changes in the different cases:

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(mean(vertcat(detrended_data.stim)'),[BMModel.gain],[BMModel.gain_err])
xlabel('Stimulus Mean (V)')
ylabel('Gain (Hz/V)')
set(gca,'XScale','log','YScale','log')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% The following figure shows the co-efficient of variation of the stimulus and the response. 


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(cv(vertcat(detrended_data.stim)),cv(horzcat(detrended_data.resp)),'+-k')
xlabel('CV (stimulus)')
ylabel('CV (response)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%      ########    ###     ######  ########       ###    ########     ###    ########  ######## 
%      ##         ## ##   ##    ##    ##         ## ##   ##     ##   ## ##   ##     ##    ##    
%      ##        ##   ##  ##          ##        ##   ##  ##     ##  ##   ##  ##     ##    ##    
%      ######   ##     ##  ######     ##       ##     ## ##     ## ##     ## ########     ##    
%      ##       #########       ##    ##       ######### ##     ## ######### ##           ##    
%      ##       ##     ## ##    ##    ##       ##     ## ##     ## ##     ## ##           ##    
%      ##       ##     ##  ######     ##       ##     ## ########  ##     ## ##           ##    


%% Neuron Responses: Fast Adaptation 
% How well do LN models predict the response of neurons in these paradigms? Do we see evidence of fast gain adaptation in these data sets? Now we perform our gain analysis of the prediction using standard methods described elsewhere in this repo. We do the analysis for each stimulus set. 


for i = 1 % only doing first because there is a trend in all others
	ph = [];

	history_lengths = (3*floor(1000*logspace(-1.5,1,30)/3))/1e3;
	example_history_length = 0.135;

	f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	ph(3) = subplot(1,2,1); hold on 
	axis square
	ph(4) = subplot(1,2,2); hold on


	clear x
	x.response = LNModel(i).this_orn(1:end-32); % the 32 is to account for the acausal part of the filter
	x.prediction = LNModel(i).LNFit(1:end-32);
	x.stimulus = detrended_data(i).stim(1:end-32); 
	x.time = detrended_data(i).time(1:end-32);
	x.filter_length = 299;

	if redo
		[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
		s=abs(l-h);
		s(p_LN(1,:)>0.05)=NaN;
		[~,loc]=max(s);

		% save it for later
		LNModel(i).ehl = history_lengths(loc);
		LNModel(i).pb = p_LN;

	else
		GainAnalysis4(x,history_lengths,LNModel(i).ehl,ph,LNModel(i).pb);
	end

	xlabel(ph(3),'LN Prediction (Hz)')
	set(ph(4),'XScale','log')
	title(ph(4),paradigm_names{i})

	if being_published

		snapnow;
		delete(f2);
	end
end

%         ########     ###       ##     ##  #######  ########  ######## ##       
%         ##     ##   ## ##      ###   ### ##     ## ##     ## ##       ##       
%         ##     ##  ##   ##     #### #### ##     ## ##     ## ##       ##       
%         ##     ## ##     ##    ## ### ## ##     ## ##     ## ######   ##       
%         ##     ## #########    ##     ## ##     ## ##     ## ##       ##       
%         ##     ## ##     ##    ##     ## ##     ## ##     ## ##       ##       
%         ########  ##     ##    ##     ##  #######  ########  ######## ######## 


%% The DA Model
% In this section, we fit the DA model to the data. First, we fit a DA model to the first trace, which corresponds to the lowest stimulus and also the only stimulus where there is no large trend in the data. 

i=1;
plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
this_pid=mean2(combined_data.PID(plot_these,:));
this_resp=mean2(combined_data.fA(:,plot_these));

% remove baseline
this_pid = this_pid-mean(this_pid(400:1600));

clear data
data.response = this_resp;
data.time = 3e-3*(1:length(this_resp));
data.stimulus = this_pid;

%[p, Rguess,x ] = FitDAModelToData(data);

x1 = [162.9375    5.5391    0.0012    0.8672   25.9746  188.0098    0.5977 -0.0013];
p = ValidateDAParameters2(x1);
R = DA_integrate2(data.stimulus,p);

figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
plot(data.time,data.response,'k')
lh=plot(data.time,R,'r');
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
r = rsquare(data.response,R);
legend(lh,strcat('r^{2}=',oval(r,2)))
title('Best Fit DA Model Prediction (3%)')
PrettyFig;

if being_published

	snapnow;
	delete(gcf);
end

%%
% The parameters of this DA model are:
disp(p)


%%
% In the following figure, we attempt to fit a DA model the case where the mean of the stimulus is highest. For this, we do two different things: first, we fit a DA model to the data directly, as before. Second, we use DA model parameters from the first fit (to the low stimulus) and use it to predict the response to the high fit. 

%%
% In the following figure, the black curve is the response to the highest concentration, the green curve is the DA model prediction from the response to the *lowest* concentration, and the red curve is the prediction from the best-fit DA Model. 

i=6;
plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
this_pid=mean2(combined_data.PID(plot_these,:));
this_resp=mean2(combined_data.fA(:,plot_these));
% remove baseline
this_pid = this_pid-mean(this_pid(400:1600));

clear data
data.response = this_resp;
data.time = 3e-3*(1:length(this_resp));
data.stimulus = this_pid;

x6 = [  187.9277    5.4844    0.0012    8.8340    2.0156  186.6661    0.0156  -0.0013];
p = ValidateDAParameters2(x6);
R = DA_integrate2(data.stimulus,p);

figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
plot(data.time,data.response,'k')
clear lh
p = ValidateDAParameters2(x1);
R = DA_integrate2(data.stimulus,p);
lh(1)=plot(data.time,R,'g');
r1 = rsquare(data.response(300:end),R(300:end));
p = ValidateDAParameters2(x6);
R = DA_integrate2(data.stimulus,p);
lh(2)=plot(data.time,R,'r');
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
r6 = rsquare(data.response(300:end),R(300:end));
legend(lh,{ strcat('r^{2}=',oval(r1,2)),strcat('r^{2}=',oval(r6,2))})
title('Best Fit DA Model Prediction (6.25%)')
PrettyFig;

if being_published

	snapnow;
	delete(gcf);
end

%%
% Emboldened by this result, we now use the DA model fit to the lowest stimulus to predict responses to all the other stimuli, and compute the gain for each case, and compare to the actual data: 

x1 = [162.9375    5.5391    0.0012    0.8672   25.9746  188.0098    0.5977 -0.0013];
p = ValidateDAParameters2(x1);

nbins = 10;
gain = [];
gain_err = [];

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	this_pid = this_pid-mean(this_pid(400:1600));
	R = DA_integrate2(this_pid,p);
	a = floor(20/3e-3);
	z = floor(55/3e-3);
	R = R(a:z); this_pid = this_pid(a:z);
	time = 3e-3*(1:length(R));

	R = R-mean(R);
	this_pid = this_pid -mean(this_pid);

	[K,~,filtertime] = FindBestFilter(this_pid/std(this_pid),R/std(R),[],'filter_length=299;');
	filtertime  = filtertime*3e-3;

	% recompute linear fits
	LinearFit = convolve(time,this_pid,K,filtertime);

	% bin the data
	minx = min(LinearFit);
	maxx = max(LinearFit);
	r = (maxx -minx)/nbins;


	% also fit a line to the middle 1/3 of the data
	a = minx + r/3;
	z = maxx - r/3;
	x = LinearFit(LinearFit < z & LinearFit > a);
	y = R(LinearFit < z & LinearFit > a);
	ff = fit(x(:),y(:),'poly1');
	temp=confint(ff);
	gain_err = [gain_err temp(2,1)-temp(1,1)];
	gain = [gain ff.p1];

end

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(mean(vertcat(detrended_data.stim)'),[BMModel.gain],[BMModel.gain_err],'k')
errorbar(mean(vertcat(detrended_data.stim)'),gain,gain_err,'r')
xlabel('Stimulus Mean (V)')
ylabel('Gain (Hz/V)')
set(gca,'XScale','log')
legend('Data','DA Model fit to lowest dose')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%            ########     ###       ##     ##  #######  ########  ######## ##       
%            ##     ##   ## ##      ###   ### ##     ## ##     ## ##       ##       
%            ##     ##  ##   ##     #### #### ##     ## ##     ## ##       ##       
%            ########  ##     ##    ## ### ## ##     ## ##     ## ######   ##       
%            ##   ##   #########    ##     ## ##     ## ##     ## ##       ##       
%            ##    ##  ##     ##    ##     ## ##     ## ##     ## ##       ##       
%            ##     ## ##     ##    ##     ##  #######  ########  ######## ######## 


%% Fitting a Receptor-Adaptation Model
% In this section, we consider a class of models where the affinity of the receptor changes in a stimulus-driven fashion. Thierry suggested the following functional form:
% 
% $$f=N\left(K_{r}\otimes\frac{s(t)}{s(t)+K_{a}\otimes s(t)}\right)$$
%

%%
% where $K_{r}$ and $K_{a}$ are response and adaptation filters (here parametrised as double gamma functions) and $N$ is a Hill function parametrised by height $A$, offset $k$, and steepness $n$. However, this form cannot work. Consider that $K_{a}$ is simply a delta function. Then, 
%
% $$f\sim K_{r}\otimes\frac{1}{2}$$
%
% which is not what we want. 

%%
% Instead I consider the following functional form:
%
% $$f=N\left(K_{r}\otimes\frac{s(t)}{1+\beta*K_{a}\otimes s(t)}\right)$$
%
% Here, if we let $A,k\rightarrow\infty$ and set $n=1$, the non-linearity drops out. Furthermore, if we let the amplitude of the adaptation filter $K_{a}\rightarrow0$, the "adapting" part of the model drops away and we are left with the simple linear model. 

%%
% Armed with this intuition, we attempt to fit the data (the response to the lowest dose) with this reduced model, which should give us, essentially, a simpler LN fit to the data. The following figure shows such a fit: where we neglect the adaptation filter.

a = 5000;
z = 18000;

plot_these=find(strcmp(paradigm_names{1}, combined_data.paradigm));
this_orn=mean2(combined_data.fA(:,plot_these));
this_pid=mean2(combined_data.PID(plot_these,:));
time = 3e-3*(1:length(this_orn));
data.response = this_orn(a:z);
data.stimulus = this_pid(a:z);
data.time = time(a:z);


clear p
p.   tau1= 5.9218;
p.    K_n= 3.6250;
p.   tau2= 27.5625;
p.    K_A= 0.6094;
p. a_tau1= 10;
p.    a_n= 1.3750;
p. a_tau2= 51.3750;
p.    a_A= -1;
p.beta= 0;
p.      A= 51.0312;
p.      n= 2;
p.     Kd= 14;
p. offset= 11.9479;

%				 Kr-----------------|-----Ka ---------------   |----Hill
lb = mat2struct([1   1  2  1e-3 		1 	 1 2    -2      -1    1      1   1      0],fieldnames(p));
ub = mat2struct([100 5 200 1e3 			100  5 200  -1      0     1e3    5   1e3    1e2],fieldnames(p));

%p=FitModel2Data(@ReceptorAdaptationModel,data,p,lb,ub);


fp=ReceptorAdaptationModel(data.stimulus,p);

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on

lh = [];
L = {};
subplot(2,2,1:2), hold on
plot(data.time,data.response,'k')
lh=plot(data.time,fp,'r');
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend(lh,strcat('r^2=',oval(rsquare(fp,data.response),2)));
set(gca,'XLim',[35 55])

subplot(2,2,3)
t=1:300;
Kr = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
t = t*mean(diff(data.time));
plot(t,Kr,'r')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (a.u.)')

subplot(2,2,4)
x = [p.A p.Kd p.n];
f = hill(x,sort(shat2));
plot(sort(shat2),f,'r')
ylabel('Output non-linearity (Hz)')
xlabel('Filter Output (a.u.)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% OK, that looks good. Now, we remove all constraints from the model and allow the adapting filter to play a  role too. 

clear p
p.   tau1= 5.9218;
p.    K_n= 3.6250;
p.   tau2= 27.5625;
p.    K_A= 0.6094;
p. a_tau1= 10;
p.    a_n= 1.3750;
p. a_tau2= 51.3750;
p.    a_A= -1;
p.beta= 1e-4;
p.      A= 51.0312;
p.      n= 2;
p.     Kd= 14;
p. offset= 11.9479;

%				 Kr-----------------|-----Ka ---------------   |----Hill
lb = mat2struct([1   1  2  0 		1 	 1 2    0        -1        1      1   1      0],fieldnames(p));
ub = mat2struct([100 5 200 1        100  5 200  1e-6     1        1e3    5   1e3    1e2],fieldnames(p));
%p=FitModel2Data(@ReceptorAdaptationModel,data,p,lb,ub);

fp=ReceptorAdaptationModel(data.stimulus,p);

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
lh = [];
L = {};
subplot(2,2,1:2), hold on
plot(data.time,data.response,'k')
lh=plot(data.time,fp,'r');
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend(lh,strcat('r^2=',oval(rsquare(fp,data.response),2)));
set(gca,'XLim',[35 55])

subplot(2,2,3), hold on
t=1:300;
Kr = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
Ka = p.beta*filter_gamma2(p.a_tau1,p.a_n,p.a_tau2,p.a_A,t);
t = t*mean(diff(data.time));
plot(t,Kr,'r')
plot(t,Ka,'b')
legend('Response','Adaptation')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (a.u.)')

subplot(2,2,4)
x = [p.A p.Kd p.n];
f = hill(x,sort(shat2));
plot(sort(shat2),f)

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end
