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

% shorten paradigm names by throwing out 'MFC'
short_paradigm_names = paradigm_names;
for i = 1:length(paradigm_names)
	short_paradigm_names{i} = paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end);
end


% some global parameters
nbins = 50;
histx = [];
histy = [];
dt = 1e-3;
all_pid = [];


%% Mean Shifted Gaussians
% How to ORNs respond to mean shifted gaussians? Specifcally, how do they vary their input-output curve? Is the adaptation to this mean optimal (a la Laughlin etc)? Can we find evidence for fast gain adaptation in the curves themselves? 

%         ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%        ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%        ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%         ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%              ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%        ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%         ######     ##    #### ##     ##  #######  ########  #######   ######  

%% Stimulus 
% The stimulus we present in this experiment consists of a Gaussian with various means. The following figure shoes the mean of the each Gaussian we present. 


a = floor(15/dt);
z = floor(55/dt);

clear ax
figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
ax(1) = subplot(1,5,1:4); hold on
xlabel('Time (s)')
ylabel('PID (V)')
ax(2) = subplot(1,5,5); hold on
xlabel('p.d.f.')
c = parula(length(paradigm_names));

mean_pid = NaN(length(c),1);

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	plot_hist = (combined_data.PID(plot_these,a:z));
	[~,~,hx,hy]  = splinehist(plot_hist);
	if i ~= 2
		% sometheing weird here causes splinehist to crash
	 	plot(ax(2),hy,hx,'Color',c(i,:))
	end

	plot_this = mean2(combined_data.PID(plot_these,:));
	time = dt*(1:length(plot_this));
	plot(ax(1),time,plot_this,'Color',c(i,:))

	mean_pid(i) = mean(plot_this);
end

PrettyFig;


if being_published
	snapnow
	delete(gcf)
end



%% Stimulus Reproducibility 
% In this section, we look at the reproducibility of the stimulus. The following figure shows the stimulus for all the data we look at here, plotted on top of each other, colour-coded by experimental paradigm. 

fig_handle=figure('Units','pixels','outerposition',[100 302 1400 498],'PaperUnits','points','PaperSize',[1400 498]); hold on
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position', [12.45 12.45 149.4 149.4]);
axes_handles(2)=axes('Units','pixels','Position', [174.3 12.45 149.4 149.4]);
axes_handles(3)=axes('Units','pixels','Position', [336.15 12.45 149.4 149.4]);
axes_handles(4)=axes('Units','pixels','Position', [498 12.45 149.4 149.4]);
axes_handles(5)=axes('Units','pixels','Position', [659.85 12.45 149.4 149.4]);
axes_handles(6)=axes('Units','pixels','Position', [821.7 12.45 149.4 149.4]);
axes_handles(7)=axes('Units','pixels','Position', [983.55 12.45 149.4 149.4]);
axes_handles(8)=axes('Units','pixels','Position', [1145.4 12.45 100+149.4 149.4]);
axes_handles(9)=axes('Units','pixels','Position', [12.45 211.65 149.4 149.4]);
axes_handles(10)=axes('Units','pixels','Position',[174.3 211.65 149.4 149.4]);
axes_handles(11)=axes('Units','pixels','Position',[336.15 211.65 149.4 149.4]);
axes_handles(12)=axes('Units','pixels','Position',[498 211.65 149.4 149.4]);
axes_handles(13)=axes('Units','pixels','Position',[659.85 211.65 149.4 149.4]);
axes_handles(14)=axes('Units','pixels','Position',[821.7 211.65 149.4 149.4]);
axes_handles(15)=axes('Units','pixels','Position',[983.55 211.65 149.4 149.4]);
axes_handles(16)=axes('Units','pixels','Position',[1145.4 211.65 100+149.4 149.4]);

ti = 1;
for i = 1:length(paradigm_names)
	x = combined_data.PID(find(strcmp(paradigm_names{i},combined_data.paradigm)),a:z);
	[r2,s] = rsquare(x);
	axes(axes_handles(i))
	hold(axes_handles(i) ,'on')
	imagescnan(r2)
	caxis([0 1])
	axis image
	if i == length(paradigm_names)
		colorbar
	end
	axis off
	tx = (strcat('r^2 = ',oval(mean(r2(~isnan(r2))),2)));
	text(1,length(r2)*.8,tx)

	axes(axes_handles(i+length(paradigm_names)))
	hold(axes_handles(i+length(paradigm_names)) ,'on')
	imagescnan(s)
	caxis([0 1.4])
	if i == length(paradigm_names)
		colorbar
	end
	axis image
	axis off
	tx = (strcat('slope = ',oval(mean(s(~isnan(s))),2)));
	text(1,length(r2)*.8,tx)

end

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


a = floor(15/dt);
z = floor(55/dt);

clear ax
figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
ax(1) = subplot(1,5,1:4); hold on

xlabel('Time (s)')
ylabel('Firing Rate (Hz) (V)')
ax(2) = subplot(1,5,5); hold on
xlabel('p.d.f.')
c = parula(length(paradigm_names));

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	plot_hist = (combined_data.fA(a:z,plot_these));
	[~,~,hx,hy]  = splinehist(mean2(plot_hist));
	% if i ~= 2
	% 	% sometheing weird here causes splinehist to crash
		plot(ax(2),hy,hx,'Color',c(i,:))
	% end

	plot_this = mean2(combined_data.fA(:,plot_these));
	time = dt*(1:length(plot_this));
	plot(ax(1),time,plot_this,'Color',c(i,:))
end

PrettyFig('EqualiseY=1;');


if being_published
	snapnow
	delete(gcf)
end

%      ########     ###    ########    ###     ######  ######## ######## 
%      ##     ##   ## ##      ##      ## ##   ##    ## ##          ##    
%      ##     ##  ##   ##     ##     ##   ##  ##       ##          ##    
%      ##     ## ##     ##    ##    ##     ##  ######  ######      ##    
%      ##     ## #########    ##    #########       ## ##          ##    
%      ##     ## ##     ##    ##    ##     ## ##    ## ##          ##    
%      ########  ##     ##    ##    ##     ##  ######  ########    ##    

%% Dataset Details
% Because this experiment was performed on different days on different neurons, all direct comparisons are not possible. The following figure shows the number of trials of data we have for each neuron, for each experimental paradigm. 

d = NaN(length(unique(combined_data.neuron)),length(unique(combined_data.paradigm)));
for i = 1:length(unique(combined_data.neuron))
	for j = 1:length(unique(combined_data.paradigm))
		d(i,j) = length(intersect(find(combined_data.neuron == i), find(strcmp(paradigm_names{j},combined_data.paradigm))));
	end
end
d(d==0) =NaN;

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1000 700]); hold on

imagescnan(d)
colorbar
xlabel('Experimental Paradigm')
set(gca,'XTick',[1:length(unique(combined_data.paradigm))],'XTickLabel',short_paradigm_names,'XTickLabelRotation',45)
ylabel('Neuron #')
title('# trials/neuron/experiment')

PrettyFig();


if being_published
	snapnow
	delete(gcf)
end


%% 
% As can be seen, the data exists in three different subsets, corresponding to three different days. 

         
%          ######## ########  ######## ##    ## ########   ######  
%             ##    ##     ## ##       ###   ## ##     ## ##    ## 
%             ##    ##     ## ##       ####  ## ##     ## ##       
%             ##    ########  ######   ## ## ## ##     ##  ######  
%             ##    ##   ##   ##       ##  #### ##     ##       ## 
%             ##    ##    ##  ##       ##   ### ##     ## ##    ## 
%             ##    ##     ## ######## ##    ## ########   ######  


%% Trends in Data
% Are there any trends in the data? In the following figure, we coarse-grain the data by binning everything along 5-second bins to look at long-term trends in the data. The various colors correspond to various stimulus means, and correspond to other figures in this document. 

group = [ones(1,4) 2*ones(1,4) 3*ones(1,5)];

% plot_data is indexed by where we start
all_start = [15:5:50];
all_end = all_start+5;

for k = 1:3 % over groups
	for i = 1:length(paradigm_names)
		plot_data(i,k).stim_slope = [];
		plot_data(i,k).stim_slope_err = [];
		plot_data(i,k).stim_mean = [];
		plot_data(i,k).stim_mean_err = [];
		plot_data(i,k).resp_slope = [];
		plot_data(i,k).resp_slope_err = [];
		plot_data(i,k).resp_mean = [];
		plot_data(i,k).resp_mean_err = [];

		for j = 1:length(all_start)
			a = floor(all_start(j)/dt);
			z = floor(all_end(j)/dt);
			n = sqrt(z-a);

			plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
			plot_these = intersect(plot_these,find(ismember(combined_data.neuron,find(group==k))));
			these_pid=mean2(combined_data.PID(plot_these,:));
			these_resp=mean2(combined_data.fA(:,plot_these));

			cropped_pid = these_pid(a:z);
			cropped_resp = these_resp(a:z);

			plot_data(i,k).stim_mean = 		[plot_data(i,k).stim_mean mean(cropped_pid)];
			plot_data(i,k).stim_mean_err = 	[plot_data(i,k).stim_mean_err std(cropped_pid)/n];
			plot_data(i,k).resp_mean = 		[plot_data(i,k).resp_mean mean(cropped_resp)];
			plot_data(i,k).resp_mean_err = 	[plot_data(i,k).resp_mean_err std(cropped_resp)/n];


		end
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(plot_data)
	for k = 1:3
		errorbar(all_start+2.5,plot_data(i,k).stim_mean,plot_data(i,k).stim_mean_err,'Color',c(i,:))
	end
end
xlabel('Time (s)')
ylabel('PID (V)')

subplot(1,2,2), hold on
c = parula(length(paradigm_names));
for i = 1:length(plot_data)
	for k = 1:3
		errorbar(all_start+2.5,plot_data(i,k).resp_mean,plot_data(i,k).resp_mean_err,'Color',c(i,:))
	end
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

% remove trend -- and combine by experimental group
b = floor(5/dt);
a = floor(35/dt);
z = floor(55/dt);
detrended_data = cache('detrended_data');
if isempty(detrended_data)
	detrended_data.time = [];
	detrended_data.stim = [];
	detrended_data.resp = [];
	for k = 1:3
		for i = 1:length(paradigm_names)
			plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
			plot_these = intersect(plot_these,find(ismember(combined_data.neuron,find(group==k))));
			if ~isempty(plot_these)
				this_pid=mean2(combined_data.PID(plot_these,:));
				this_resp=mean2(combined_data.fA(:,plot_these));
				time = dt*(1:length(this_resp));
				baseline = mean(this_pid(1:b));
				time = time(a:z);
				this_pid = this_pid(a:z);
				this_resp = this_resp(a:z);

				detrended_data(i,k).time = time;
				ff = fit(time(:),this_pid(:),'poly2');
				detrended_data(i,k).stim = this_pid - ff(time)' + mean(ff(time)) - baseline;

				ff = fit(time(:),this_resp(:),'poly2');
				detrended_data(i,k).resp = this_resp - ff(time) + mean(ff(time));
			end
		end
	end
	cache('detrended_data',detrended_data)
end

% remove trend -- and combine by neuron
b = floor(5/dt);
a = floor(35/dt);
z = floor(55/dt);
if ~exist('MSG_per_neuron.mat','file')
	MSG_data.time = [];
	MSG_data.stim = [];
	MSG_data.resp = [];
	for k = 1:max(combined_data.neuron)
		for i = 1:length(paradigm_names)
			plot_these = find(strcmp(paradigm_names{i}, combined_data.paradigm));
			plot_these = intersect(plot_these,find(combined_data.neuron==k));
			if ~isempty(plot_these)
				MSG_data(i,k).stim = NaN(20001,length(plot_these));
				MSG_data(i,k).resp = NaN(20001,length(plot_these));
				for n = 1:length(plot_these)
					this_pid=(combined_data.PID(plot_these(n),:));
					this_resp=(combined_data.fA(:,plot_these(n)));
					time = dt*(1:length(this_resp));
					baseline = mean(this_pid(1:b));
					time = time(a:z);
					this_pid = this_pid(a:z);
					this_resp = this_resp(a:z);
					MSG_data(i,k).time = time;
					ff = fit(time(:),this_pid(:),'poly2');
					this_stim = this_pid - ff(time)' + mean(ff(time)) - baseline;
					MSG_data(i,k).stim(:,n) = this_stim(:);

					ff = fit(time(:),this_resp(:),'poly2');
					this_resp  = this_resp - ff(time) + mean(ff(time));
					MSG_data(i,k).resp(:,n) = this_resp(:);
				end
			end
		end
	end
	save('MSG_per_neuron.mat','MSG_data')
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
peak_loc = NaN(length(detrended_data),3);
mean_stim = NaN(length(detrended_data),3);
clear l 
l = zeros(8,1);
for k = 1:3
	for i = 1:length(detrended_data)
		mean_stim(i,k) = mean(detrended_data(i,k).stim);
		a = detrended_data(i,k).resp - mean(detrended_data(i,k).resp);
		if ~isempty(a)
			a = a/std(a);
			b = detrended_data(i,k).stim - mean(detrended_data(i,k).stim);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = dt*(1:length(x));
			t = t-mean(t);
			x = x/max(x);

			l(i) = plot(t,x,'Color',c(i,:));

			[~,loc] = max(x);
			peak_loc(i,k) = t(loc);
		end
	end
end
set(gca,'XLim',[-.1 .4])
xlabel('Lag (s)')
ylabel('Cross Correlation (norm)')
L = paradigm_names;
for i = 1:length(L)
	L{i} = L{i}(strfind(L{i},'-')+1:end);
end
legend(l,L)
subplot(1,3,3), hold on
for k = 1:3
	plot(mean_stim(:,k),peak_loc(:,k)*1000,'+k')
end
xlabel('Mean Stimulus (V)')
set(gca,'XLim',[0 4])
ylabel('Peak of cross correlation (ms)')

[ff,gof]=fit(mean_stim(~isnan(mean_stim)),peak_loc(~isnan(mean_stim)),'poly1');
l = plot(0:5,1e3*ff(0:5),'--k');
legend(l,strcat('r^2=',oval(gof.rsquare)))

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% Ok, we show through the cross-correlation function that responses speed up on increasing stimulus mean. Can we see the same effect when we back out filters from the data? 

c = parula(length(detrended_data));
figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
subplot(1,3,1:2), hold on
peak_loc = NaN(length(detrended_data),3);
min_loc = NaN(length(detrended_data),3);
mean_stim = NaN(length(detrended_data),3);
allfilters = cache('allfilters');
if isempty(allfilters)
	for k = 1:3
		for i = 1:length(detrended_data)
			allfilters(i,k).K = [];
			allfilters(i,k).p = [];
		end
	end
	clear l 
	l = zeros(8,1);
	for k = 1:3
		for i = 1:length(detrended_data)
			mean_stim(i,k) = mean(detrended_data(i,k).stim);
			a = detrended_data(i,k).resp - mean(detrended_data(i,k).resp);
			if ~isempty(a)
				if isempty(allfilters(i,k).K)
					a = a/std(a);
					b = detrended_data(i,k).stim - mean(detrended_data(i,k).stim);
					b = b/std(b);

					[thisK, ~, filtertime_full] = FindBestFilter(b,a,[],'regmax=10;','regmin=1;','filter_length=1099;','offset=300;');
					thisK = thisK(100:1000);
					allfilters(i,k).K = thisK;
					filtertime = filtertime_full(100:1000);

					% fit a parametric filter to this
					clear d p
					p.   n= 1.4766;
					p.   A= 0.7921;
					p.tau1= 52.3750;
					p.tau2= 35.6094;

					d.stimulus = thisK(200:end);
					d.stimulus = d.stimulus/max(d.stimulus);
					d.response = d.stimulus;
					p = FitModel2Data(@FitFilter,d,p);
					p = FitModel2Data(@FitFilter,d,p);
					allfilters(i,k).p = FitModel2Data(@FitFilter,d,p);
				end
			end
		end
	end
	cache('allfilters',allfilters);
end

for k = 1:3
	for i = 1:length(detrended_data)
		if ~isempty(allfilters(i,k).p)
			K2 = FitFilter(allfilters(i,k).K(200:end),allfilters(i,k).p);
			filtertime = dt*(1:length(K2));
			l(i) = plot(filtertime,K2,'Color',c(i,:));

			[~,loc] = max(K2);
			peak_loc(i,k) = filtertime(loc);
			[~,loc] = min(K2);
			min_loc(i,k) = filtertime(loc);
		end
	end
end


% debug -- check that the fit works. 
% for k = 1:3
% 	for i = 1:8
% 		if ~isempty(allfilters(i,k).K)
% 			figure, hold on
% 			thisK = allfilters(i,k).K(200:end);
% 			thisK = thisK/max(thisK);
% 			plot(thisK)
% 			hold on
% 			K2 = FitFilter(thisK,allfilters(i,k).p);
% 			plot(K2)
% 			pause(2)
% 			close(gcf)
% 		end

% 	end
% end

set(gca,'XLim',[-.1 .5])
xlabel('Lag (s)')
ylabel('Fitler amplitude (norm)')
L = paradigm_names;
for i = 1:length(L)
	L{i} = L{i}(strfind(L{i},'-')+1:end);
end
legend(l,L)

subplot(1,3,3), hold on
plot(mean_stim,peak_loc/dt,'k+')
xlabel('Mean Stimulus (V)')
ylabel('Filter Peak (ms)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%  ######      ###    #### ##    ##       ###    ##    ##    ###    ##       ##    ##  ######  ####  ######  
% ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##        ##  ##  ##    ##  ##  ##    ## 
% ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##         ####   ##        ##  ##       
% ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##          ##     ######   ##   ######  
% ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##          ##          ##  ##        ## 
% ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##          ##    ##    ##  ##  ##    ## 
%  ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ########    ##     ######  ####  ######  

%% Gain Analysis
% In this section, we investigate how the gain of the neuron changes with stimulus mean. Since we have already backed out linear filters from the data, we convolve them with the stimulus to get a linear prediction of the response. We plot the actual neuron response vs. the linear prediciton, to see how the slope (thus, the gain) varies with the different stimuli (here, shown in different colours). 

filtertime = dt*(1:length(allfilters(1,1).K));
filtertime = filtertime - .2;

% first figure out the trivial scaling using the lowest dose
clear trival_scaling
trival_scaling = struct;
for j = 1:3
	detrended_data(1,j).fp = convolve(detrended_data(1,j).time,detrended_data(1,j).stim,allfilters(1,j).K,filtertime);
	x = detrended_data(1,j).fp;
	y = detrended_data(1,j).resp;
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = [];
	y(rm_this) = [];
	trival_scaling(j).cf = fit(x,y,'poly1');
end

for i = 1:length(detrended_data)
	for j = 1:width(detrended_data)
		if ~isempty(allfilters(i,j).K)
			detrended_data(i,j).fp = convolve(detrended_data(i,j).time,detrended_data(i,j).stim,allfilters(i,j).K,filtertime);
			% measure the gain
			x = detrended_data(i,j).fp;
			y = detrended_data(i,j).resp;
			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];
			temp = fit(x,y,'poly1');
			detrended_data(i,j).gain = temp.p1;

			% account for some trivial scaling
			detrended_data(i,j).fp = trival_scaling(j).cf(detrended_data(i,j).fp);
			

		end
	end
end

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
ss = 20;
for j = 1:3
	subplot(1,3,j), hold on
	title(strcat('Group # ',oval(j)))
	set(gca,'XLim',[0 45],'YLim',[0 45])
	xlabel('Linear Prediction (Hz)')
	ylabel('Neuron Response (Hz)')
	for i = 1:length(detrended_data)
		if ~isempty(allfilters(i,j).K)
			plot(detrended_data(i,j).fp(1:ss:end),detrended_data(i,j).resp(1:ss:end),'.','Color',c(i,:))
		end
	end
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% That's quite convincing. Now, we quantify this effect by plotting gain vs. the mean of the stimulus. We also plot the coefficient of variation of the response as a function of the coefficient of variation of the stimulus.  
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for j = 1:3
	for i = 1:length(detrended_data)
		if ~isempty(allfilters(i,j).K)
			plot(mean(detrended_data(i,j).stim),detrended_data(i,j).gain,'k+')
		end
	end
end
set(gca,'XScale','log','YScale','log'); %,'YLim',[5 60],'XLim',[0.4 6])
xlabel('Mean Stimulus (V)')
ylabel('Neuron Gain (Hz/V')

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
