% Paper Figures
% makes all the figures for the paper
% 
% created by Srinivas Gorur-Shandilya at 12:57 , 21 January 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%    ######## ####  ######   ##     ## ########  ########       ##   
%    ##        ##  ##    ##  ##     ## ##     ## ##           ####   
%    ##        ##  ##        ##     ## ##     ## ##             ##   
%    ######    ##  ##   #### ##     ## ########  ######         ##   
%    ##        ##  ##    ##  ##     ## ##   ##   ##             ##   
%    ##        ##  ##    ##  ##     ## ##    ##  ##             ##   
%    ##       ####  ######    #######  ##     ## ########     ###### 


%% Figure 1: ORNs decrease gain on increasing stimulus mean
clearvars -except being_published
fig_handle=figure('Units','pixels','outerposition',[81 5 599 871],'PaperUnits','points','PaperSize',[599 871],'Color','w','Toolbar','none');
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position',[63.825 744.625 489.325 85.1]);
axes_handles(2)=axes('Units','pixels','Position',[63.825 616.975 489.325 106.375]);
axes_handles(3)=axes('Units','pixels','Position',[63.825 425.5 170.2 170.2]);
axes_handles(4)=axes('Units','pixels','Position',[382.95 425.5 170.2 170.2]);
axes_handles(5)=axes('Units','pixels','Position',[63.825 276.575 212.75 106.375]);
axes_handles(6)=axes('Units','pixels','Position',[340.4 276.575 212.75 106.375]);
axes_handles(7)=axes('Units','pixels','Position',[63.825 42.55 212.75 191.475]);
axes_handles(8)=axes('Units','pixels','Position',[340.4 42.55 212.75 191.475]);
for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end

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


c = parula(length(paradigm_names)+1);

mean_pid = NaN(length(c),1);

% plot lowest dose stimulus
plot_these=find(strcmp(paradigm_names{1}, combined_data.paradigm));
plot_this = mean2(combined_data.PID(plot_these,:));
err = std(combined_data.PID(plot_these,:));
err = err/sqrt(length(plot_these));
time = dt*(1:length(plot_this));
plot(axes_handles(1),time,plot_this,'Color',c(1,:))
% axes(axes_handles(1))
% shadedErrorBar(time,plot_this,err,{'Color',c(1,:)})
ylabel(axes_handles(1),'Stimulus (V)')


% plot lowest dose response
plot_this = mean2(combined_data.fA(:,plot_these));
err = std(combined_data.fA(:,plot_these)');
err = err/sqrt(length(plot_these));
time = dt*(1:length(plot_this));
plot(axes_handles(2),time,plot_this,'Color',c(1,:))
ylabel(axes_handles(2),'ORN Response (Hz)')
% axes(axes_handles(2))
% shadedErrorBar(time,plot_this,err,{'Color',c(1,:)})
set(axes_handles(1),'XLim',[15 55])
set(axes_handles(2),'XLim',[15 55])
xlabel(axes_handles(2),'Time (s)')


% remove trend
load('MSG_per_neuron.mat','MSG_data')

% back out all filters
% for i = 1:8
% 	for j = 1:13
% 		if width(MSG_data(i,j).stim) > 1
% 			this_stim = mean2(MSG_data(i,j).stim);
% 			this_resp = mean2(MSG_data(i,j).resp);
% 			[K, ~, filtertime_full] = FindBestFilter(this_stim,this_resp,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
% 			filtertime_full = filtertime_full;
% 			filtertime = (-200:900);
% 			K = interp1(filtertime_full,K,filtertime);
% 			MSG_data(i,j).K = K/max(K);
% 		end
% 	end
% end
% save('MSG_per_neuron.mat','MSG_data')

% plot the filter for the lowest dose 
filtertime = (-200:900)*1e-3;
K = 0*filtertime;
for i = 1:13 % get all the filters for the lowest dose
	K = [K; MSG_data(1,i).K];
end
err = std(K);
err = err/sqrt(width(K));
axes(axes_handles(3))
shadedErrorBar(filtertime,mean2(K),err,{'Color',c(1,:)})
set(axes_handles(3),'XLim',[min(filtertime) max(filtertime)])
xlabel(axes_handles(3),'Lag (s)')
ylabel(axes_handles(3),'Filter K (norm)')

% make linear predictions everywhere
for i = 1:8
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			this_stim = mean2(MSG_data(i,j).stim);
			this_resp = mean2(MSG_data(i,j).resp);
			MSG_data(i,j).fp = convolve(MSG_data(i,j).time,mean2(MSG_data(i,j).stim),MSG_data(i,j).K,filtertime);
		end
	end
end


% first figure out the trivial scaling using the lowest dose
% clear trival_scaling
% trival_scaling = struct;
% for j = 1:3
% 	detrended_data(1,j).fp = convolve(detrended_data(1,j).time,detrended_data(1,j).stim,allfilters(1,j).K,filtertime);
% 	x = detrended_data(1,j).fp;
% 	y = detrended_data(1,j).resp;
% 	rm_this = isnan(x) | isnan(y);
% 	x(rm_this) = [];
% 	y(rm_this) = [];
% 	trival_scaling(j).cf = fit(x,y,'poly1');
% end

% for i = 1:length(detrended_data)
% 	for j = 1:width(detrended_data)
% 		if ~isempty(allfilters(i,j).K)
% 			detrended_data(i,j).fp = convolve(detrended_data(i,j).time,detrended_data(i,j).stim,allfilters(i,j).K,filtertime);
% 			% measure the gain
% 			x = detrended_data(i,j).fp;
% 			y = detrended_data(i,j).resp;
% 			rm_this = isnan(x) | isnan(y);
% 			x(rm_this) = [];
% 			y(rm_this) = [];
% 			temp = fit(x,y,'poly1');
% 			detrended_data(i,j).gain = temp.p1;

% 			% account for some trivial scaling
% 			detrended_data(i,j).fp = trival_scaling(j).cf(detrended_data(i,j).fp);
% 		end
% 	end
% end

% plot linear prediction vs. data for the lowest dose. 
ss = 25;
y = mean2([MSG_data(1,:).resp]);
x = mean2([MSG_data(1,:).fp]);
plot(axes_handles(4),x(1:ss:end),y(1:ss:end),'.','Color',c(1,:))


xlabel(axes_handles(4),'K \otimes s')
ylabel(axes_handles(4),'Response (Hz)')


% plot the stimulus distributions 
a = floor(15/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	plot_hist = (combined_data.PID(plot_these,a:z));
	[hy,hx]  = hist(plot_hist(:),50);
	hy = hy/sum(hy);
	plot(axes_handles(5),hx,hy,'Color',c(i,:));
end

xlabel(axes_handles(5),'PID (V)')
ylabel(axes_handles(5),'p.d.f')


% plot the response distributions 
a = floor(15/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	temp =  [MSG_data(i,:).resp];
	temp = mean2(temp);
	[hy,hx]  = hist(temp,50);
	hy = hy/sum(hy);
	plot(axes_handles(6),hx,hy,'Color',c(i,:));
end

xlabel(axes_handles(6),'Response (Hz)')
ylabel(axes_handles(6),'p.d.f')


% show gain changes for all paradigms -- average over neurons 
ss = 50;
for i = 1:8 % iterate over all paradigms 
	y = ([MSG_data(i,:).resp]);
	x = ([MSG_data(i,:).fp]);
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end 
	plot(axes_handles(7),x(1:ss:end),y(1:ss:end),'.','Color',c(i,:))
end

xlabel(axes_handles(7),'K\otimes s')
ylabel(axes_handles(7),'Neuron Response (Hz)')


% compute gain changes on a per-neuron basis

gain = NaN(8,13);
mean_stim = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	this_x = [];
	this_y = [];
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			x = MSG_data(i,j).fp;
			if ~isvector(x)
				x = mean2(x);
			end
			if ~isvector(y)
				y = mean2(y);
			end 
			
			gain(i,j) = EstimateGain2(x,y);
			mean_stim(i,j) = mean(mean([MSG_data(i,:).stim]));
		end
	end	
end

% show gain changes -- gain vs. mean stimulus
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		plot(axes_handles(8),mean_stim(i,j),gain(i,j),'+','Color',c(i,:));
	end
end

mean_stim = mean_stim(~isnan(mean_stim));
gain = gain(~isnan(gain));


cf = fit(mean_stim(:),gain(:),'power1');
set(axes_handles(8),'XScale','log','YScale','log','YLim',[1e-3 2e-2])
xlabel(axes_handles(8),'Mean Stimulus (V)')
ylabel(axes_handles(8),'Neuron Gain (Hz/V)')
l=plot(axes_handles(8),sort(mean_stim),cf(sort(mean_stim)),'k');
legend(l,strcat('\alpha=',oval(cf.b)))

PrettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')


if being_published
	snapnow
	delete(gcf)
end

return

%      ######## ####  ######   ##     ## ########  ########     #######  
%      ##        ##  ##    ##  ##     ## ##     ## ##          ##     ## 
%      ##        ##  ##        ##     ## ##     ## ##                 ## 
%      ######    ##  ##   #### ##     ## ########  ######       #######  
%      ##        ##  ##    ##  ##     ## ##   ##   ##          ##        
%      ##        ##  ##    ##  ##     ## ##    ##  ##          ##        
%      ##       ####  ######    #######  ##     ## ########    ######### 


%% Figure 2: ORNs speed up responses on increasing stimulus mean

figure('outerposition',[0 0 1400 900],'PaperUnits','points','PaperSize',[1400 900]); hold on

% filters for mean shifted gaussians
clear l 
l = zeros(8,1);
peak_loc_K = NaN(length(detrended_data),3);
subplot(2,3,1), hold on
for k = 1:3
	for i = 1:length(detrended_data)
		if ~isempty(allfilters(i,k).p)
			K2 = FitFilter(allfilters(i,k).K(200:end),allfilters(i,k).p);
			filtertime = dt*(1:length(K2));
			l(i)=plot(filtertime,K2,'Color',c(i,:));

			[~,loc] = max(K2);
			peak_loc_K(i,k) = filtertime(loc);
		end
	end
end
set(gca,'XLim',[-.01 .5])
xlabel('Lag (s)')
ylabel('Filter')
L = paradigm_names;
for i = 1:length(L)
	L{i} = L{i}(strfind(L{i},'-')+1:end);
end
legend(l,L)


% xcorr for mean shifted gaussians
subplot(2,3,2), hold on
peak_loc_xcorr = NaN(length(detrended_data),3);
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

			plot(t,x,'Color',c(i,:));

			[~,loc] = max(x);
			peak_loc_xcorr(i,k) = t(loc);
		end
	end
end

set(gca,'XLim',[-.1 .3])
xlabel('Lag (s)')
ylabel('Cross Correlation (norm)')


subplot(2,3,3), hold on
clear l
l(1) = plot(mean_stim(:),peak_loc_xcorr(:)/dt,'k+');
l(2) = plot(mean_stim(:),peak_loc_K(:)/dt,'ko');

ff=fit(mean_stim(~isnan(mean_stim)),peak_loc_K(~isnan(mean_stim))/dt,'poly1');
plot(mean_stim(~isnan(mean_stim)),ff(mean_stim(~isnan(mean_stim))),'k')
ff=fit(mean_stim(~isnan(mean_stim)),peak_loc_xcorr(~isnan(mean_stim))/dt,'poly1');
plot(mean_stim(~isnan(mean_stim)),ff(mean_stim(~isnan(mean_stim))),'k')

ylabel('Peak time (ms)')
xlabel('Mean Stimulus (V)')
legend(l,{'Cross correlation','Filter'})

% now we make the plots for the pulse responses.
clearvars -except being_published
load('/local-data/DA-paper/carlotta-martelli/fig3/abc.mat')

% combine all spikes
all_spikes = [];
paradigm = [];
for i = 1:length(spikes)
	if isempty(spikes(i).discard)
		all_spikes = vertcat(all_spikes,spikes(i).A);
		paradigm = [paradigm i*ones(1,width(spikes(i).A))];
	else
		temp = spikes(i).A;
		rm_this = find(spikes(i).discard);
		rm_this(rm_this>width(spikes(i).discard)) = [];
		temp(rm_this,:) = [];
		all_spikes = vertcat(all_spikes,temp);
		paradigm = [paradigm i*ones(1,width(temp))];
	end
end

time = 1e-4*(1:length(all_spikes));

hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

t = 1e-3*(1:length(fA));


% remove bizzare outliers
rm_this(1) = 8+ find(paradigm==6,1,'first');
rm_this(2) = 2+ find(paradigm==8,1,'first');
rm_this = [rm_this find(max(fA)==0)];
fA(:,rm_this) = [];
paradigm(rm_this) = [];

% combine all PID traces. the data is very spotty, so we have to be careful. 
all_PID = [];
paradigm_PID = [];
for i = 1:length(spikes)
	temp = data(i).PID;
	temp = temp(find(temp(:,1)),:);
	all_PID = vertcat(all_PID,temp);
	paradigm_PID = [paradigm_PID i*ones(1,width(temp))];
end
all_PID = all_PID';

% filter and subsample
for i = 1:width(all_PID)
	all_PID(:,i) = filter(ones(30,1)/30,1,all_PID(:,i));
end
all_PID = all_PID(1:10:end,:);

% remove bizzare outliers
rm_this = [];
rm_this(1) = find(paradigm_PID==9,1,'first');
rm_this(2) = find(paradigm_PID==10,1,'first');
all_PID(:,rm_this) = [];
paradigm_PID(rm_this) = [];

% plot PID
subplot(2,3,4), hold on
for i = 1:max(paradigm)
	plot_these = find(paradigm_PID==i);
	plot(t,mean2(all_PID(:,plot_these)));
end

set(gca,'XLim',[.5 3],'YScale','log','YLim',[7e-2 15])
xlabel('Time (s)')
ylabel('Stimulus (V)')


% plot neuron responses
subplot(2,3,5), hold on
for i = 1:max(paradigm)
	plot_these = find(paradigm==i);
	plot(t,mean2(fA(:,plot_these)));
end

set(gca,'XLim',[.5 3])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')


load('2ac_timing.mat')
a = background_stim == 0;
foreground_stim(foreground_stim<1e-2) = NaN; % not reliable

% plot time to peak for the data
subplot(2,3,6), hold on
plot(foreground_stim(a),resp_half_time(a)-stim_half_time(a),'k+')
set(gca,'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
xlabel('Mean Stimulus (V)')
ylabel('\tau_{ORN}-\tau_{PID} (ms)','interpreter','tex')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%       ######## ####  ######   ##     ## ########  ########     #######  
%       ##        ##  ##    ##  ##     ## ##     ## ##          ##     ## 
%       ##        ##  ##        ##     ## ##     ## ##                 ## 
%       ######    ##  ##   #### ##     ## ########  ######       #######  
%       ##        ##  ##    ##  ##     ## ##   ##   ##                 ## 
%       ##        ##  ##    ##  ##     ## ##    ##  ##          ##     ## 
%       ##       ####  ######    #######  ##     ## ########     #######  


%% Figure 3: ORNs can change gain on a fast time scale 

clearvars -except being_published

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

fig_handle=figure('Units','pixels','outerposition',[100 143 1241 677],'Color','w','PaperUnits','points','PaperSize',[1241 677],'Toolbar','none','Menubar','none');
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position',[295.65 476.325 558.45 147.825]);
axes_handles(2)=axes('Units','pixels','Position',[295.65 328.5 558.45 131.4]);
axes_handles(3)=axes('Units','pixels','Position',[295.65 180.675 558.45 131.4]);
axes_handles(4)=axes('Units','pixels','Position',[49.275 328.5 180.675 180.675]);
axes_handles(5)=axes('Units','pixels','Position',[936.225 361.35 279.225 262.8]);
axes_handles(6)=axes('Units','pixels','Position',[936.225 49.275 279.225 180.675]);
axes_handles(7)=axes('Units','pixels','Position',[295.65 49.275 558.45 114.975]);
axes_handles(8)=axes('Units','pixels','Position',[49.275 49.275 180.675 180.675]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end

% set up a colour map
c = parula(8);


% plot stimulus
plot(axes_handles(1),tA,mean2(PID),'k')
set(axes_handles(1),'XLim',[10 60],'XTickLabel',{})
ylabel(axes_handles(1),'Stimulus (V)')

% plot response
plot(axes_handles(2),tA,mean2(fA),'k')
set(axes_handles(2),'XLim',[10 60],'XTickLabel',{})
ylabel(axes_handles(2),'Firing Rate (Hz)')

% extract Linear model
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);
K = K/max(K);
fp = convolve(tA,mean2(PID),K,filtertime);
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;

% plot linear filter
plot(axes_handles(4),filtertime,K,'Color',c(3,:))
xlabel(axes_handles(4),'Lag (s)')
ylabel(axes_handles(4),'Filter (norm)')

% plot prediction and prediction quality
l=plot(axes_handles(3),tA,fp,'Color',c(3,:));
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2)))
set(axes_handles(3),'XTickLabel',{})
ylabel(axes_handles(3),'K\otimes stimulus (Hz)')


% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles(5:6);
history_lengths = (3*floor(1000*logspace(-1,1,30)/3))/1e3;
[~,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',mean2(fA),'prediction',fp,'stimulus',mean2(PID),'time',tA,'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(11));
set(axes_handles(6),'XLim',[.09 11]) % to show .1 and 10 on the log scale

% plot gain vs preceding stimulus
[x,y] = MakeFig6G(mean2(PID),mean2(fA),fp,history_lengths(11)*1e3);
gain_time = mean(diff(tA))*(1:length(x));
rm_this = (isnan(x) | isnan(y));
x(rm_this) = [];
y(rm_this) = [];
gain_time(rm_this) = [];
ss = 50;
plot(axes_handles(8),x(1:ss:end),y(1:ss:end),'k.')
xlabel(axes_handles(8),'Stimulus in preceding 489ms')
ylabel(axes_handles(8),'Relative gain')

% plot gain
gain = y;
plot(axes_handles(7),gain_time,gain,'k')
ylabel(axes_handles(7),'Relative Gain')
set(axes_handles(7),'XLim',[10 60],'YLim',[0 7])
xlabel(axes_handles(7),'Time (s)')


% link some axes
linkaxes(axes_handles([1:3 7]),'x')
linkaxes(axes_handles([2:3]),'y')

% fix some labels
ylabel(axes_handles(6),'Relative Gain')
set(axes_handles(2),'YLim',[0 100])
ylabel(axes_handles(5),'Firing Rate (Hz)')
xlabel(axes_handles(5),'K\otimes stimulus (Hz)')
title(axes_handles(5),'')

% Indicate regions of high and low stimulus on the stimulus
hl = history_lengths(11)*1e3;
shat = ComputeSmoothedStimulus(mean2(PID),hl);

n = floor(sum(~isnan(mean2(fA)))*.33);
shat(1:hl) = Inf; % the initial segment where we can't estimate shat is excluded
shat(isnan(shat)) = Inf;
shat(isnan(mean2(fA))) = Inf;
[~, t_low] = sort(shat,'ascend');
t_low = t_low(1:n); % this is an index
t_low = tA(t_low); % t_low is now a time. 
 
shat = ComputeSmoothedStimulus(mean2(PID),hl);
shat(1:hl) = -Inf;
shat(isinf(shat)) = -Inf;
shat(isnan(mean2(fA))) = -Inf;
[~, t_high] = sort(shat,'descend');
t_high = t_high(1:n);
t_high  = tA(t_high);

plot(axes_handles(1),t_low,1+0*t_low,'g.')
plot(axes_handles(1),t_high,1+0*t_low,'r.')


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end

%      ######## ####  ######   ##     ## ########  ########    ##        
%      ##        ##  ##    ##  ##     ## ##     ## ##          ##    ##  
%      ##        ##  ##        ##     ## ##     ## ##          ##    ##  
%      ######    ##  ##   #### ##     ## ########  ######      ##    ##  
%      ##        ##  ##    ##  ##     ## ##   ##   ##          ######### 
%      ##        ##  ##    ##  ##     ## ##    ##  ##                ##  
%      ##       ####  ######    #######  ##     ## ########          ##  

%% Figure 4: Gain Control is widely observed

clearvars -except being_published
load('CMData_Gain.mat')
load('CM_Data_filters.mat')
combined_data_file = ('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat');
load(combined_data_file)
filtertime = -200:700;
filtertime = filtertime*1e-3;
history_lengths = (3*floor(1000*logspace(-1,1,30)/3))/1e3;


s = 860;
figure('outerposition',[0 0 s s],'PaperUnits','points','PaperSize',[s s]); hold on, clear s
axes_handles = [];
for i = 1:9
	axes_handles(i) = subplot(3,3,i); hold on
end

% first row: experimental replicates
do_these = [2 7 9 13 15 16];

for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	plot(axes_handles(1),filtertime,mean2(K))
end
clear ph
ph(3:4) = axes_handles(2:3);
for i = do_these
	[~,ehl]=max(gain_data(i).low_slopes - gain_data(i).high_slopes);
	time = data(i).dt*(1:length(data(i).PID(1e4:end,1)));
	response = mean2(data(i).fA(1e4:end,:));
	stimulus = mean2(data(i).PID(1e4:end,:));
	prediction = mean2(data(i).LinearFit(1e4:end,:));
	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph);
end


% find the right place to clip the x axis
c = [];
for i = do_these
	c = [c find(gain_data(i).low_gof > .85 & gain_data(i).high_gof > .85,1,'first')];
end
set(ph(4),'XLim',[history_lengths(floor(mean(c))) 11],'YLim',[.5 1.5])

% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')


%    ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%    ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%    ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%    ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%    ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
   
%    ##    ## ######## ##     ## ########   #######  ##    ##  ######  
%    ###   ## ##       ##     ## ##     ## ##     ## ###   ## ##    ## 
%    ####  ## ##       ##     ## ##     ## ##     ## ####  ## ##       
%    ## ## ## ######   ##     ## ########  ##     ## ## ## ##  ######  
%    ##  #### ##       ##     ## ##   ##   ##     ## ##  ####       ## 
%    ##   ### ##       ##     ## ##    ##  ##     ## ##   ### ##    ## 
%    ##    ## ########  #######  ##     ##  #######  ##    ##  ######  


do_these = [11 17];
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(7),filtertime,mean2(K))];
end
legend(l,{'pb1A','ab3A'})

clear ph
ph(3:4) = axes_handles(8:9);
for i = do_these
	[~,ehl]=max(gain_data(i).low_slopes - gain_data(i).high_slopes);
	time = data(i).dt*(1:length(data(i).PID(1e4:end,1)));
	response = mean2(data(i).fA(1e4:end,:));
	stimulus = mean2(data(i).PID(1e4:end,:));
	prediction = mean2(data(i).LinearFit(1e4:end,:));
	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph);
end

% find the right place to clip the x axis
c = [];
for i = do_these
	c = [c find(gain_data(i).low_gof > .85 & gain_data(i).high_gof > .85,1,'first')];
end
set(ph(4),'XLim',[history_lengths(floor(mean(c))) 11],'YLim',[.5 1.5])

% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')

%     ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%     ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%     ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%     ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%     ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
    
%      #######  ########   #######  ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##    ## 
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##       
%     ##     ## ##     ## ##     ## ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##   ##         ## 
%     ##     ## ##     ## ##     ## ##     ## ##    ##  ##    ## 
%      #######  ########   #######   #######  ##     ##  ######  


do_these = [7 8 10 14 17 18];
odours = {'1but','1o3ol','dsucc','2ac','2but','5ol'};
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(4),filtertime,mean2(K))];
end
legend(l,odours)

clear ph
ph(3:4) = axes_handles(5:6);
for i = do_these
	[~,ehl]=max(gain_data(i).low_slopes - gain_data(i).high_slopes);
	time = data(i).dt*(1:length(data(i).PID(1e4:end,1)));
	response = mean2(data(i).fA(1e4:end,:));
	stimulus = mean2(data(i).PID(1e4:end,:));
	prediction = mean2(data(i).LinearFit(1e4:end,:));
	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph);
end


% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')

% find the right place to clip the x axis
c = [];
for i = do_these
	c = [c find(gain_data(i).low_gof > .85 & gain_data(i).high_gof > .85,1,'first')];
end
set(ph(4),'XLim',[history_lengths(floor(mean(c))) 11],'YLim',[.5 1.5])


% cosmetics
xlabel(axes_handles(1),'Lag (s)')
xlabel(axes_handles(4),'Lag (s)')
title(axes_handles(1),'Filters')
title(axes_handles(4),'Filters')
title(axes_handles(7),'Filters')
ylabel(axes_handles(3),'Relative Gain')
ylabel(axes_handles(6),'Relative Gain')
ylabel(axes_handles(9),'Relative Gain')
title(axes_handles(2),strcat('T_H=',oval(history_lengths(10))))
title(axes_handles(5),strcat('T_H=',oval(history_lengths(10))))
title(axes_handles(8),strcat('T_H=',oval(history_lengths(10))))

for j = 1:length(axes_handles)
	% remove all the scatter points
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		if strcmp(get(h(i),'Marker'),'.')
			rm_this = [rm_this i];
		end
	end
	delete(h(rm_this))

	% remove the line indicating the example history plot
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		if  strcmp(get(h(i),'LineStyle'),'-.')
			rm_this = [rm_this i];
		end
	end
	delete(h(rm_this))
end


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

ylabel(axes_handles(1),'Exp. Replicates','FontSize',20)
ylabel(axes_handles(7),'Diff. ORNs','FontSize',20)
ylabel(axes_handles(4),'Diff. odors','FontSize',20)

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

t = toc;
%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))



