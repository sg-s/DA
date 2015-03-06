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


a = floor(15/dt);
z = floor(55/dt);

c = parula(length(paradigm_names));

mean_pid = NaN(length(c),1);

% plot lowest dose stimulus
plot_these=find(strcmp(paradigm_names{1}, combined_data.paradigm));
plot_this = mean2(combined_data.PID(plot_these,:));
time = dt*(1:length(plot_this));
plot(axes_handles(1),time,plot_this,'Color',c(1,:))
ylabel(axes_handles(1),'Stimulus (V)')


% plot lowest dose response
plot_this = mean2(combined_data.fA(:,plot_these));
time = dt*(1:length(plot_this));
plot(axes_handles(2),time,plot_this,'Color',c(1,:))
ylabel(axes_handles(2),'ORN Response (Hz)')

set(axes_handles(1),'XLim',[-1 61])
set(axes_handles(2),'XLim',[-1 61])
xlabel(axes_handles(2),'Time (s)')


% remove trend
b = floor(5/dt);
a = floor(35/dt);
z = floor(55/dt);
detrended_data = cache('detrended_data'); % needs this in the cache. run MeanShiftedGaussians.m to generate this

% load the filters
allfilters = cache('allfilters');

% plot the filter for the lowest dose 
filtertime = dt*(1:length(allfilters(1,1).K));
filtertime = filtertime - 200*dt;
plot(axes_handles(3),filtertime,allfilters(1,1).K,'Color',c(1,:))
set(axes_handles(3),'XLim',[min(filtertime) max(filtertime)])
xlabel(axes_handles(3),'Lag (s)')
ylabel(axes_handles(3),'Filter')

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

% plot linear prediction vs. data
ss = 20;
plot(axes_handles(4),detrended_data(1,1).fp(1:ss:end),detrended_data(1,1).resp(1:ss:end),'.','Color',c(1,:))
xlabel(axes_handles(4),'K \otimes s')
ylabel(axes_handles(4),'Response (Hz)')


% plot the stimulus distributions 
a = floor(15/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	plot_hist = (combined_data.PID(plot_these,a:z));
	[hy,hx]  = hist(plot_hist(:),50);
	hy = hy/max(hy);
	plot(axes_handles(5),hx,hy,'Color',c(i,:));
end

xlabel(axes_handles(5),'PID (V)')
ylabel(axes_handles(5),'count (norm)')


% plot the response distributions 
a = floor(15/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	temp =  [detrended_data(i,:).resp];
	temp = temp(:);
	[hy,hx]  = hist(temp,50); % this is being 
	hy = hy/max(hy);
	plot(axes_handles(6),hx,hy,'Color',c(i,:));
end

xlabel(axes_handles(6),'Response (Hz)')
ylabel(axes_handles(6),'count (norm)')


% show gain changes -- change in slope of scatter plot
ss = 50;
for j = 1:3
	for i = 1:length(detrended_data)
		if ~isempty(allfilters(i,j).K)
			plot(axes_handles(7),detrended_data(i,j).fp(1:ss:end),detrended_data(i,j).resp(1:ss:end),'.','Color',c(i,:))
		end
	end
end

set(axes_handles(7),'XLim',[0 45],'YLim',[0 45])
xlabel(axes_handles(7),'Linear Prediction (Hz)')
ylabel(axes_handles(7),'Neuron Response (Hz)')

% show gain changes -- gain vs. mean stimulus

x = []; y = [];
for j = 1:3
	for i = 1:length(detrended_data)
		if ~isempty(allfilters(i,j).K)
			plot(axes_handles(8),mean(detrended_data(i,j).stim),detrended_data(i,j).gain,'+','Color',c(i,:));
			x = [x mean(detrended_data(i,j).stim)];
			y = [y detrended_data(i,j).gain];
		end
	end
end
cf = fit(x(:),y(:),'power1');
set(axes_handles(8),'XScale','log','YScale','log'); %,'YLim',[5 60],'XLim',[0.4 6])
xlabel(axes_handles(8),'Mean Stimulus (V)')
ylabel(axes_handles(8),'Neuron Gain (Hz/V')
l=plot(axes_handles(8),x(:),cf(x(:)),'k');
legend(l,strcat('\alpha=',oval(cf.b)))


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end

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


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end


%% Figure 3

clearvars -except being_published
fig_handle=figure('Units','pixels','outerposition',[66 5 661 871],'PaperUnits','points','PaperSize',[661 871],'Color','w','Toolbar','none');
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position',[63.825 744.625 574.425 85.1]);
axes_handles(2)=axes('Units','pixels','Position',[63.825 595.7 574.425 127.65]);
axes_handles(3)=axes('Units','pixels','Position',[63.825 489.325 574.425 85.1]);
axes_handles(4)=axes('Units','pixels','Position',[63.825 276.575 191.475 170.2]);
axes_handles(5)=axes('Units','pixels','Position',[382.95 276.575 255.3 170.2]);
axes_handles(6)=axes('Units','pixels','Position',[63.825 42.55 191.475 191.475]);
axes_handles(7)=axes('Units','pixels','Position',[382.95 42.55 212.75 191.475]);
for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end


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

% plot prediction and prediction quality
l=plot(axes_handles(2),tA,fp,'r');
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2)))

% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles(4:5);
GainAnalysisWrapper(mean2(fA),fp,mean2(PID),tA,0.4290,ph);

% plot gain vs preceding stimulus
[x,y] = MakeFig6G(mean2(PID),mean2(fA),fp,400);
%c = MakeFig6H(mean2(PID),mean2(fA),400);
gain_time = mean(diff(tA))*(1:length(x));
rm_this = (isnan(x) | isnan(y));
x(rm_this) = [];
y(rm_this) = [];
gain_time(rm_this) = [];
ss = 50;
plot(axes_handles(6),x(1:ss:end),y(1:ss:end),'k.')
xlabel(axes_handles(6),'Stimulus in preceding 400ms')
ylabel(axes_handles(6),'Instantaneous gain')

% plot gain
plot(axes_handles(3),gain_time,y,'r')
ylabel(axes_handles(3),'Gain')
set(axes_handles(3),'XLim',[10 60],'YLim',[0 7])
xlabel(axes_handles(3),'Time (s)')

% show adaptation of dynamics

% fix some labels
ylabel(axes_handles(5),'Gain')
xlabel(axes_handles(4),'Linear Prediction (Hz)')
ylabel(axes_handles(4),'Neuron Response (Hz)')

PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

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



