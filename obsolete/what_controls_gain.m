% What controls gain?
% 
% created by Srinivas Gorur-Shandilya at 3:16 , 31 March 2015. Contact me at http://srinivas.gs/contact/
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

%% What controls gain?
% In this document, we analyse what the fast gain modulation depends on purely from the phenomenology. For example, does it depend on the stimulus? Or the derivative of the stimulus? Or the response? 

%%
% First, we demonstrate that gain is indeed controlled on a fast time scale:


clearvars -except being_published fig1 fig2 fig3 fig4 fig5

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

fig_handle=figure('Units','pixels','outerposition',[100 100 1240 720],'Color','w','PaperUnits','points','PaperSize',[1240 720],'Toolbar','none');
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position',[70 490 402.5 105]);
axes_handles(2)=axes('Units','pixels','Position',[542.5 437.5 105 105]);
axes_handles(3)=axes('Units','pixels','Position',[735 560 472.5 105]);
axes_handles(4)=axes('Units','pixels','Position',[735 437.5 472.5 105]);
axes_handles(5)=axes('Units','pixels','Position',[70 70 262.5 262.5]);
axes_handles(6)=axes('Units','pixels','Position',[437.5 70 262.5 262.5]);
axes_handles(7)=axes('Units','pixels','Position',[805 70 350 262.5]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end

% set up a colour map
c = parula(8);


% plot stimulus
plot(axes_handles(1),tA,mean2(PID),'k')
set(axes_handles(1),'XLim',[10 60])
ylabel(axes_handles(1),'Stimulus (V)')

% plot response
plot(axes_handles(3),tA,mean2(fA),'k')
set(axes_handles(3),'XLim',[10 60],'XTickLabel',{});
ylabel(axes_handles(3),'Firing Rate (Hz)')

% extract Linear model for each trial 
K = [];
for i = [3:10 13:20]
	[this_K, ~, filtertime_full] = FindBestFilter(PID(:,i),fA(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	K = [K;this_K];
end

% plot linear filter
axes(axes_handles(2))
plot_this = mean2(K);
err = std(K)/sqrt(width(K));
shadedErrorBar(filtertime,plot_this,err,{'Color',c(3,:)})
xlabel(axes_handles(2),'Lag (s)')
ylabel(axes_handles(2),'Filter K')


% make a linear prediction using a filter fit to the mean data (this is almost exactly the same)
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
l=plot(axes_handles(4),tA,fp,'Color',c(3,:));
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2)))
ylabel(axes_handles(4),'K\otimes stimulus (Hz)')


% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles([5 7]);

hl_min = .1;
hl_max = 10;
history_lengths = logspace(log10(hl_min),log10(hl_max),30);

[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',mean2(fA),'prediction',fp,'stimulus',mean2(PID),'time',tA,'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(9),'use_cache',1,'engine',@GainAnalysis5);
set(axes_handles(7),'XLim',[.09 11]) % to show .1 and 10 on the log scale

% show the p-value
axes(axes_handles(5))
text(10,60,'p < 0.01')

% plot gain vs preceding stimulus


% link some axes
linkaxes(axes_handles([3:4]))

% fix some labels
xlabel(axes_handles(4),'Time (s)')
xlabel(axes_handles(1),'Time (s)')
ylabel(axes_handles(7),'Relative Gain')
set(axes_handles(7),'YScale','log','YTick',[0.5 1 2],'YLim',[0.4 3.5])
set(axes_handles(4),'YLim',[0 100])
ylabel(axes_handles(5),'Firing Rate (Hz)')
xlabel(axes_handles(5),'K\otimes stimulus (Hz)')
title(axes_handles(5),'')

% Indicate regions of high and low stimulus on the stimulus
hl = round(history_lengths(9)*1e3);
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


%% Is gain a function of the preceding stimulus? 
% In this section, we plot the gain as a function of the stimulus in the preceding 489ms. We compute the gain in three ways: on the left, as the instantaneous gain (defined vs. the linear prediction), in the middle, by fitting lines to these segments, and on the right, as the slope of the 1st principal component to snippets of data vs. linear prediction over 489ms. 

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,3,1), hold on

gain_inst = mean2(fA)./fp; 
plot(shat(1e4:end-1e4),gain_inst(1e4:end-1e4),'k.');
ylabel('Instantaneous relative gain')
xlabel('Mean Stimulus')

subplot(1,3,2), hold on
[~,gain_slopes] = MakeFig6G(mean2(PID(1e4:end-1e4,:)),mean2(fA(1e4:end-1e4,:)),fp(1e4:end-1e4),round(history_lengths(9)*1e3));
ylabel('Relative gain (slopes)')
xlabel('Mean Stimulus')
plot(shat(1e4:end-1e4),gain_slopes,'k.');


subplot(1,3,3), hold on
[gain_pca,r2,t] = EstimateGain(fp(1e4:end-1e4),mean2(fA(1e4:end-1e4,:)),hl,1);
t = t*1e-3;
gain_pca = interp1(t,gain_pca,tA);
plot(shat(1e4:end-1e4),gain_pca(1e4:end-1e4),'k.');

ylabel('Relative gain (PCA)')
xlabel('Mean Stimulus')
set(gca,'YLim',[0 15])

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% As can be seen from these plots, there seems to be more than just the box-filtered stimulus that is responsible for the gain control. 

%% Does stimulus or stimulus derivative control gain? 
% Does the ORN change gain as a function of preceding stimulus, or preceding stimulus derivative? To determine this, we attempt to back out a linear filters from the gain (from the slopes) and the stimulus. The filter is parametrised so that it can range from an integrator (meaning gain is controlled by stimulus alone) to a differentiator (meaning gain is controlled by the stimulus derivative alone). Thus, the shape of the filter (and how well it explains the data) is informative of what the ORN cares about in controlling gain. 

%%
% To fit, we take the log of the gain, and fit a parametric filter from the stimulus to log(gain). 

clear d
d.response = gain_slopes;
d.response(d.response<0) = NaN;
d.response = log10(d.response);
d.stimulus = mean2(PID(1e4:end-1e4,:));
d.stimulus = d.stimulus - mean(d.stimulus);

clear p
p.     n= 0.8001;
p.     A= 0.3828;
p.  tau1= 198;
p.  tau2= 198;
p. scale= -0.0072;
p.offset= -35.2434;

[gain_pred,K_g] = pLinearModel(d.stimulus,p);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,3,1:2), hold on
plot(tA(1e4:end-1e4),d.response,'k')
l=plot(tA(1e4:end-1e4),gain_pred,'r');
xlabel('Time (s)')
ylabel('log (gain)')
set(gca,'YLim',[-1 1])
r2 = rsquare(gain_pred,d.response);
legend(l,strcat('r^2=',oval(r2)))

subplot(1,3,3), hold on
plot(K_g,'r')
ylabel('Gain Filter')
xlabel('Lag (ms)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%% Does response or response derivative control gain? 
% Does the ORN change gain as a function of preceding response, or preceding response derivative? To determine this, we attempt to back out a linear filters from the gain (from the slopes) and the response. The filter is parametrised so that it can range from an integrator (meaning gain is controlled by response alone) to a differentiator (meaning gain is controlled by the response derivative alone). Thus, the shape of the filter (and how well it explains the data) is informative of what the ORN cares about in controlling gain. 

%%
% To fit, we take the log of the gain, and fit a parametric filter from the response to log(gain). 

clear d
d.response = gain_slopes;
d.response(d.response<0) = NaN;
d.response = log10(d.response);
d.stimulus = mean2(fA(1e4:end-1e4,:));
d.stimulus = d.stimulus - mean(d.stimulus);

clear p

p.     n= 0.6201;
p.     A= 1.1021;
p.  tau1= 197.3970;
p.  tau2= 199.9998;
p. scale= 3.5249e-05;
p.offset= 7.3670e+03;

[gain_pred,K_g] = pLinearModel(d.stimulus,p);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,3,1:2), hold on
plot(tA(1e4:end-1e4),d.response,'k')
l=plot(tA(1e4:end-1e4),gain_pred,'r');
xlabel('Time (s)')
ylabel('log (gain)')
set(gca,'YLim',[-1 1])
r2 = rsquare(gain_pred,d.response);
legend(l,strcat('r^2=',oval(r2)))

subplot(1,3,3), hold on
plot(K_g,'r')
ylabel('Gain Filter')
xlabel('Lag (ms)')

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

t = toc;
%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))
