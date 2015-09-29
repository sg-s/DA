% Fast Gain Control
% 
% created by Srinivas Gorur-Shandilya at 3:28 , 24 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
    path1 = [path1 ':usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called 
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic


%         ##    ##    ###    ######## ##     ## ########     ###    ##       
%         ###   ##   ## ##      ##    ##     ## ##     ##   ## ##   ##       
%         ####  ##  ##   ##     ##    ##     ## ##     ##  ##   ##  ##       
%         ## ## ## ##     ##    ##    ##     ## ########  ##     ## ##       
%         ##  #### #########    ##    ##     ## ##   ##   ######### ##       
%         ##   ### ##     ##    ##    ##     ## ##    ##  ##     ## ##       
%         ##    ## ##     ##    ##     #######  ##     ## ##     ## ######## 

%          ######  ######## #### ##     ## ##     ## ##       #### 
%         ##    ##    ##     ##  ###   ### ##     ## ##        ##  
%         ##          ##     ##  #### #### ##     ## ##        ##  
%          ######     ##     ##  ## ### ## ##     ## ##        ##  
%               ##    ##     ##  ##     ## ##     ## ##        ##  
%         ##    ##    ##     ##  ##     ## ##     ## ##        ##  
%          ######     ##    #### ##     ##  #######  ######## ####




%% Figure 4: ORNs can change gain on a fast time scale 
% Sensors that control gain in response to changing stimulus statistics must do so relevant timescales: a gain control mechanism that is too slow risks saturation or insensitivity, while a gain control mechanism that is too quick to change leads to fluctuating responses. To mimic an odor signal that changes its statistics rapidly, we recorded ORN responses to a flickering stimulus with a large variance (A-B).  The linear filter extracted from this data (C) convolved with the stimulus yields a linear prediction (D). Comparing the responses (y-axis) to the linear prediction (x-axis) when the stimulus is locally low (green) vs. locally high (red) shows that neuron gain at these two times is significantly different (E). Within this single stimulus presentation, gain varies inversely with the recent stimulus (F). Repeating the analysis over a range of time scales identifies a region of sub-second timescales over which ORNs significantly change their gain in a stimulus-driven manner. The odour used is ethyl acetate, and recordings are from ab3A. The vertical line in G indicates the history length used in E. 

clearvars -except being_published 

% make figure placeholders 
fig_handle=figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
clf(fig_handle);
axes_handles(1) = subplot(10,3,1:9); 
axes_handles(2) = subplot(10,3,10:18); 
axes_handles(3) = subplot(10,3,[19:3:28]);
axes_handles(4) = subplot(10,3,[20:3:29]);
axes_handles(5) = subplot(10,3,[21:3:30]);

movePlot(axes_handles(1),'up',.03)
movePlot(axes_handles(2),'up',.04)
movePlot(axes_handles(3),'down',.02)
movePlot(axes_handles(4),'down',.02)
movePlot(axes_handles(5),'down',.02)

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end

load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;
B_spikes = spikes(2).B;


% A spikes --> firing rate
hash = dataHash(full(all_spikes));
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

axes(axes_handles(1)); hold on
errorShade(tA(1:10:end),mean2(PID(1:10:end,:)),sem(PID(1:10:end,:)),'Color',[0.2 .2 .2]);
set(axes_handles(1),'XLim',[0 70])
ylabel('Stimulus (V)')

axes(axes_handles(3));
y = zeros(300,width(PID));
for i = 1:width(PID)
	[y(:,i),x] = histcounts(PID(:,i),300);x(1) = [];
end
errorShade(x,mean2(y),sem(y),'Color',[.2 .2 .2]);
warning off % because there are some -ve values on the log scale
set(axes_handles(3),'XScale','log','YScale','log','XLim',[.1 10],'YTick',[1e1 1e3 1e5])
xlabel('Stimulus (V)')
ylabel('count')
warning on

% make a linear filter
R = mean2(fA);
[K, filtertime_full] = fitFilter2Data(mean2(PID),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean2(PID),K,filtertime);

% plot the response and the prediction
clear l
[ax,plot1,plot2] = plotyy(axes_handles(2),tA,R,tA,fp);
set(ax(1),'XLim',[0 70])
set(ax(2),'XLim',[0 70],'YLim',[min(fp) max(fp)])
set(plot1,'Color','k')
set(plot2,'Color','r')
ylabel(ax(1),'ORN Response (Hz)')
ylabel(ax(2),'Projected Stimulus')
set(axes_handles(2),'box','off')
legend([plot1 plot2],'Neuron Response',strcat('Projected Stimulus, r^2=',oval(rsquare(mean2(fA),fp))),'Location','northwest');
xlabel(axes_handles(2),'Time (s)')

shat = computeSmoothedStimulus(mean2(PID),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
axes(axes_handles(4))
ss = 1;
cc = parula(100);
c= cc(shat,:);
scatter(fp(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')
xlabel(axes_handles(4),'Projected Stimulus')
ylabel(axes_handles(4),'Actual response (Hz)')
shat = computeSmoothedStimulus(mean2(PID),500);
colorbar;
caxis([min(shat) max(shat)]);

% plot gain vs stimulus for all these whiffs
axes(axes_handles(5))
fp = convolve(tA,mean2(PID),K,filtertime);
shat = computeSmoothedStimulus(mean2(PID),500);

% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = computeOnsOffs(R>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff=fit(fp(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];

errorbar(axes_handles(5),mean_stim,gain,gain_err,'k+')


options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./gain_err;
ff = fit(mean_stim(:),gain(:),'power1',options);
l(1) = plot(axes_handles(5),sort(mean_stim),ff(sort(mean_stim)),'r');
L{1} = ['y=\alpha x^{-1}',char(10),'r^2=',oval(rsquare(ff(mean_stim),gain))];
legend(l,L)
xlabel(axes_handles(5),'Stimulus in preceding 500ms (V)')
ylabel(axes_handles(5),'Gain (Hz/V)')

% cosmetic fixes
set(axes_handles(1),'XTick',[])
set(axes_handles(4),'XLim',[min(fp) max(fp)])
set(axes_handles(5),'YLim',[0 150])

prettyFig('plw=2;','lw=2;','fs=12;')

if being_published
	snapnow
	delete(gcf)
end



%         ########    ###     ######  ########     ######      ###    #### ##    ## 
%         ##         ## ##   ##    ##    ##       ##    ##    ## ##    ##  ###   ## 
%         ##        ##   ##  ##          ##       ##         ##   ##   ##  ####  ## 
%         ######   ##     ##  ######     ##       ##   #### ##     ##  ##  ## ## ## 
%         ##       #########       ##    ##       ##    ##  #########  ##  ##  #### 
%         ##       ##     ## ##    ##    ##       ##    ##  ##     ##  ##  ##   ### 
%         ##       ##     ##  ######     ##        ######   ##     ## #### ##    ## 

%          ######   #######  ##    ## ######## ########   #######  ##       
%         ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%         ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
%         ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%         ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%         ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%          ######   #######  ##    ##    ##    ##     ##  #######  ######## 


fig_handle=figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
clf(fig_handle);
axes_handles(1) = subplot(10,3,1:9); 
axes_handles(2) = subplot(10,3,10:18); 
axes_handles(3) = subplot(10,3,[19:3:28]);
axes_handles(4) = subplot(10,3,[20:3:29]);
axes_handles(5) = subplot(10,3,[21:3:30]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end

movePlot(axes_handles(1),'up',.03)
movePlot(axes_handles(2),'up',.04)
movePlot(axes_handles(3),'down',.02)
movePlot(axes_handles(4),'down',.02)
movePlot(axes_handles(5),'down',.02)

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);

% A spikes --> firing rate
hash = dataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	disp('Computing firing rate for A neuron...')
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end


tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = PID(i,1:10:end);
end
PID = PID2; clear PID2


% set up a colour map
c = parula(8);

% plot stimulus
plot(axes_handles(1),tA,mean2(PID),'k')
set(axes_handles(1),'XLim',[0 60],'XTick',[])
ylabel(axes_handles(1),'Stimulus (V)')

% make linear predictions for the whole experiment
K = NaN(1e3,1);
[this_K, filtertime_full] = fitFilter2Data(mean2(PID(10e3:end,:)),mean2(fA(10e3:end,:)),'reg',1,'filter_length',1200,'offset',200);
K = this_K(100:end-101);
filtertime = filtertime_full(100:end-101);

fp = NaN*fA;
for i = 1:width(fA)
	fp(:,i) = convolve(tA,PID(:,i),K,filtertime);
end

clear l
R = mean2(fA);
[ax,plot1,plot2] = plotyy(axes_handles(2),tA,R,tA,mean2(fp));
set(ax(1),'XLim',[0 60],'YLim',[min(mean2(fA)) 2*max(mean2(fA))])
set(ax(2),'XLim',[0 60],'YLim',[min(mean2(fp)) max(mean2(fp))])
set(plot1,'Color','k')
set(plot2,'Color','r')
ylabel(ax(1),'ORN Response (Hz)')
ylabel(ax(2),'Projected Stimulus')
set(axes_handles(2),'box','off')
xlabel(axes_handles(2),'Time (s)')

% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles([4 5]);

hl_min = .1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

resp = mean2(fA(10e3:55e3,[3:10 13:20]));
pred = mean2(fp(10e3:55e3,[3:10 13:20]));
time = 1e-3*(1:length(resp));
stim = mean2(PID(10e3:55e3,[3:10 13:20]));

[p,~,~,~,~,history_lengths]=gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);
set(axes_handles(5),'XLim',[.09 11]) % to show .1 and 10 on the log scale
set(axes_handles(4),'XLim',[min(pred) max(pred)])
xlabel(axes_handles(4),'Projected Stimulus')
ylabel(axes_handles(4),'ORN Response (Hz)')

% show the p-value
axes(axes_handles(4))
text(-.2,60,'p < 0.01')
title(axes_handles(4),[])

% plot gain vs preceding stimulus
[x,y] = makeFig6G(mean2(PID),mean2(fA),mean2(fp),500);
gain_time = mean(diff(tA))*(1:length(x));
rm_this = (isnan(x) | isnan(y));
x(rm_this) = [];
y(rm_this) = [];
gain_time(rm_this) = [];
ss = 50;
plot(axes_handles(3),x(1:ss:end),y(1:ss:end)/mean(y),'k.')
xlabel(axes_handles(3),'Stimulus in preceding 500ms (V)')
ylabel(axes_handles(3),'Relative gain')
set(axes_handles(3),'YLim',[0 3.5])

% fix some labels
set(axes_handles(5),'YScale','log','YTick',[0.5 1 2],'YLim',[0.4 2.5],'XLim',[0.09 10.1],'YMinorTick','on')

% Indicate regions of high and low stimulus on the stimulus
hl = round(history_lengths(9)*1e3);
shat = computeSmoothedStimulus(mean2(PID),hl);

n = floor(sum(~isnan(mean2(fA)))*.33);
shat(1:hl) = Inf; % the initial segment where we can't estimate shat is excluded
shat(isnan(shat)) = Inf;
shat(isnan(mean2(fA))) = Inf;
[~, t_low] = sort(shat,'ascend');
t_low = t_low(1:n); % this is an index
t_low = tA(t_low); % t_low is now a time. 
 
shat = computeSmoothedStimulus(mean2(PID),hl);
shat(1:hl) = -Inf;
shat(isinf(shat)) = -Inf;
shat(isnan(mean2(fA))) = -Inf;
[~, t_high] = sort(shat,'descend');
t_high = t_high(1:n);
t_high  = tA(t_high);

plot(axes_handles(1),t_low,1+0*t_low,'.','Color',c(1,:))
plot(axes_handles(1),t_high,1+0*t_low,'.','Color',c(7,:))


prettyFig('plw=1.5;','lw=1.5;','fs=14;','FixLogX=1;')

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
disp(dataHash(strcat(mfilename,'.m'),Opt))

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

% tag the file as being published 
if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
