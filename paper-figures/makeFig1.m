% makeFig1.m
% makes figure 1 for the paper
% 
% created by Srinivas Gorur-Shandilya at 6:52 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
    path1 = [path1 ':/usr/local/bin'];
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


%% Figure 1

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
% scatter(fp(1:ss:end),R(1:ss:end),[],'k','filled')
scatter(fp(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')
xlabel(axes_handles(4),'Projected Stimulus')
ylabel(axes_handles(4),'Actual response (Hz)')
shat = computeSmoothedStimulus(mean2(PID),500);
colorbar;
caxis([min(shat) max(shat)]);

% plot gain vs stimulus for all these whiffs
axes(axes_handles(5))

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
