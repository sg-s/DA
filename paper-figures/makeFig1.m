% makeFig1.m
% makes figure 1 for the paper
% 
% created by Srinivas Gorur-Shandilya at 6:52 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
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


%% Figure 1: Gain changes with a naturalistic stimulus
% 

clearvars -except being_published 

scatter_size = 12;

% make figure placeholders 
fig_handle=figure('outerposition',[0 0 800 950],'PaperUnits','points','PaperSize',[800 950]); hold on
clf(fig_handle);
axes_handles(1) = subplot(16,6,[1:2 7:8 13:14 19:20]); 
axes_handles(2) = subplot(16,6,2+[1:2 7:8 13:14 19:20]); 
axes_handles(3) = subplot(16,6,4+[1:2 7:8 13:14 19:20]); 

axes_handles(4) = subplot(16,6,25:42); 
axes_handles(5) = subplot(16,6,43:60); 

axes_handles(6) = subplot(16,6,[61:63 67:69 73:75 79:81 85:87 91:93]);
axes_handles(7) = subplot(16,6,3+[61:63 67:69 73:75 79:81 85:87 91:93]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end

load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;


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

% remove the baseline from the PID, and remember the error
PID_baseline = mean(mean(PID(1:5e3,:)));
PID = PID - PID_baseline;


plot(axes_handles(4),tA(1:10:end),mean(PID(1:10:end,:),2),'Color',[0.2 .2 .2]);
set(axes_handles(4),'XLim',[0 70],'YLim',[0 7],'XTick',[])
ylabel(axes_handles(4),'Stimulus (V)')

% % make an inset showing details of stimulus reproducibility 
% inset(1) = axes(); % to show natural flickering vs. linear prediction
% set(inset(1),'Position',[.58 .62 .18 .08],'box','on','XTickLabel',{},'YTickLabel',{})
% plot(inset(1),PID(27e3:29e3,:),'Color',[.5 .5 .5])
% set(inset(1),'XLim',[1 1e3],'XTick',[],'YTick',[],'YLim',[0 7])

axes(axes_handles(1));
y = zeros(300,width(PID));
for i = 1:width(PID)
	[y(:,i),x] = histcounts(PID(:,i),300);x(1) = [];
	y(:,i) = y(:,i)/sum(y(:,i));
end
errorShade(x,mean(y,2),sem(y'),'Color',[.2 .2 .2]);
warning off % because there are some -ve values on the log scale
set(axes_handles(1),'XScale','log','YScale','log','XLim',[min(x) 10],'YLim',[1e-5 1],'YTick',logspace(-5,0,6))
xlabel(axes_handles(1),'Stimulus (V)')
ylabel(axes_handles(1),'Probability')
warning on


% show the whiff durations 
whiff_durations = []; 
for i = 1:width(PID)
	[ons,offs] = computeOnsOffs(PID(:,i) > .024);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1)  =[];
y = y/sum(y);
a = 1; m = fittype('a*(x).^n');
ff = fit(x(a:end)',y(a:end)',m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(axes_handles(2),x,y,'k+')
plot(axes_handles(2),x,ff(x),'r')
set(axes_handles(2),'YScale','log','XScale','log')
xlabel(axes_handles(2),'Whiff duration (ms)')


% show the blank durations 
whiff_durations = [];
for i = 1:width(PID)
	[ons,offs] = computeOnsOffs(PID(:,i) < .024);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1)  =[];
y = y/sum(y);
a = 1; m = fittype('a*(x).^n');
ff = fit(x(a:end)',y(a:end)',m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(axes_handles(3),x,y,'k+')
plot(axes_handles(3),x,ff(x),'r')
set(axes_handles(3),'YScale','log','XScale','log')
xlabel(axes_handles(3),'Blank duration (ms)')


% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

% plot the response and the prediction
clear l
[ax,plot1,plot2] = plotyy(axes_handles(5),tA,R,tA,fp);
set(ax(1),'XLim',[0 70],'YLim',[0 120])
set(ax(2),'XLim',[0 70],'YLim',[min(fp) max(fp)])
set(plot1,'Color','b')
set(plot2,'Color','r')
ylabel(ax(1),'ab3A Response (Hz)')
ylabel(ax(2),'Projected Stimulus (V)')
set(axes_handles(5),'box','off')

shat = computeSmoothedStimulus(mean(PID,2),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
axes(axes_handles(6))
cc = parula(100);
c = cc(shat,:);
ab3.c = c;
scatter(fp,R,scatter_size,c,'filled')
xlabel(axes_handles(6),'Projected Stimulus (V)')
ylabel(axes_handles(6),'ab3A response (Hz)')
shat = computeSmoothedStimulus(mean(PID,2),500);
ch = colorbar('east');
set(ch,'Position',[.44 .1 .02 .1])
caxis([min(shat) max(shat)]);

% plot gain vs stimulus for all these whiffs
axes(axes_handles(7))

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

l(1) = plot(axes_handles(7),mean_stim,gain,'k+');


% save this for a later fit
ab3.mean_stim = mean_stim;
ab3.gain = gain;
ab3.gain_err = gain_err;
ab3.R = R;
ab3.fp = fp;
ab3.K = K;
ab3.PID = mean(PID,2);

% now also add ab2 data
load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab2_1_1_all.mat')
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

% remove the baseline from the PID, and remember the error
PID_baseline = mean(mean(PID(1:5e3,:)));
PID_baseline_err = std(mean(PID(1:5e3,:)));
PID = PID - PID_baseline;

% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

shat = computeSmoothedStimulus(mean(PID,2),500);

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
rm_this = (abs(gain_err./gain)) > .5 | gain < 0; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];

l(2) = plot(axes_handles(7),mean_stim,gain,'ko');

% fit a inverse relationship to all the data
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./[gain_err; ab3.gain_err];
ff = fit([mean_stim; ab3.mean_stim],[gain; ab3.gain],'power1',options);
plot(axes_handles(7),sort([mean_stim; ab3.mean_stim]),ff(sort([mean_stim; ab3.mean_stim])),'r');

legend(l,{'ab3A','ab2A'},'Location','southwest')
xlabel(axes_handles(7),'Stimulus in preceding 500ms (V)')
ylabel(axes_handles(7),'Gain (Hz/V)')

% cosmetic fixes
movePlot(axes_handles(1),'up',.05)
movePlot(axes_handles(2),'up',.05)
movePlot(axes_handles(3),'up',.05)

movePlot(axes_handles(1),'left',.02)
movePlot(axes_handles(3),'right',.02)

movePlot(axes_handles(4),'up',.025)
movePlot(axes_handles(5),'up',.025)

movePlot(axes_handles(6),'down',.025)
movePlot(axes_handles(7),'down',.025)
movePlot(axes_handles(6),'left',.025)
movePlot(axes_handles(7),'right',.025)



set(axes_handles(2),'XLim',[50 5000])
set(axes_handles(1),'XLim',[1e-2 10])
xlabel(axes_handles(5),'Time (s)')
set(axes_handles(6),'XLim',[-.1 5])
set(axes_handles(7),'YLim',[0 max(gain)],'YScale','log','XScale','log','XLim',[.001 10])


% more cosmetics
set(ax(2),'YColor','r')
set(ax(1),'YColor','b')
set(axes_handles(6),'YColor','b','XColor','r')

prettyFig('plw=1.5;','lw=1.5;','fs=12;')

legend('boxoff')

if being_published
	snapnow
	delete(gcf)
end


%% Supplmentary Figure
% This supplementary figure shows that the choice of filter doesn't affect our results of large, rapid gain control. 

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on

% first row: ab3A filters and projections
subplot(2,3,1), hold on
plot(filtertime,ab3.K,'b');

% load MSG filter for ab3
load('../data/MSG_per_neuron.mat');
ab3.K_white = mean(reshape([MSG_data(1,:).K],1001,11),2);
t = 1e-3*(1:length(ab3.K_white))- .2;
plot(t,ab3.K_white,'r')
legend({'Nat. Stimulus Filter','White Noise filter'})
ylabel('ab3A filter')
xlabel('Lag (s)')

ab3.fp2 = convolve(1e-3*(1:length(ab3.PID)),ab3.PID,ab3.K_white,t);

subplot(2,3,3), hold on
scatter(ab3.fp2,ab3.R,scatter_size,ab3.c,'filled')
xlabel('Projected Stimulus (V)')
title('using white-noise filter')
ylabel('ab3A Response (Hz)')
set(gca,'XColor','r','XLim',[-1 14])

subplot(2,3,2), hold on
scatter(ab3.fp,ab3.R,scatter_size,ab3.c,'filled')
xlabel('Projected Stimulus (V)')
ylabel('ab3A Response (Hz)')
title('using nat. stimulus filter')
set(gca,'XColor','b','XLim',[-.5 6])

subplot(2,3,4), hold on
plot(filtertime,K,'b')
ylabel('ab2A filter')
xlabel('Lag (s)')

subplot(2,3,5), hold on

shat = computeSmoothedStimulus(mean(PID,2),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
cc = parula(100);
c= cc(shat,:);
scatter(fp,R,scatter_size,c,'filled')
xlabel('Projected Stimulus (V)')
ylabel('ab2A response (Hz)')
shat = computeSmoothedStimulus(mean(PID,2),500);
ch = colorbar('east');
set(ch,'Position',[.58 .15 .02 .1])
caxis([min(shat) max(shat)]);
set(gca,'XLim',[-.15 4],'XColor','b')

prettyFig('fs=12;')

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
