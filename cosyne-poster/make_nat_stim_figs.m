% makeFig1.m
% makes figure 1 for the paper
% 
% created by Srinivas Gorur-Shandilya at 6:52 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;


%% Figure 1: Gain changes with a naturalistic stimulus
% In this figure, we show that the gain of ab3A and ab2A ORNs changes dramatically in response to a naturalistic stimulus, and that this gain change can be correlated to the mean or the variance of the stimulus in the last 500ms. 

clearvars -except being_published 

scatter_size = 12;

% make figure placeholders 
fig_handle=figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
clf(fig_handle);

axes_handles(1) = subplot(2,1,1);
for i = 2:5
	axes_handles(i) = subplot(2,4,4+i-1);
	hold(axes_handles(i),'on');
end



load('/local-data/DA-paper/natural-flickering/without-lfp/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
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

% plot the sitmulus
plot(axes_handles(1),tA(1:10:end),mean(PID(1:10:end,:),2),'Color',[0.2 .2 .2]);
set(axes_handles(1),'XLim',[0 70],'YLim',[-.1 7])



% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);


% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);


shat = computeSmoothedStimulus(mean(PID,2),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
axes(axes_handles(2))
cc = parula(100);
c = cc(shat,:);
ab3.c = c;
scatter(fp,R,scatter_size,c,'filled')

shat = computeSmoothedStimulus(mean(PID,2),500);
ch = colorbar('east');
set(ch,'Position',[0.2182    0.1462    0.0115    0.1358])
caxis([min(shat) max(shat)]);

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


l(1) = plot(axes_handles(3),mean_stim,gain,'k+');




% save this for a later fit
ab3.mean_stim = mean_stim;
ab3.gain = gain;
ab3.gain_err = gain_err;
ab3.R = R;
ab3.fp = fp;
ab3.K = K;
ab3.PID = PID;

% now also add ab2 data
load('/local-data/DA-paper/natural-flickering/without-lfp/2014_07_11_EA_natflick_non_period_CFM_1_ab2_1_1_all.mat')
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

l(2) = plot(axes_handles(3),mean_stim,gain,'ko');

% % fit a inverse relationship to all the data
% options = fitoptions(fittype('power1'));
% options.Lower = [-Inf -1];
% options.Upper = [Inf -1];
% options.Weights = 1./[gain_err; ab3.gain_err];
% ff = fit([mean_stim; ab3.mean_stim],[gain; ab3.gain],'power1',options);
% plot(axes_handles(7),sort([mean_stim; ab3.mean_stim]),ff(sort([mean_stim; ab3.mean_stim])),'r');
legend(l,{'ab3A','ab2A'},'Location','southwest')



% show the gain as function of convolution with a differentiating filter. 
Kdiff = [ones(250,1) ; -ones(250,1)];
shat = abs(filter(Kdiff,length(Kdiff),mean(ab3.PID,2)));
% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = computeOnsOffs(ab3.R>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff = fit(fp(whiff_starts(i):whiff_ends(i)),ab3.R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5 | gain < 0; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];
clear l
l(1) = plot(axes_handles(4),mean_stim,gain,'k+');

% save this so we can fit a line to all the data
ab3.std_stim = mean_stim;
ab3.std_stim_gain = gain;
ab3.std_stim_gain_err = gain_err;

% now add the ab2 data
shat = abs(filter(Kdiff,length(Kdiff),mean(PID,2)));

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
l(2) = plot(axes_handles(4),mean_stim,gain,'ko');
legend(l,{'ab3A','ab2A'})


% now show that the variance and the mean are correlated 
all_block_sizes = factor2(length(ab3.PID));
all_block_sizes = all_block_sizes(6:end-1);
all_block_sizes = all_block_sizes(1:41);
clear l r2
r2 = NaN*all_block_sizes;

for i = 1:length(all_block_sizes)
	temp = ab3.PID(:);
	temp = reshape(temp,all_block_sizes(i),length(temp)/all_block_sizes(i));
	if all_block_sizes(i) == 5e2
		plot(axes_handles(5),mean(temp),std(temp),'Marker','.','Color',[.5 .5 .5],'LineStyle','none')
	end
	r2(i) = rsquare(mean(temp),std(temp));
end
plot(axes_handles(5),[1e-3 2],[1e-3 2],'k--')

% label all the axes, set the scales, etc.
set(axes_handles(5),'YScale','log','XScale','log','XLim',[1e-3 1e1],'XTick',[1e-3 1e-2 1e-1 1e0 1e1],'YLim',[1e-3 1e1])

set(axes_handles(1),'box','off')

set(axes_handles(2),'XLim',[-.1 5.1],'YColor','k','XColor','r','box','off')

set(axes_handles(3),'YLim',[10 1e4],'YScale','log','XScale','log','XLim',[.0009 10],'XTick',[1e-4 1e-3 1e-2 1e-1 1 10])
set(axes_handles(4),'YLim',[10 1e4],'YScale','log','XScale','log','XLim',[.0009 10],'XTick',[1e-4 1e-3 1e-2 1e-1 1 10])


prettyFig('plw=1.5;','lw=1.5;','fs=12;','FixLogX=true;')


% move plots a bit to make room for labels
movePlot(axes_handles(2),'left',.025)

movePlot(axes_handles(4),'right',.0125)
movePlot(axes_handles(5),'right',.025)

