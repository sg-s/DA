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
fig_handle=figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
clf(fig_handle);

axes_handles(2) = subplot(7,4,[1 5 9 13]);   % cartoon showing prep
axes_handles(3) = subplot(7,4,[2 3 6 7]); % stimulus 
axes_handles(4) = subplot(7,4,[10 11 14 15]); % response + linear prediction 
axes_handles(5) = subplot(7,4,[4 8]);  % linear filter

axes_handles(6) = subplot(7,4,[17:4:25]);
axes_handles(7) = subplot(7,4,1+[17:4:25]);
axes_handles(8) = subplot(7,4,2+[17:4:25]);
axes_handles(9) = subplot(7,4,3+[17:4:25]);

for i = 2:length(axes_handles)
	hold(axes_handles(i),'on');
end

axes_handles(2).Position  = [0.1300    0.4738    0.12    0.34];
o = imread('../images/prep.png');
axes(axes_handles(2));
imagesc(o);
axis ij
axis image
axis off


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
plot(axes_handles(3),tA(1:10:end),mean(PID(1:10:end,:),2),'Color',[0.2 .2 .2]);
set(axes_handles(3),'XLim',[0 70],'YLim',[0 7],'XTick',[])



% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% show the filter
plot(axes_handles(5),filtertime,K,'r');

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

% plot the response and the prediction
clear l
[ax,plot1,plot2] = plotyy(axes_handles(4),tA,R,tA,fp);
set(ax(1),'XLim',[0 70],'YLim',[0 120],'YColor','k')
set(ax(2),'XLim',[0 70],'YLim',[min(fp) max(fp)])
set(plot1,'Color','k')
set(plot2,'Color','r')
ylabel(ax(1),'ab3A Response (Hz)')
ylabel(ax(2),'Projected Stimulus (V)')
set(axes_handles(4),'box','off')

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

shat = computeSmoothedStimulus(mean(PID,2),500);
ch = colorbar('east');
set(ch,'Position',[0.2582    0.1462    0.0115    0.1358])
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


l(1) = plot(axes_handles(7),mean_stim,gain,'k+');


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

l(2) = plot(axes_handles(7),mean_stim,gain,'ko');

% fit a inverse relationship to all the data
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./[gain_err; ab3.gain_err];
ff = fit([mean_stim; ab3.mean_stim],[gain; ab3.gain],'power1',options);
plot(axes_handles(7),sort([mean_stim; ab3.mean_stim]),ff(sort([mean_stim; ab3.mean_stim])),'r');
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
l(1) = plot(axes_handles(8),mean_stim,gain,'k+');

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
l(2) = plot(axes_handles(8),mean_stim,gain,'ko');
legend(l,{'ab3A','ab2A'})

% fit a line to this
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./[ab3.std_stim_gain_err; gain_err];
ff = fit([ab3.std_stim; mean_stim],[ab3.std_stim_gain; gain],'power1',options);
plot(axes_handles(8),[1e-4 1],ff([1e-4 1]),'r')

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
		plot(axes_handles(9),mean(temp),std(temp),'Marker','.','Color',[.5 .5 .5],'LineStyle','none')
	end
	r2(i) = rsquare(mean(temp),std(temp));
end
plot(axes_handles(9),[1e-3 2],[1e-3 2],'k--')

% label all the axes, set the scales, etc.
ylabel(axes_handles(3),'Stimulus (V)')

xlabel(axes_handles(4),'Time (s)')

xlabel(axes_handles(6),'Projected Stimulus (V)')
ylabel(axes_handles(6),'ab3A response (Hz)')
set(axes_handles(6),'XLim',[-.1 5.1],'YColor','k','XColor','r','box','off')

xlabel(axes_handles(9),'\mu_{stimulus} (V)')
ylabel(axes_handles(9),'\sigma_{stimulus} (V)')
set(axes_handles(9),'XLim',[0 2],'YLim',[0 2])

xlabel(axes_handles(7),'\mu_{Stimulus} in preceding 500ms (V)')
ylabel(axes_handles(7),'Gain (Hz/V)')
set(axes_handles(7),'YLim',[10 1e4],'YScale','log','XScale','log','XLim',[.001 10],'XTick',[.001 .01 .1 1 10])

xlabel(axes_handles(8),'\sigma_{Stimulus} in preceding 500ms (V)')
set(axes_handles(8),'XTick',[1e-4 1e-3 1e-2 1e-1 1 10],'XScale','log','YScale','log','XLim',[1e-4 10],'YLim',[10 1e4])

% shrink the bottom row a little bit
for i = 6:9
	temp = get(axes_handles(i),'Position');
	temp(4) = .27;
	set(axes_handles(i),'Position',temp);
end

% fix position of filter plot
set(axes_handles(5),'Position',[.8 .65 .08 .14],'box','on')
xlabel(axes_handles(5),'Filter Lag (s)')
ylabel(axes_handles(5),'Filter')

prettyFig('plw',1.5,'lw',1.5,'fs',12,'FixLogX',true)

set(axes_handles(9),'YScale','log','XScale','log','XLim',[1e-3 1e1],'XTick',[1e-3 1e-2 1e-1 1e0 1e1],'YLim',[1e-3 1e1])
set(axes_handles(6),'box','off')

legend('boxoff')
labelFigure

if being_published
	snapnow
	delete(gcf)
end

return

%% Sanity Check 1
% Here, I show that my filter extraction works well by comparing it to Damon's FFT-based filter extraction function. Note that the two filters are almost identical:

K_Damon = backOutFilter(mean(PID,2),R,'offset',200);
figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,K/max(K),'k')
plot(filtertime(1:700),K_Damon/max(K_Damon),'r')
xlabel('Filter Lag (s)')
ylabel('Filter (a.u.)')
legend({'Srinivas Filter','Damon Filter'})

prettyFig;

legend('boxoff')

if being_published
	snapnow
	delete(gcf)
end



%      ######  ##     ## ########  ########     ######## ####  ######   
%     ##    ## ##     ## ##     ## ##     ##    ##        ##  ##    ##  
%     ##       ##     ## ##     ## ##     ##    ##        ##  ##        
%      ######  ##     ## ########  ########     ######    ##  ##   #### 
%           ## ##     ## ##        ##           ##        ##  ##    ##  
%     ##    ## ##     ## ##        ##           ##        ##  ##    ##  
%      ######   #######  ##        ##           ##       ####  ######   


%% Supplmentary Figure 1
% This figure shows some statistics of the stimulus. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on

% show the rsquare of the mean and variance as a function of box size
subplot(2,2,1), hold on
for i = 1:length(all_block_sizes)
	plot(all_block_sizes(i),r2(i),'k+')
end
set(gca,'XScale','log','XTick',[1 1e1 1e2 1e3 1e4])
xlabel('Window (ms)')
ylabel('r^2 (\mu, \sigma)')

% show the PDF of the stimulus
subplot(2,2,2), hold on
y = zeros(300,width(PID));
for i = 1:width(PID)
	[y(:,i),x] = histcounts(PID(:,i),300);x(1) = [];
	y(:,i) = y(:,i)/sum(y(:,i));
end
errorShade(x,mean(y,2),sem(y'),'Color',[.2 .2 .2]);
warning off % because there are some -ve values on the log scale
set(gca,'XScale','log','YScale','log','XLim',[min(x) 10],'YLim',[1e-5 1],'YTick',logspace(-5,0,6))
xlabel(gca,'Stimulus (V)')
ylabel(gca,'Probability')
warning on

% show the whiff durations 
subplot(2,2,3), hold on
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
plot(x,y,'k+')
plot(x,ff(x),'r')
ylabel(gca,'Probability')
set(gca,'YScale','log','XScale','log')
xlabel('Whiff duration (ms)')

% show the blank durations 
subplot(2,2,4), hold on
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
plot(x,y,'k+')
plot(x,ff(x),'r')
set(gca,'YScale','log','XScale','log')
xlabel('Blank duration (ms)')
ylabel(gca,'Probability')

prettyFig('fs',18,'FixLogX',true)
labelFigure

legend('boxoff')

if being_published
	snapnow
	delete(gcf)
end




%% Supplementary Figure 2
% This supplementary figure shows that the choice of filter doesn't affect our results of large, rapid gain control. 


figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on

% ab3A filters and projections
subplot(2,4,1), hold on
plot(filtertime,ab3.K,'b');

% load MSG filter for ab3
load('../data/MSG_per_neuron.mat');
ab3.K_white = mean(reshape([MSG_data(1,:).K],1001,11),2);
t = 1e-3*(1:length(ab3.K_white))- .2;
plot(t,ab3.K_white,'r')
legend({'Nat. Stimulus Filter','White Noise filter'})
ylabel('ab3A filter')
xlabel('Lag (s)')
set(gca,'YLim',[-.005 .04])

ab3.fp2 = convolve(1e-3*(1:length(mean(ab3.PID,2))),mean(ab3.PID,2),ab3.K_white,t);

subplot(2,4,3), hold on
scatter(ab3.fp2,ab3.R,scatter_size,ab3.c,'filled')
xlabel('Projected Stimulus (V)')
title('using white-noise filter')
ylabel('ab3A Response (Hz)')
set(gca,'XColor','r','XLim',[-1 14])

subplot(2,4,2), hold on
scatter(ab3.fp,ab3.R,scatter_size,ab3.c,'filled')
xlabel('Projected Stimulus (V)')
ylabel('ab3A Response (Hz)')
title('using nat. stimulus filter')
set(gca,'XColor','b','XLim',[-.5 6])

% now also show the gain-per-whiff vs. the mean stimulus in the preceding 500ms for the white-noise-projection
shat = computeSmoothedStimulus(mean(ab3.PID,2),500);
[whiff_starts,whiff_ends] = computeOnsOffs(ab3.R>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	[ff,gof]=fit(ab3.fp2(whiff_starts(i):whiff_ends(i)),ab3.R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = gof.rsquare;
end
rm_this = gain < 0 | gain_err < .8;
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];
subplot(2,4,4), hold on
plot(mean_stim,gain,'ro');
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
ff = fit(mean_stim,gain,'power1',options);
plot(gca,mean_stim,ff(mean_stim),'r')
set(gca,'XScale','log','YScale','log')
ylabel('Gain (Hz/V)')
xlabel('Mean Stimulus in preceding 500ms')
title('Gain estimation using white-noise filter')

% now show ab2A filters
subplot(2,4,5), hold on
plot(filtertime,K,'b')
ylabel('ab2A filter')
xlabel('Lag (s)')

subplot(2,4,6), hold on
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

prettyFig('fs',18,'FixLogX',true)
labelFigure

if being_published
	snapnow
	delete(gcf)
end

%% Inst. Gain Analysis of Natural Stimulus
% Now, we perform an instantaneous gain analysis of this data, using the analysis methods we developed for looking at the fast gain control. This analysis knows nothing of "whiffs", etc, and blindly plots the inst. gain vs. the mean stimulus in this data. In the following figure, we compute the inst. gain vs. the projected stimulus, and observe that there is a minimum in the Spearman correlation at a few hundred ms. 

load('/Users/sigbhu/code/da/data/nat_stim_parametric_fits.mat','od')

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
hl = [50 300 1e4];
for i = 1:3
	plot_handles = plot(od,[ax(i) ax(4)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i),'nbins',300);
	title(ax(i),['\tau_{H} = ' oval(hl(i)) 'ms'])
	if i < 3
		delete(plot_handles(2).f2)
	end
end
set(ax(1:3),'XScale','log','YScale','log')
set(ax(4),'YLim',[-1 0])

prettyFig('fs',14);
labelFigure

if being_published
	snapnow
	delete(gcf)
end

%% Gain control is truly dynamic, and cannot be explained by a nonlinearity. 
% The big question now is if the gain control we observe merely a consequence of a static nonlinearity. Based on our visualization of the response vs. the projected stimulus, we think not: we see multiple curves that cannot be fit by any one nonlinearity. However, we also observe that the Spearman correlation in the previous plot is non-zero even at very small history lengths, suggesting that the static nonlinearity contributes to some degree to the observed gain changes. 

%% 
% To show that there is a gain control mechanism even if we account for the nonlinearity, we fit a nonlinear function to the data and then compute an instantaneous gain vs. the prediction of the full LN model. 

od = computeInstGain(od,true);

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
hl = [50 300 1e4];
for i = 1:3
	plot_handles = plot(od,[ax(i) ax(4)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i),'nbins',100);
	title(ax(i),['\tau_{H} = ' oval(hl(i)) 'ms'])
	if i < 3
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (norm)')
end
set(ax(1:3),'XScale','log','YScale','log','YLim',[.5 10])
set(ax(4),'YLim',[-1 .5])

prettyFig('fs',14,'FixLogY',true);

if being_published
	snapnow
	delete(gcf)
end

%%
% In the following figure, we plot the inst. gain vs the stimulus in the preceding window for a variety of windows, to verify that what we have is actually a change in these clouds of points. 

figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
for i = 1:8
	ax(i) = subplot(2,4,i); hold on
end
hl = round(logspace(1,4,8));
dm = figure; hold on
dummy = subplot(2,2,1);

for i = 1:length(hl)
	plot_handles = plot(od,[ax(i) dummy],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i),'data_bin_type','dots','nbins',5);
	title(ax(i),['\tau_{H} = ' oval(hl(i)) 'ms'])
	if i < 3
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (norm)')
end
set(ax,'XScale','log','YScale','log','YLim',[.1 100],'XLim',[1e-2 10])

delete(dm)

prettyFig('fs',14,'FixLogY',true,'FixLogX',true);
labelFigure

if being_published
	snapnow
	delete(gcf)
end


%% Dynamic Gain Control
% We have shown through our inst. gain aanlysis that there is a gain control mechanism that is "dynamic" in the sense that it cannot be explained by a static nonlinearity. If this is true, a dynamic gain model (the DA Model) should outperform a LN model in explaining this data. Let's check. 

% fit LN model
R = nanmean(od.firing_rate,2);
S = nanmean(od.stimulus,2);
LN = nanmean(od.firing_projected,2);

temp1 = LN(:); temp2 = R(:);
rm_this = isnan(temp1) | isnan(temp2);
temp1(rm_this) = []; temp2(rm_this) = [];
ft = fittype('hillFit(x,A,k,n,x_offset)');
ff = fit(temp1,temp2,ft,'StartPoint',[max(temp2) mean(temp1) 1 0],'Upper',[1e3 max(temp1) 10 Inf],'Lower',[min(temp2)/2 0 0 -Inf]);
LN = ff(LN);

% fit DA model
clear p
p.   s0 = 7.8242e-04;
p.  n_z = 2;
p.tau_z = 151.1249;
p.  n_y = 2;
p.tau_y = 26.7002;
p.    C = 0.5457;
p.    A = 163.2252;
p.    B = 2.4703;
DA = DAModelv2(S,p);

figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
subplot(2,1,1), hold on
plot(tA,R,'k')
l = plot(tA,LN,'r');
legend(l,['LN Model, r^2 = ' oval(rsquare(LN,R))]);
ylabel('Firing Rate (Hz)')

subplot(2,1,2), hold on
plot(tA,R,'k')
l = plot(tA,DA,'r');
legend(l,['LN Model, r^2 = ' oval(rsquare(DA,R))]);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

prettyFig('fs',20,'FixLogY',true);
labelFigure

if being_published
	snapnow
	delete(gcf)
end

%%
% Finally, we can check that the timescale of gain control as reported by the best-fit DA model agrees with the timescale we determined using the fast gain control. The timescale of fast gain control in the DA model is (in ms):

disp(p.tau_z*p.n_z)

%% Visualizing the changing gain 
% In this section, we try to visualize how the input-output curve changes over the course of this experiment by binning the data into regions where the stimulus is high or low and visualizing input-output curves for those regions. In the following figure, we plot the input-output curves for the data segregated into five quintiles. Brighter colours indicate higher stimulus in the preceding 500ms.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(od,gca,'dynamicIO.firing_rate','nbins',19,'min_inst_gain_firing',2)

subplot(1,2,2), hold on
da_model = ORNData;
da_model.stimulus = S;
da_model.firing_projected = nanmean(od.firing_projected,2);
da_model.firing_rate = R;

plot(da_model,gca,'dynamicIO.firing_rate','nbins',19,'min_inst_gain_firing',2)
ylabel('DA Model Prediction (Hz)')
prettyFig('fs',20,'FixLogY',true);
labelFigure

if being_published
	snapnow
	delete(gcf)
end





%% Version Info
% The file that generated this document is called:
pFooter;