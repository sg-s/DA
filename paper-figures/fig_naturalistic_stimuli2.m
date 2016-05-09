% makeFig1.m
% makes figure 1 for the paper
% 
% created by Srinivas Gorur-Shandilya at 6:52 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

% this script uses dataManager to ensure data integrity
dm = dataManager;


%% Figure 1: Gain changes with a naturalistic stimulus
% In this figure, we show that the gain of ab3A and ab2A ORNs changes dramatically in response to a naturalistic stimulus, and that this gain change can be correlated to the mean or the variance of the stimulus in the last 500ms. 

clearvars -except being_published dm

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
axes(axes_handles(2));
imagesc(imread('../images/fig-1-cartoon.png'));
axis ij
axis image
axis off

% first, grab the ab3 and ab2 data

%         ######   ######## ########       ###    ########   #######     ###    
%        ##    ##  ##          ##         ## ##   ##     ## ##     ##   ## ##   
%        ##        ##          ##        ##   ##  ##     ##        ##  ##   ##  
%        ##   #### ######      ##       ##     ## ########   #######  ##     ## 
%        ##    ##  ##          ##       ######### ##     ##        ## ######### 
%        ##    ##  ##          ##       ##     ## ##     ## ##     ## ##     ## 
%         ######   ########    ##       ##     ## ########   #######  ##     ## 

clear ab3 ab2
load(dm.getPath('5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);

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

% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

% add to ab3 and remove everything else
ab3.PID = PID;
ab3.fA = fA;
ab3.K = K;
ab3.fp = fp;
clear PID fA K R filtertime_full PID_baseline all_spikes data fp spikes time 


% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = computeOnsOffs(mean(ab3.fA,2)>10);

% filter the stimulus using a box filter
shat_mean = computeSmoothedStimulus(mean(ab3.PID,2),500);


% filter the stimulus using a diff. filter
Kdiff = [ones(250,1) ; -ones(250,1)];
shat_std = abs(filter(Kdiff,length(Kdiff),mean(ab3.PID,2)));


% for each excursion, estimate mean, std stimulus, gain, etc. 
mean_stim = NaN*whiff_ends;
std_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
R = mean(ab3.fA,2);
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat_mean(whiff_starts(i):whiff_ends(i)));
	std_stim(i) = mean(shat_std(whiff_starts(i):whiff_ends(i)));
	ff=fit(ab3.fp(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
clear R
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];
std_stim(rm_this) = [];

% save this for a later fit
ab3.mean_stim = mean_stim;
ab3.std_stim = std_stim;
ab3.gain = gain;
ab3.gain_err = gain_err;

clear ff gain gain_err Kdiff mean_stim rm_this shat_mean shat_std std_stim temp whiff_starts whiff_ends

%  ######   ######## ########       ###    ########   #######     ###    
% ##    ##  ##          ##         ## ##   ##     ## ##     ##   ## ##   
% ##        ##          ##        ##   ##  ##     ##        ##  ##   ##  
% ##   #### ######      ##       ##     ## ########   #######  ##     ## 
% ##    ##  ##          ##       ######### ##     ## ##        ######### 
% ##    ##  ##          ##       ##     ## ##     ## ##        ##     ## 
%  ######   ########    ##       ##     ## ########  ######### ##     ## 


% now also add ab2 data
load(dm.getPath('8af556aa49c4af116c7f66e8417c0dc2'),'data','spikes')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);

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

% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

% add to ab2 and remove everything else
ab2.PID = PID;
ab2.fA = fA;
ab2.K = K;
ab2.fp = fp;
clear PID fA K R filtertime_full PID_baseline all_spikes data fp spikes time 


% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = computeOnsOffs(mean(ab2.fA,2)>10);

% filter the stimulus using a box filter
shat_mean = computeSmoothedStimulus(mean(ab2.PID,2),500);


% filter the stimulus using a diff. filter
Kdiff = [ones(250,1) ; -ones(250,1)];
shat_std = abs(filter(Kdiff,length(Kdiff),mean(ab2.PID,2)));


% for each excursion, estimate mean, std stimulus, gain, etc. 
mean_stim = NaN*whiff_ends;
std_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
R = mean(ab2.fA,2);
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat_mean(whiff_starts(i):whiff_ends(i)));
	std_stim(i) = mean(shat_std(whiff_starts(i):whiff_ends(i)));
	ff=fit(ab2.fp(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
clear R
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];
std_stim(rm_this) = [];

% save this for a later fit
ab2.mean_stim = mean_stim;
ab2.std_stim = std_stim;
ab2.gain = gain;
ab2.gain_err = gain_err;

clear ff gain gain_err Kdiff mean_stim rm_this shat_mean shat_std std_stim temp whiff_starts whiff_ends

% ##     ##    ###    ##    ## ########    ########  ##        #######  ######## 
% ###   ###   ## ##   ##   ##  ##          ##     ## ##       ##     ##    ##    
% #### ####  ##   ##  ##  ##   ##          ##     ## ##       ##     ##    ##    
% ## ### ## ##     ## #####    ######      ########  ##       ##     ##    ##    
% ##     ## ######### ##  ##   ##          ##        ##       ##     ##    ##    
% ##     ## ##     ## ##   ##  ##          ##        ##       ##     ##    ##    
% ##     ## ##     ## ##    ## ########    ##        ########  #######     ##    

% plot the stimulus
plot(axes_handles(3),tA(1:10:end),mean(ab3.PID(1:10:end,:),2),'Color',[0.2 .2 .2]);
set(axes_handles(3),'XLim',[0 70],'YLim',[0 7],'XTick',[])

% show the filter
plot(axes_handles(5),filtertime,ab3.K,'r');

% plot the response and the prediction
clear l
[ax,plot1,plot2] = plotyy(axes_handles(4),tA,mean(ab3.fA,2),tA,ab3.fp);
set(ax(1),'XLim',[0 70],'YLim',[0 150],'YColor','k')
set(ax(2),'XLim',[0 70],'YLim',[0 6],'YColor','r','YTick',[0 2 4 6])
set(plot1,'Color','k')
set(plot2,'Color','r')

set(axes_handles(4),'box','off')

shat = computeSmoothedStimulus(mean(ab3.PID,2),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
axes(axes_handles(6))
cc = parula(100);
c = cc(shat,:);
ab3.c = c;
scatter(ab3.fp,mean(ab3.fA,2),scatter_size,c,'filled')

shat = computeSmoothedStimulus(mean(ab3.PID,2),500);
ch = colorbar('east');
set(ch,'Position',[0.2582    0.1462    0.0115    0.1358])
caxis([min(shat) max(shat)]);


% plot whiff gain vs. mean stimulus preceding whiff
l(1) = plot(axes_handles(7),ab3.mean_stim,ab3.gain,'k+');
l(2) = plot(axes_handles(7),ab2.mean_stim,ab2.gain,'ko');
legend(l,{'ab3A','ab2A'},'Location','southwest')
set(axes_handles(7),'XScale','log','YScale','log')

% fit a inverse relationship to all the data
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./[ab2.gain_err; ab3.gain_err];
ff = fit([ab2.mean_stim; ab3.mean_stim],[ab2.gain; ab3.gain],'power1',options);
plot(axes_handles(7),sort([ab2.mean_stim; ab3.mean_stim]),ff(sort([ab2.mean_stim; ab3.mean_stim])),'r');

% plot whiff gain vs. std. dev. stimulus preceding whiff
l(1) = plot(axes_handles(8),ab3.std_stim,ab3.gain,'k+');
l(2) = plot(axes_handles(8),ab2.std_stim,ab2.gain,'ko');
set(axes_handles(8),'XScale','log','YScale','log')
legend(l,{'ab3A','ab2A'},'Location','southwest')

% fit a inverse relationship to all the data
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
ff = fit([ab2.std_stim; ab3.std_stim],[ab2.gain; ab3.gain],'power1',options);
ff.a = 10;
plot(axes_handles(8),sort([ab2.std_stim; ab3.std_stim]),ff(sort([ab2.std_stim; ab3.std_stim])),'r');

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

ylabel(ax(1),'ab3A Firing Rate (Hz)')
ylabel(ax(2),'Projected Stimulus (V)')

ylabel(axes_handles(3),'Stimulus (V)')

xlabel(axes_handles(4),'Time (s)')


xlabel(axes_handles(6),'Projected Stimulus (V)')
ylabel(axes_handles(6),'ab3A Firing Rate (Hz)')
set(axes_handles(6),'XLim',[-.1 5.1],'YColor','k','XColor','r','box','off')



xlabel(axes_handles(7),'\mu_{Stimulus} in preceding 500ms (V)')
ylabel(axes_handles(7),'Firing Gain (Hz/V)')
set(axes_handles(7),'YLim',[10 1e4],'YScale','log','XScale','log','XLim',[.001 10],'XTick',[.001 .01 .1 1 10])

xlabel(axes_handles(8),'\sigma_{Stimulus} in preceding 500ms (V)')
set(axes_handles(8),'XTick',[1e-4 1e-3 1e-2 1e-1 1 10],'XScale','log','YScale','log','XLim',[1e-3 10],'YLim',[10 1e4])

% shrink the bottom row a little bit
for i = 6:9
	temp = get(axes_handles(i),'Position');
	temp(4) = .27;
	set(axes_handles(i),'Position',temp);
end

% fix position of some plots
set(axes_handles(5),'Position',[.8 .65 .08 .14],'box','on')
xlabel(axes_handles(5),'Filter Lag (s)')
ylabel(axes_handles(5),'Filter')

set(axes_handles(6),'box','off')
axes_handles(6).Position(1) = .1;
axes_handles(9).Position(1) = .8;
axes_handles(8).Position(1) = .57;
axes_handles(7).Position(1) = .35;

prettyFig('plw',1.5,'lw',1.5,'fs',18,'FixLogX',true)


set(axes_handles(6),'XTick',[0 2 4 6],'box','off')

xlabel(axes_handles(9),'\mu_{stimulus} (V)')
ylabel(axes_handles(9),'\sigma_{stimulus} (V)')
set(axes_handles(9),'YScale','log','XScale','log','XLim',[1e-3 1e1],'XTick',[1e-3 1e-2 1e-1 1e0 1e1],'YLim',[1e-3 1e1],'YTick',[1e-3 1e-2 1e-1 1 10])


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
y = zeros(300,width(ab3.PID));
for i = 1:width(ab3.PID)
	[y(:,i),x] = histcounts(ab3.PID(:,i),300);x(1) = [];
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
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) > .024);
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
set(gca,'YScale','log','XScale','log','XTick',[1e1 1e2 1e3 1e4])
xlabel('Whiff duration (ms)')

% show the blank durations 
subplot(2,2,4), hold on
whiff_durations = [];
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) < .024);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1)  =[];
y = y/sum(y);
a = 1; m = fittype('a*(x).^n');
ff = fit(x(a:end)',y(a:end)',m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(x,y,'k+')
plot(x,ff(x),'r')
set(gca,'YScale','log','XScale','log','XTick',[1e2 1e3 1e4 1e5])
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
load(dm.getPath('b934ad74a78cb20df400b67e107037e3'));
ab3.K_white = mean(reshape([MSG_data(1,:).K],1001,11),2);
t = 1e-3*(1:length(ab3.K_white))- .2;
plot(t,ab3.K_white,'r')
legend({'Nat. Stimulus Filter','White Noise filter'})
ylabel('ab3A filter')
xlabel('Lag (s)')
set(gca,'YLim',[-.005 .04])

ab3.fp2 = convolve(1e-3*(1:length(mean(ab3.PID,2))),mean(ab3.PID,2),ab3.K_white,t);

subplot(2,4,3), hold on
scatter(ab3.fp2,mean(ab3.fA,2),scatter_size,ab3.c,'filled')
xlabel('Projected Stimulus (V)')
title('using white-noise filter')
ylabel('ab3A Response (Hz)')
set(gca,'XColor','r','XLim',[-1 14])

subplot(2,4,2), hold on
scatter(ab3.fp,mean(ab3.fA,2),scatter_size,ab3.c,'filled')
xlabel('Projected Stimulus (V)')
ylabel('ab3A Response (Hz)')
title('using nat. stimulus filter')
set(gca,'XColor','b','XLim',[-.5 6])

% now also show the gain-per-whiff vs. the mean stimulus in the preceding 500ms for the white-noise-projection
shat = computeSmoothedStimulus(mean(ab3.PID,2),500);
[whiff_starts,whiff_ends] = computeOnsOffs(mean(ab3.fA,2)>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
R = mean(ab3.fA,2);
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	[ff,gof]=fit(ab3.fp2(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
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
ff.a = 6;
plot(gca,mean_stim,ff(mean_stim),'r')
set(gca,'XScale','log','YScale','log','XTick',[1e-2 1e-1 1e0 1e1])
ylabel('ab3A Firing Gain (Hz/V)')
xlabel('Mean Stimulus in preceding 500ms')
title(['Gain estimation using' char(10) 'white-noise filter'])

% now show ab2A filters
subplot(2,4,5), hold on
plot(filtertime,ab2.K,'b')
ylabel('ab2A filter')
xlabel('Lag (s)')

subplot(2,4,6), hold on
shat = computeSmoothedStimulus(mean(ab2.PID,2),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
cc = parula(100);
c= cc(shat,:);
scatter(ab2.fp,mean(ab2.fA,2),scatter_size,c,'filled')
xlabel('Projected Stimulus (V)')
ylabel('ab2A response (Hz)')
shat = computeSmoothedStimulus(mean(ab2.PID,2),500);
ch = colorbar('east');
set(ch,'Position',[.58 .15 .02 .1])
caxis([min(shat) max(shat)]);
set(gca,'XLim',[-.15 4],'XColor','b')

prettyFig('fs',18,'FixLogX',true)

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
% The file that generated this document is called:
pFooter;