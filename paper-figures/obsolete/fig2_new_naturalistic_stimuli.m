% makeFig2.m
% makes figure 2 for the paper
% 
% created by Srinivas Gorur-Shandilya at 6:52 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

% this script uses dataManager to ensure data integrity

history_length = 200; % ms


%% Figure 1: Gain changes with a naturalistic stimulus
% In this figure, we show that the gain of ab3A and ab2A ORNs changes dramatically in response to a naturalistic stimulus, and that this gain change can be correlated to the mean or the variance of the stimulus in the last 200ms. 

clearvars -except being_published dm history_length od whiffs

scatter_size = 12;

% make figure placeholders 
fig_handle = figure('PaperUnits','centimeters','PaperSize',[20 12],'Position',[100 100 1200 720],'toolbar','figure'); hold on
clf(fig_handle);

axes_handles(3) = subplot(7,4,[2 3 4 6 7 8]); % stimulus 
axes_handles(4) = subplot(7,4,[10 11 12 14 15 16]); % response + linear prediction 
axes_handles(5) = subplot(7,4,[9 13]);  % linear filter

axes_handles(6) = subplot(7,4,[17:4:25]);
axes_handles(7) = subplot(7,4,1+[17:4:25]);
axes_handles(8) = subplot(7,4,2+[17:4:25]);
axes_handles(9) = subplot(7,4,3+[17:4:25]);

for i = 3:length(axes_handles)
	hold(axes_handles(i),'on');
end

% analyse kinetics of LFP and firing rate during the naturalistic stimulus presentation
load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat');
example_orn = 4;

% remove the baseline from all the stimulus traces
for i = 1:length(od)
	for j = 1:od(i).n_trials
		od(i).stimulus(:,j) = od(i).stimulus(:,j) - nanmean(od(i).stimulus(1:5e3,j));
	end
end

% compute gains, mean, std. at each whiff
if ~exist('whiffs','var')
	use_this_segment = false(length(od(1).firing_rate),1);
	use_this_segment(5e3:end-5e3) = true;
	for i = [2 3 5 6]
		temp = od(i);
		pred = nanmean(temp.firing_projected,2); pred = pred(use_this_segment);
		resp = nanmean(temp.firing_rate,2);  resp = resp(use_this_segment);
		stim = nanmean(temp.stimulus,2); stim = stim - mean(stim(1:5e3));
		stim = stim(use_this_segment);

		% find when the valve opens
		[ons,offs] = findWhiffs(stim);

		% plot the gain in each of these windows
		[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

		rm_this = gain < 0 | gain_err < .8;
		gain(rm_this) = [];
		gain_err(rm_this) = [];
		ons(rm_this) = [];
		offs(rm_this) = [];

		% find the mean stimulus in the preceding X ms in these windows
		shat = computeSmoothedStimulus(stim,history_length);
		mu = findMeanInWindows(ons,offs,shat);
		std_devs = NaN*mu;
		for j = 1:length(ons)
			std_devs(j) = std(shat(ons(j):offs(j)));
		end

		% combine
		whiffs(i).mean_s = mu(:);
		whiffs(i).gain = gain(:);
		whiffs(i).gain_err = gain_err(:);
		whiffs(i).std_s = std_devs(:);
	end
end

fp = nanmean(od(example_orn).firing_projected,2);
fA = nanmean(od(example_orn).firing_rate,2);
S = nanmean(od(example_orn).stimulus,2);

time = 1e-3*(1:length(od(1).stimulus));

% ##     ##    ###    ##    ## ########    ########  ##        #######  ######## 
% ###   ###   ## ##   ##   ##  ##          ##     ## ##       ##     ##    ##    
% #### ####  ##   ##  ##  ##   ##          ##     ## ##       ##     ##    ##    
% ## ### ## ##     ## #####    ######      ########  ##       ##     ##    ##    
% ##     ## ######### ##  ##   ##          ##        ##       ##     ##    ##    
% ##     ## ##     ## ##   ##  ##          ##        ##       ##     ##    ##    
% ##     ## ##     ## ##    ## ########    ##        ########  #######     ##    

% plot the stimulus
plot(axes_handles(3),time,S,'k')
set(axes_handles(3),'XLim',[0 70],'YLim',[0 2],'XTick',[])

% show the filter
plot(axes_handles(5),od(example_orn).filtertime_firing,nanmean([od.K_firing],2),'r');

% plot the response and the prediction
clear l
[ax,plot1,plot2] = plotyy(axes_handles(4),time,fA,time,fp);
set(ax(1),'XLim',[0 70],'YLim',[0 150],'YColor','k')
set(ax(2),'XLim',[0 70],'YLim',[0 1.5],'YColor','r')
set(plot1,'Color','k')
set(plot2,'Color','r')

set(axes_handles(4),'box','off')

shat = computeSmoothedStimulus(S,history_length);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;



% make the output analysis plot
axes(axes_handles(6))
cc = parula(100);
c = cc(shat,:);
ab3.c = c;
scatter(fp,mean(fA,2),scatter_size,c,'filled')

shat = computeSmoothedStimulus(S,history_length);
ch = colorbar('east');
set(ch,'Position',[0.2582    0.13    0.0115    0.12])
caxis([min(shat) max(shat)]);

% plot whiff gain vs. mean stimulus preceding whiff
plot(axes_handles(7),vertcat(whiffs.mean_s),vertcat(whiffs.gain),'k+');
set(axes_handles(7),'XScale','log','YScale','log','YLim',[1e1 1e3],'XLim',[1e-2 1])

% fit a inverse relationship to all the data
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./vertcat(whiffs.gain_err);
ff = fit(vertcat(whiffs.mean_s),vertcat(whiffs.gain),'power1',options);
plot(axes_handles(7),[1e-2 1],ff([1e-2 1]),'r');


% plot whiff gain vs. std. dev. stimulus preceding whiff
plot(axes_handles(8),vertcat(whiffs.std_s),vertcat(whiffs.gain),'k+');
set(axes_handles(8),'XScale','log','YScale','log')

return

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
	if all_block_sizes(i) == history_length
		plot(axes_handles(9),mean(temp),std(temp),'Marker','.','Color',[.5 .5 .5],'LineStyle','none')
	end
	r2(i) = rsquare(mean(temp),std(temp));
end
plot(axes_handles(9),[1e-3 2],[1e-3 2],'k--')




% label all the axes, set the scales, etc.

ylabel(ax(1),'ab3A firing rate (Hz)')
ylabel(ax(2),'Projected stimulus (V)')

ylabel(axes_handles(3),'Stimulus (V)')

xlabel(axes_handles(4),'Time (s)')


xlabel(axes_handles(6),'Projected stimulus (V)')
ylabel(axes_handles(6),'ab3A firing rate (Hz)')
set(axes_handles(6),'XLim',[-.1 5.1],'YColor','k','XColor','r','box','off')



xlabel(axes_handles(7),['\mu_{Stimulus} in preceding ' oval(history_length) 'ms (V)'])
ylabel(axes_handles(7),'ORN gain (Hz/V)')
set(axes_handles(7),'YLim',[10 1e4],'YScale','log','XScale','log','XLim',[.001 10],'XTick',[.001 .01 .1 1 10])

xlabel(axes_handles(8),['\sigma_{Stimulus} in preceding ' oval(history_length) 'ms (V)'])
set(axes_handles(8),'XTick',[1e-4 1e-3 1e-2 1e-1 1 10],'XScale','log','YScale','log','XLim',[1e-3 10],'YLim',[10 1e4])

% shrink the bottom row a little bit
for i = 6:9
	temp = get(axes_handles(i),'Position');
	temp(4) = .27;
	set(axes_handles(i),'Position',temp);
end

% fix position of some plots
axes_handles(5).Position(1) = .1;
xlabel(axes_handles(5),'Filter lag (s)')
ylabel(axes_handles(5),'Filter')

set(axes_handles(6),'box','off')
axes_handles(6).Position(1) = .1;
axes_handles(9).Position(1) = .8;
axes_handles(8).Position(1) = .57;
axes_handles(7).Position(1) = .35;

set(axes_handles(6),'XTick',[0 2 4 6],'box','off')

xlabel(axes_handles(9),'\mu_{stimulus} (V)')
ylabel(axes_handles(9),'\sigma_{stimulus} (V)')
set(axes_handles(9),'YScale','log','XScale','log','XLim',[1e-3 1e1],'XTick',[1e-3 1e-2 1e-1 1e0 1e1],'YLim',[1e-3 1e1],'YTick',[1e-3 1e-2 1e-1 1 10])

% add some annotation
% a(1) = annotation('arrow','Position',[0.6194 0.5532 0.0090 -0.0285]);
% a(2) = annotation('arrow','Position',[0.6278 0.6237 0.0090 -0.0255]);

% add insets showing gain control
inset(1) = axes('Position',[0.6774 0.7907 0.1012 0.1040]); hold on
inset(2) = axes('Position',[0.6774 0.56 0.1012 0.1040]); hold on
S = mean(ab3.PID,2);
R = mean(ab3.fA,2);
t = 1e-3*(1:1000)-.8;
plot(inset(1),t,S(6042-800:6042+199))
plot(inset(1),t,S(24990-800:24990+199))
plot(inset(2),t,R(6042-800:6042+199))
plot(inset(2),t,R(24990-800:24990+199))

prettyFig('plw',1.5,'lw',1.5,'fs',.5,'FixLogX',true,'font_units','centimeters')



% compress some plots to make small EPS files
shrinkDataInPlot(axes_handles(3),1);
shrinkDataInPlot(axes_handles(4),2);
shrinkDataInPlot(ax(2),1);
shrinkDataInPlot(axes_handles(5),2)
shrinkDataInPlot(axes_handles(6),1)

axes_handles(9).XLim(2) = 11;
deintersectAxes(axes_handles(9))

axes_handles(8).XLim(2) = 11;
deintersectAxes(axes_handles(8))

axes_handles(7).XLim(2) = 11;
deintersectAxes(axes_handles(7))
deintersectAxes(axes_handles(5))

if being_published
	snapnow
	delete(gcf)
end

return

%      ######  ##     ## ########  ########     ######## ####  ######   
%     ##    ## ##     ## ##     ## ##     ##    ##        ##  ##    ##  
%     ##       ##     ## ##     ## ##     ##    ##        ##  ##        
%      ######  ##     ## ########  ########     ######    ##  ##   #### 
%           ## ##     ## ##        ##           ##        ##  ##    ##  
%     ##    ## ##     ## ##        ##           ##        ##  ##    ##  
%      ######   #######  ##        ##           ##       ####  ######   


%% Supplementary Figure 1
% This figure shows some statistics of the stimulus. 


clear ax
figure('PaperUnits','centimeters','PaperSize',[20 5],'Position',[100 100 1300 300],'toolbar','none'); hold on

% show the rsquare of the mean and variance as a function of box size
ax(1) = subplot(1,4,1); hold on
for i = 1:length(all_block_sizes)
	plot(ax(1),all_block_sizes(i),r2(i),'k+')
end
set(ax(1),'XScale','log','XTick',[1 1e1 1e2 1e3 1e4],'XLim',[1 1.1e4])
xlabel(ax(1),'Window (ms)')
ylabel(ax(1),'r^2 (\mu, \sigma)')

% show the PDF of the stimulus
ax(2) = subplot(1,4,2); hold on
y = zeros(300,width(ab3.PID));
for i = 1:width(ab3.PID)
	[y(:,i),x] = histcounts(ab3.PID(:,i),300);x(1) = [];
	y(:,i) = y(:,i)/sum(y(:,i));
end
errorShade(x,mean(y,2),sem(y'),'Color',[.2 .2 .2]);
warning off % because there are some -ve values on the log scale
set(ax(2),'XScale','log','YScale','log','XLim',[min(x) 11],'YLim',[1e-5 1],'YTick',logspace(-5,0,6),'XTick',[1e-2 1e-1 1 10])
xlabel(ax(2),'Stimulus (V)')
ylabel(ax(2),'Probability')
warning on

odour_thresh = 0.024;

% show the whiff durations 
ax(3) = subplot(1,4,3); hold on
whiff_durations = []; 
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) > odour_thresh);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1) = [];
y = y/sum(y);
a = 1; m = fittype('a + n*x');
xx = vectorise(log(x)); yy = vectorise(log(y));
ff = fit(xx(yy>-Inf),yy(yy>-Inf),m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(ax(3),x,y,'k+')
plot(ax(3),x,exp(ff(log(x))),'r')
ylabel(ax(3),'Probability')
set(ax(3),'YScale','log','XScale','log','XTick',[1e1 1e2 1e3 1e4],'XLim',[10 1.1e4])
xlabel(ax(3),'Whiff duration (ms)')

% show the blank durations 
ax(4) = subplot(1,4,4); hold on
whiff_durations = [];
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) < odour_thresh);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1)  =[];
y = y/sum(y);
a = 1; m = fittype('a*(x).^n');
ff = fit(x(a:end)',y(a:end)',m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(ax(4),x,y,'k+')
plot(ax(4),x,ff(x),'r')
set(ax(4),'YScale','log','XScale','log','XTick',[1e2 1e3 1e4],'XLim',[30 3e4])
xlabel('Blank duration (ms)')
ylabel(ax(4),'Probability')

prettyFig('fs',.5,'FixLogX',false,'font_units','centimeters')



% fix some plots
ax(1).Position(1)=.11;
ax(2).Position(2) = ax(1).Position(2);
ax(2).Position(4) = ax(1).Position(4);

for i = 1:4
	ax(i).Position(2) = .175;
	ax(i).Position(4) = .725;
	deintersectAxes(ax(i))
end

if being_published
	snapnow
	delete(gcf)
end

return


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

% now also show the gain-per-whiff vs. the mean stimulus in the preceding 200ms for the white-noise-projection
shat = computeSmoothedStimulus(mean(ab3.PID,2),history_length);
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
ylabel('ab3A firing gain (Hz/V)')
xlabel(['Mean stimulus in preceding ' oval(history_length) 'ms'])
title(['Gain estimation using' char(10) 'white-noise filter'])

% now show ab2A filters
subplot(2,4,5), hold on
plot(filtertime,ab2.K,'b')
ylabel('ab2A filter')
xlabel('Lag (s)')

subplot(2,4,6), hold on
shat = computeSmoothedStimulus(mean(ab2.PID,2),history_length);
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
shat = computeSmoothedStimulus(mean(ab2.PID,2),history_length);
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