% makeFig2.m
% makes figure 2 for the paper
% 
% created by Srinivas Gorur-Shandilya at 6:52 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

% this script uses dataManager to ensure data integrity
dm = dataManager;
history_length = 300; % ms


%% Figure 1: Gain changes with a naturalistic stimulus
% In this figure, we show that the gain of ab3A and ab2A ORNs changes dramatically in response to a naturalistic stimulus, and that this gain change can be correlated to the mean or the variance of the stimulus in the last 200ms. 

clearvars -except being_published dm history_length

scatter_size = 12;

% make figure placeholders 
fig_handle = figure('PaperUnits','centimeters','PaperSize',[14 14],'Position',[100 100 800 801],'toolbar','figure'); hold on
clf(fig_handle);

clear ax
ax.stim = subplot(3,3,1:3); hold on
ax.resp = subplot(3,3,4:6); hold on
ax.L = subplot(3,3,7); hold on
ax.N = subplot(3,3,8); hold on
ax.gain = subplot(3,3,9); hold on


% first, grab the ab3 and ab2 data

%         ######   ######## ########       ###    ########   #######     ###    
%        ##    ##  ##          ##         ## ##   ##     ## ##     ##   ## ##   
%        ##        ##          ##        ##   ##  ##     ##        ##  ##   ##  
%        ##   #### ######      ##       ##     ## ########   #######  ##     ## 
%        ##    ##  ##          ##       ######### ##     ##        ## ######### 
%        ##    ##  ##          ##       ##     ## ##     ## ##     ## ##     ## 
%         ######   ########    ##       ##     ## ########   #######  ##     ## 

clear ab3 ab2
load(getPath(dataManager,'5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
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
shat_mean = computeSmoothedStimulus(mean(ab3.PID,2),history_length);


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
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];
std_stim(rm_this) = [];
whiff_ends(rm_this) = [];
whiff_starts(rm_this) = [];

% save this for a later fit
ab3.mean_stim = mean_stim;
ab3.std_stim = std_stim;
ab3.gain = gain;
ab3.gain_err = gain_err;
ab3.ex1 = 36.3280;
ab3.ex2 = 37.2890;

clear ff gain gain_err Kdiff mean_stim rm_this shat_mean shat_std std_stim temp whiff_starts whiff_ends R



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
shat_mean = computeSmoothedStimulus(mean(ab2.PID,2),history_length);


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
plot(ax.stim,tA(1:10:end),mean(ab3.PID(1:10:end,:),2),'Color',[0.2 .2 .2]);
set(ax.stim,'XLim',[0 70],'YLim',[0 7],'XTick',[])

% show the filter
plot(ax.L,filtertime,1e3*ab3.K,'r');

% plot the response and the prediction
clear l
[axyy,plot1,plot2] = plotyy(ax.resp,tA,mean(ab3.fA,2),tA,ab3.fp);
set(axyy(1),'XLim',[0 70],'YLim',[0 150],'YColor','k')
set(axyy(2),'XLim',[0 70],'YLim',[0 6],'YColor','r','YTick',[0 2 4 6])
set(plot1,'Color','k')
set(plot2,'Color','r')

set(ax.resp,'box','off')

shat = computeSmoothedStimulus(mean(ab3.PID,2),history_length);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
axes(ax.N)
cc = parula(100);
c = cc(shat,:);
ab3.c = c;
scatter(ab3.fp,mean(ab3.fA,2),scatter_size,c,'filled')

shat = computeSmoothedStimulus(mean(ab3.PID,2),history_length);
ch = colorbar('east');
ch.Position = [0.582    0.14    0.01    0.06];
caxis([min(shat) max(shat)]);


% plot whiff gain vs. mean stimulus preceding whiff
l(1) = plot(ax.gain,ab3.mean_stim,ab3.gain,'k+');
l(2) = plot(ax.gain,ab2.mean_stim,ab2.gain,'ko');
legend1 = legend(l,{'ab3A','ab2A'},'Location','southwest');
set(ax.gain,'XScale','log','YScale','log')

% fit a inverse relationship to all the data
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./[ab2.gain_err; ab3.gain_err];
ff = fit([ab2.mean_stim; ab3.mean_stim],[ab2.gain; ab3.gain],'power1',options);
plot(ax.gain,sort([ab2.mean_stim; ab3.mean_stim]),ff(sort([ab2.mean_stim; ab3.mean_stim])),'r');


% now show that the variance and the mean are correlated 
all_block_sizes = factor2(length(ab3.PID));
all_block_sizes = all_block_sizes(6:end-1);
all_block_sizes = all_block_sizes(1:41);
clear l r2
all_block_sizes = unique([all_block_sizes history_length]);
r2 = NaN*all_block_sizes;

% save this data for later
nat_stim_stats.all_block_sizes = all_block_sizes;

for i = 1:length(all_block_sizes)
	temp = ab3.PID(:);
	temp = reshape(temp,all_block_sizes(i),length(temp)/all_block_sizes(i));
	if all_block_sizes(i) == history_length
		nat_stim_stats.history_length = history_length;
		nat_stim_stats.mu = mean(temp);
		nat_stim_stats.sigma = std(temp);
		% plot(axes_handles(9),mean(temp),std(temp),'Marker','.','Color',[.5 .5 .5],'LineStyle','none')
	end
	r2(i) = rsquare(mean(temp),std(temp));
end
nat_stim_stats.r2 = r2;

% label all the axes, set the scales, etc.

ylabel(axyy(1),'ab3A firing rate (Hz)')
ylabel(axyy(2),'Projected stimulus (V)')
ylabel(ax.stim,'Stimulus (V)')
xlabel(ax.resp,'Time (s)')
xlabel(ax.N,'Projected stimulus (V)')
ylabel(ax.N,'ab3A firing rate (Hz)')
set(ax.N,'XLim',[-.1 5.1],'YColor','k','XColor','r','box','off')
xlabel(ax.gain,['\mu_{Stimulus} in preceding ' oval(history_length) 'ms (V)'])
ylabel(ax.gain,'ORN gain (Hz/V)')
set(ax.gain,'YLim',[10 1e4],'YScale','log','XScale','log','XLim',[.001 10],'XTick',[.001 .01 .1 1 10])
xlabel(ax.L,'Filter lag (s)')
ylabel(ax.L,'Filter (a.u.)')
set(ax.N,'box','off','XTick',[0 2 4 6],'box','off')
set(ax.L,'XLim',[-.2 1])


% fix position of top two rows
ax.resp.Position = [0.1300    0.4096    0.7750    0.19];
ax.stim.Position = [0.1300    0.64    0.7750    0.19];

prettyFig('plw',1.5,'lw',1.5,'fs',.5,'FixLogX',true,'font_units','centimeters','x_minor_ticks',false,'y_minor_ticks',false)

% add some labels
labelFigure;

% compress some plots to make small EPS files
shrinkDataInPlot(ax.stim,1);
shrinkDataInPlot(ax.resp,2);
shrinkDataInPlot(axyy(2),2);
shrinkDataInPlot(ax.L,2)
shrinkDataInPlot(ax.N,1)

ax.gain.XLim(2) = 11;
deintersectAxes(ax.gain)
deintersectAxes(ax.L)

if being_published
	snapnow
	delete(gcf)
end

% save some stuff for later
save('nat_stim_stats','nat_stim_stats')
save('nat_stim_data','ab3','ab2')


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
figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on

fs = 25; % font size

% show the PDF of the whiff intensities 
odour_thresh = 0.024;

% whiff intensity distribution  
ax(2) = subplot(2,3,1); hold on
whiff_intensities = []; 
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) > odour_thresh);
	whiff_intensities =  [whiff_intensities; findMeanInWindows(ons,offs,ab3.PID(:,i));];
end

[y,x] = histcounts(whiff_intensities,50); x(1) = [];
y = y/sum(y);
plot(ax(2),x,y,'k+')
set(ax(2),'YScale','log','XScale','log','XTick',[1e-2 1e-1 1 10],'XLim',[1e-2 11])
ylabel(ax(2),'Probability')
xlabel(ax(2),'Whiff intensity (V)')

m = fittype('log((a./x).*exp(-x./b))');
ff = fit(x(y>0)',log(y(y>0))',m,'Upper',[10 100],'Lower',[0 1],'StartPoint',[7e-3 4]);
plot(ax(2),sort(x),exp(ff(sort(x))),'r')

th(2) = text(.1, .2,'$\sim\frac{1}{c}\exp\left(-\frac{c}{C}\right)$','interpreter','latex','Color','r','FontSize',fs);

% show the whiff durations 
ax(3) = subplot(2,3,2); hold on
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

axes(ax(3))
th(3) = text(1e3, .1,'$\sim t_{w}^{-\frac{3}{2}}$','interpreter','latex','Color','r','FontSize',fs);

% show the blank durations 
ax(4) = subplot(2,3,3); hold on
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

axes(ax(4))
th(4) = text(1e3, .1,'$\sim t_{b}^{-\frac{3}{2}}$','interpreter','latex','Color','r','FontSize',fs);


% ab3A filters and projections
ax(5) = subplot(2,3,4); hold on
plot(ax(5),filtertime,ab3.K/max(max(ab3.K)),'r');

% load MSG filter for ab3
load(getPath(dataManager,'b934ad74a78cb20df400b67e107037e3'));
ab3.K_white = mean(reshape([MSG_data(1,:).K],1001,11),2);

% also show the Ky filter from the DA model
clear p
p.   s0 = 0;
p.  n_z = 2;
p.tau_z = 147.3750;
p.  n_y = 2;
p.tau_y = 27.2500;
p.    C = 0.5000;
p.    A = 170.4375;
p.    B = 2.7656;
[~,~,~,Ky,Kz] = DAModelv2(ones(1e4,1),p);
plot(ax(5),1e-3*(1:length(Ky)),Ky/max(Ky),'b')


t = 1e-3*(1:length(ab3.K_white))- .2;
plot(ax(5),t,ab3.K_white/max(ab3.K_white),'k')
legend({'Natural stimulus','K_y from DA model','Gaussian stimulus'})
ylabel(ax(5),'ab3A filter (norm)')
xlabel(ax(5),'Lag (s)')
set(ax(5),'YLim',[-.2 1.5])

ab3.fp2 = convolve(1e-3*(1:length(mean(ab3.PID,2))),mean(ab3.PID,2),ab3.K_white,t);

ax(6) = subplot(2,3,5); hold on
scatter(ax(6),ab3.fp2,mean(ab3.fA,2),scatter_size,ab3.c,'filled')
xlabel(ax(6),['Stimulus projected' char(10) ' using Gaussian filter (V)'])
ylabel(ax(6),'ab3A Response (Hz)')
set(ax(6),'XColor','k','XLim',[-1 14])

t = ['r^2=' oval(rsquare(ab3.fp2,nanmean(ab3.fA,2)))]; 
axes(ax(6))
th(5) = text(7,40,t);

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

ax(7) = subplot(2,3,6); hold on
plot(ax(7),mean_stim,gain,'ko');
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
ff = fit(mean_stim,gain,'power1',options);
ff.a = 8;
plot(ax(7),mean_stim,ff(mean_stim),'r')
set(ax(7),'XScale','log','YScale','log','XTick',[1e-2 1e-1 1e0 1e1])
ylabel(ax(7),'ab3A firing gain (Hz/V)')
xlabel(ax(7),['Mean stimulus in preceding ' oval(history_length) 'ms'])

prettyFig('fs',.6,'FixLogX',true,'font_units','centimeters')

labelFigure('delete_all',true)
labelFigure('font_size',fs)

for i = 2:7
	deintersectAxes(ax(i))
end

for i = 2:5
	th(i).FontSize = fs;
end

if being_published
	snapnow
	delete(gcf)
end

return

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
ch.Position = [.58 .15 .02 .1];
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