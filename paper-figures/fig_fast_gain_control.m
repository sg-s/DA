% Fast Gain Control
% 
% created by Srinivas Gorur-Shandilya at 3:28 , 24 September 2015. Contact me at http://srinivas.gs/contact/
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


%% Fast Gain Control in ORNs
% 

fig_handle=figure('outerposition',[0 0 900 1000],'PaperUnits','points','PaperSize',[900 1000]); hold on
clf(fig_handle);
axes_handles(1) = subplot(9,3,1:6); 
axes_handles(2) = subplot(9,3,7:12); 
axes_handles(3) = subplot(9,3,13:18);
axes_handles(4) = subplot(9,3,[19 22 25]);
axes_handles(5) = subplot(9,3,[20 23 26]);
axes_handles(6) = subplot(9,3,[21 24 27]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end

movePlot(axes_handles(1),'up',.03)
movePlot(axes_handles(2),'up',.04)
movePlot(axes_handles(3),'up',.02)

movePlot(axes_handles(4),'down',.02)
movePlot(axes_handles(5),'down',.02)
movePlot(axes_handles(6),'down',.02)

load('/local-data/DA-paper/large-variance-flicker/ab3/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/ab3/2015_01_28_CS_ab3_3_EA.mat')
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
plot(axes_handles(1),tA,mean(PID,2),'k')
set(axes_handles(1),'XLim',[0 60],'XTick',[])
ylabel(axes_handles(1),'Stimulus (V)')

% make linear predictions for the whole experiment
K = NaN(1e3,1);
[this_K, filtertime_full] = fitFilter2Data(mean(PID(10e3:end,:),2),mean(fA(10e3:end,:),2),'reg',1,'filter_length',1200,'offset',200);
K = this_K(100:end-101);
filtertime = filtertime_full(100:end-101);

fp = NaN*fA;
for i = 1:width(fA)
	fp(:,i) = convolve(tA,PID(:,i),K,filtertime);
end

clear l
R = mean(fA,2);
[ax,plot1,plot2] = plotyy(axes_handles(2),tA,R,tA,mean(fp,2));
set(ax(1),'XLim',[0 60],'YLim',[min(mean(fA,2)) 1.3*max(mean(fA,2))])
set(ax(2),'XLim',[0 60],'YLim',[min(mean(fp,2)) max(mean(fp,2))])
set(ax(2),'YTick',-.1:.1:.5)
set(plot1,'Color','k')
set(plot2,'Color','r')
set(ax(1),'YColor',[0 0 0])
ylabel(ax(1),'ORN Response (Hz)')
ylabel(ax(2),'Projected Stimulus')
set(axes_handles(2),'box','off')



hl_min = .1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

resp = mean(fA(10e3:55e3,[3:10 13:20]),2);
pred = mean(fp(10e3:55e3,[3:10 13:20]),2);
time = 1e-3*(1:length(resp));
stim = mean(PID(10e3:55e3,[3:10 13:20]),2);


% plot gain vs preceding stimulus
global inst_gain
[mean_stim,inst_gain,e] = makeFig6G(mean(PID,2),mean(fA,2),mean(fp,2),500);
rm_this = (isnan(mean_stim) | isnan(inst_gain)) | mean_stim < .2 | e < .8;
mean_stim(rm_this) = NaN;
inst_gain(rm_this) = NaN;


plot(axes_handles(3),tA(1:10:end),inst_gain(1:10:end),'k.')
ylabel(axes_handles(3),'Inst. Gain (Hz/V)')
xlabel(axes_handles(3),'Time (s)')
set(axes_handles(3),'XLim',[0 60])

% show the backed-out gain filter
% these parameters come from fitting a gain filter from the PID to the inst. gain
clear p
p.   A = 0.0386;
p.tau1 = 299.8428;
p.tau2 = 300.3500;
p.   n = 0.7266;

K = filter_gamma2(1:2e3,p);
plot(axes_handles(4),1e-3*(1:length(K)),K/sum(K),'b')
xlabel(axes_handles(4),'Filter Lag (s)')
ylabel(axes_handles(4),'Gain filter')

shat = filter(K,sum(K),mean(PID,2));
shat(isnan(inst_gain)) = NaN;

axis(axes_handles(5)); hold on
[~,data] = plotPieceWiseLinear(shat,inst_gain,'nbins',50,'Color','b','use_std',true,'make_plot',false);
l = errorShade(axes_handles(5),data.x,data.y,data.ye,'Color',[0 0 1]);
temp1 = shat(~isnan(inst_gain));
temp2 = inst_gain(~isnan(inst_gain));
legend(l(1),['\rho =' oval(spear(temp1(1:10:end),temp2(1:10:end)),3)])
xlabel(axes_handles(5),'Stimulus projected by gain filter (V)')
ylabel(axes_handles(5),'Inst. Gain (Hz/V)')
set(axes_handles(5),'YScale','log','XLim',[min(shat) max(shat)])
set(axes_handles(5),'YLim',[40 1000])

% now vary the gain filter and show it is the best possible
all_A = linspace(0,1,30);
all_rho = NaN*all_A;
diff_degree = NaN*all_A;
all_tau = ceil(logspace(1,3,30));
all_rho_tau = NaN*all_tau;


for i = 1:length(all_A)
	q = p;
	q.A = all_A(i);
	q.tau1 = p.tau1*(1-(all_A(i)/2));
	q.tau2 = p.tau1*(1+(all_A(i)/2));
	[~,all_rho(i)] = findBestGainFilter(mean(PID,2),q);
	K = filter_gamma2(1:2e3,q);
	diff_degree(i) = sum(K)/sum(abs(K));

	q = p;
	q.A = 0;
	q.tau1 = all_tau(i);
	q.tau2 = all_tau(i);
	[~,all_rho_tau(i)] = findBestGainFilter(mean(PID,2),q);
end


ax = plotyy(axes_handles(6),all_rho,diff_degree,all_rho_tau,all_tau);

ylabel(ax(1),'\int {K} / \int {|K|}','interpreter','tex')
xlabel(axes_handles(6),'\rho') 
ylabel(ax(2),'\tau_{gain} (ms)')


prettyFig('plw=1.5;','lw=1.5;','fs=14;','FixLogX=1;')

if being_published
	snapnow
	delete(gcf)
end


%% Supplementary Figure
% 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,2,1), hold on


% [~,~,~,~,~,history_lengths]=gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);
% set(axes_handles(5),'XLim',[.09 11]) % to show .1 and 10 on the log scale
% set(axes_handles(4),'XLim',[min(pred) max(pred)])
% xlabel(axes_handles(4),'Projected Stimulus (V)')
% ylabel(axes_handles(4),'ORN Response (Hz)')

% % show the p-value
% axes(axes_handles(4))
% text(-.1,60,'p < 0.01')
% title(axes_handles(4),[])

axes(axes_handles(3))
plotPieceWiseLinear(mean_stim,inst_gain,'nbins',50,'use_std',true);
xlabel(axes_handles(3),'Stimulus in preceding 500ms (V)')
ylabel(axes_handles(3),'Inst. Gain (Hz/V)')
set(axes_handles(3),'YScale','log','XScale','log','YTick',[10 100 1000],'YLim',[10 1000])

% fix some labels
set(axes_handles(5),'YScale','log','YTick',[0.5 1 2],'YLim',[0.4 2.5],'XLim',[0.09 10.1],'YMinorTick','on')

% Indicate regions of high and low stimulus on the stimulus
hl = round(history_lengths(9)*1e3);
shat = computeSmoothedStimulus(mean(PID,2),hl);

n = floor(sum(~isnan(mean(fA,2)))*.33);
shat(1:hl) = Inf; % the initial segment where we can't estimate shat is excluded
shat(isnan(shat)) = Inf;
shat(isnan(mean(fA,2))) = Inf;
[~, t_low] = sort(shat,'ascend');
t_low = t_low(1:n); % this is an index
t_low = tA(t_low); % t_low is now a time. 
 
shat = computeSmoothedStimulus(mean(PID,2),hl);
shat(1:hl) = -Inf;
shat(isinf(shat)) = -Inf;
shat(isnan(mean(fA,2))) = -Inf;
[~, t_high] = sort(shat,'descend');
t_high = t_high(1:n);
t_high  = tA(t_high);

plot(axes_handles(1),t_low,1+0*t_low,'.','Color',c(1,:))
plot(axes_handles(1),t_high,1+0*t_low,'.','Color',c(7,:))




%% Test Figure: How does instantaneous gain depend on various history lengths? 
% In this section, we look at a generalized version of the plot of instate nous gain vs. mean history in the last 500ms, where we now vary the length of the history length, and see how the slope of the gain plot changes, as well as how well it is fit by a straight line. First, we plot the inst. gain vs. stimulus projected by a box filter, for various lengths of the box filter.


clear inst_gain
global inst_gain
[~,inst_gain] = makeFig6G(mean(PID,2),mean(fA,2),mean(fp,2),500);


mean_pid = mean(PID,2);
history_lengths = logspace(-2,1,40);

gain_K1_slope = NaN*history_lengths;
gain_K1_rho = NaN*history_lengths;

gain_K2_slope = NaN*history_lengths;
gain_K2_rho = NaN*history_lengths;

gain_K3_slope = NaN*history_lengths;
gain_K3_rho = NaN*history_lengths;

c = parula(length(history_lengths)+1);
for i = 1:length(history_lengths)
	% first do the simple box filter
	temp = floor(history_lengths(i)*1e3);
	proj_stim = filter(ones(temp,1),temp,mean_pid);

	rm_this = (isnan(proj_stim) | isnan(inst_gain) | inst_gain < 0);
	proj_stim(rm_this) = [];
	y = inst_gain;
	y(rm_this) = [];

	gain_K1_rho(i) = spear(proj_stim(1:10:end),y(1:10:end));

	ff = fit(proj_stim(:),y(:),'power1','StartPoint',[0 0]);
	gain_K1_slope(i) = ff.b;


	% now do a squared differentiating filter
	temp = floor(history_lengths(i)*1e3/2);
	proj_stim = filter([ones(temp,1); -ones(temp,1)],2*temp,mean_pid);
	proj_stim = proj_stim.^2;

	rm_this = (isnan(proj_stim) | isnan(inst_gain) | inst_gain < 0);
	proj_stim(rm_this) = [];
	y = inst_gain;
	y(rm_this) = [];

	gain_K2_rho(i) = spear(proj_stim(1:10:end),y(1:10:end));

	ff = fit(proj_stim(:),y(:),'power1','StartPoint',[0 0]);
	gain_K2_slope(i) = ff.b;

	% now do a absolute value of differentiating filter
	temp = floor(history_lengths(i)*1e3/2);
	proj_stim = filter([ones(temp,1); -ones(temp,1)],2*temp,mean_pid);
	proj_stim = abs(proj_stim);

	rm_this = (isnan(proj_stim) | isnan(inst_gain) | inst_gain < 0);
	proj_stim(rm_this) = [];
	y = inst_gain;
	y(rm_this) = [];

	gain_K3_rho(i) = spear(proj_stim(1:10:end),y(1:10:end));

	ff = fit(proj_stim(:),y(:),'power1','StartPoint',[0 0]);
	gain_K3_slope(i) = ff.b;

end


%%
% We now plot the slopes of these plots as a function of history length, and also plot the Spearman rank correlation, showing where these slopes are meaningful (and where there is little correlation between inst. gain and the projected stimulus). 

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
h1 = subplot(1,3,1); hold on
ax = plotyy(history_lengths,gain_K1_slope,history_lengths,gain_K1_rho);
ylabel(ax(1),'Gain Exponent')
ylabel(ax(2),'\rho')
xlabel('History Length (s)')
title('Integrating Filter')
set(ax,'XScale','log')
set(ax(2),'YLim',[-1 0],'YTick',-1:.2:0)
set(ax(1),'YLim',[-3 1],'YTick',-3:.5:1)

subplot(1,3,2), hold on
ax = plotyy(history_lengths,gain_K2_slope,history_lengths,gain_K2_rho);
ylabel(ax(1),'Gain Exponent')
ylabel(ax(2),'\rho')
xlabel('History Length (s)')
title('Sq. Differentiating Filter')
set(ax,'XScale','log')
set(ax(2),'YLim',[-1 0],'YTick',-1:.2:0)
set(ax(1),'YLim',[-3 1],'YTick',-3:.5:1)

h3= subplot(1,3,3); hold on
ax = plotyy(history_lengths,gain_K3_slope,history_lengths,gain_K3_rho);
ylabel(ax(1),'Gain Exponent')
ylabel(ax(2),'\rho')
xlabel('History Length (s)')
title('Abs. Differentiating Filter')
set(ax,'XScale','log')
set(ax(2),'YLim',[-1 0],'YTick',-1:.2:0)
set(ax(1),'YLim',[-3 1],'YTick',-3:.5:1)

prettyFig('FixLogX=1;');
pause(1)
movePlot(h1,'left',.05)
movePlot(h3,'right',.025)

if being_published
	snapnow
	delete(gcf)
end

%% Reconstructing the gain filter
% In the previous section, we have used a square filter whose length we varied. In this section, we generalize the analysis by parameterizing the filter shape, and finding the best parameters that can account for the observed variation in instantaneous gain. The following figure shows the best-fit filter from this case (parameterized by a sum of two gamma functions). As we can see, a purely integrating filter is chosen, with a time-scale roughly comparable to the timescale predicted by the previous analysis. 




%%
% The results of our optimization technique show that the "best" gain filter is one that integrates the stimulus, and doesn't differentiate it at all. However, are we sure that our optimization has found the global minimum? To verify that the gain filters with differentiating components cannot do a better job at explaining observed gain changes, we manually change the degree to which the filter differentiates, and then determine what fraction of the variance this manipulations has. 

%%
% In the following figure, in the first panel, the Spearman rank-correlation between the projected stimulus and the instatnenous gain is plotted on the Y-axis. Lower (larger absolute) values mean that it can explain more of the data. On the x-axis is a measure of how differentiating the filter is. Values close to 0 mean filters that perform perfect differentiation, while values close to 1 mean filters that are integrators. 

%% 
% In the second panel, we repeat the analysis, but now we vary the timescale of the filter. 

all_A = linspace(0,1,30);
all_rho = NaN*all_A;
diff_degree = NaN*all_A;
all_tau = ceil(logspace(1,3,30));
all_rho_tau = NaN*all_tau;


for i = 1:length(all_A)
	q = p;
	q.A = all_A(i);
	q.tau1 = p.tau1*(1-(all_A(i)/2));
	q.tau2 = p.tau1*(1+(all_A(i)/2));
	[~,all_rho(i)] = findBestGainFilter(mean_pid,q);
	K = filter_gamma2(1:2e3,q);
	diff_degree(i) = sum(K)/sum(abs(K));

	q = p;
	q.A = 0;
	q.tau1 = all_tau(i);
	q.tau2 = all_tau(i);
	[~,all_rho_tau(i)] = findBestGainFilter(mean_pid,q);
end



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(diff_degree,all_rho,'k+')
xlabel('\int {K} / \int {|K|}','interpreter','tex')
ylabel('\rho') 


subplot(1,2,2), hold on
plot(all_tau*p.n,all_rho_tau,'k+')
xlabel('\tau (ms)','interpreter','tex')
ylabel('\rho') 

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Supplementary Figure: LN models cannot account for this
% In this section we show that a static output nonlinearity cannot account for observed fast gain control. 

fp = mean(fp,2);
temp = fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;

clear p
p.A = 57.1534;
p.k = 23.6690;
p.n = 2.9341;
fp_hill = hill(p,fp);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on
ph(4) = subplot(1,2,2); hold on

hl_min = .1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

resp = mean(fA(10e3:55e3,[3:10 13:20]),2);
pred = (fp_hill(10e3:55e3));
time = 1e-3*(1:length(resp));
stim = mean(PID(10e3:55e3,[3:10 13:20]),2);

[p,~,~,~,~,history_lengths]=gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);
set(ph(4),'XLim',[.09 11]) % to show .1 and 10 on the log scale
set(ph(3),'XLim',[min(pred) mean(max(fA))],'YLim',[min(pred) mean(max(fA))])
xlabel(ph(3),'LN Prediction (Hz)')
ylabel(ph(3),'ORN Response (Hz)')
title(ph(3),'')
prettyFig('FixLogX=1;')

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

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end