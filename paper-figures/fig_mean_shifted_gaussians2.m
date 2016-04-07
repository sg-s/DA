% GainChangesWithMeanStimulus.m
% 
% created by Srinivas Gorur-Shandilya at 5:49 , 23 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;



%     ##     ## ########    ###     ######  ##     ## ########  #### ##    ##  ######   
%     ###   ### ##         ## ##   ##    ## ##     ## ##     ##  ##  ###   ## ##    ##  
%     #### #### ##        ##   ##  ##       ##     ## ##     ##  ##  ####  ## ##        
%     ## ### ## ######   ##     ##  ######  ##     ## ########   ##  ## ## ## ##   #### 
%     ##     ## ##       #########       ## ##     ## ##   ##    ##  ##  #### ##    ##  
%     ##     ## ##       ##     ## ##    ## ##     ## ##    ##   ##  ##   ### ##    ##  
%     ##     ## ######## ##     ##  ######   #######  ##     ## #### ##    ##  ######   
    
%      ######      ###    #### ##    ## 
%     ##    ##    ## ##    ##  ###   ## 
%     ##         ##   ##   ##  ####  ## 
%     ##   #### ##     ##  ##  ## ## ## 
%     ##    ##  #########  ##  ##  #### 
%     ##    ##  ##     ##  ##  ##   ### 
%      ######   ##     ## #### ##    ## 


%% Figure 1: ORN gain can be estimated by measuring responses to Gaussian inputs
% Gaussian odorant inputs with short correlation times (A), elicit flickering responses in ORNs that track the odorant stimulus well (B). A linear filter K can be extracted from the odorant input and the firing rate output of the neuron (C). The slope of the residuals in a plot of the firing response vs. the linear prediction (D) is defined as the gain of the ORN. Here, we measure the ORN gain in the linear regime: the linear filter accounts for 96% of the variance in the ORN response (red line), and adding an output nonlinearity (dotted black line), only accounts for an additional 1%. The odorant used is ethyl acetate, stimulating the ab3A neuron. Shading in all plots shows the standard error of the mean. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
clear axes_handles
axes_handles(1) = subplot(3,10,1:10);
axes_handles(2) = subplot(3,10,11:20);

axes_handles(3) = subplot(3,3,7); 
axes_handles(4) = subplot(3,3,8);
axes_handles(5) = subplot(3,3,9);
for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end  
% load data 

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
LFP(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% extract filters and find gain
a = 35e3; z = 55e3;
[K,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% some core variables
dt = 1e-3;
c = parula(max(paradigm)+1); % colour scheme

% plot some stimulus
example_dose = 1;
plot_this = nanmean(PID(:,paradigm==example_dose),2);
time = dt*(1:length(plot_this));
axes(axes_handles(1))
plot(time,plot_this,'Color',c(example_dose,:));
ylabel(axes_handles(1),'Stimulus (V)')
set(axes_handles(1),'XLim',[35 55])


% plot the filter for the lowest dose 
filtertime = (-100:599)*1e-3;
axes(axes_handles(3))
shadedErrorBar(filtertime,mean(K(:,paradigm==example_dose),2),sem(K'),{'Color',c(example_dose,:)})
set(axes_handles(3),'XLim',[-.1 .6])
xlabel(axes_handles(3),'Lag (s)')
ylabel(axes_handles(3),'Filter K (norm)')

% plot linear prediction vs. data for the example dose. 
y = mean(fA(a:z,paradigm==example_dose),2);
x = mean(fA_pred(a:z,paradigm==example_dose),2);
time = 35+dt*(1:length(x));

[ax,plot1,plot2] = plotyy(axes_handles(2),time,y,time,x);
set(ax(1),'XLim',[35 55],'YLim',[min(y(5e3:end)) max(y(5e3:end))],'YColor',c(example_dose,:))
set(ax(2),'XLim',[35 55],'YLim',[min(x(5e3:end)) max(x(5e3:end))],'YColor','r')
set(plot1,'Color',c(example_dose,:))
set(plot2,'Color','r')
ylabel(ax(1),'ORN Response (Hz)')
ylabel(ax(2),'Projected Stimulus (V)')
set(axes_handles(2),'box','off')
set(ax(2),'XMinorTick','on','YMinorTick','on','YTick',[0:0.1:0.5],'YTickLabel',{'0','.1','.2','.3','.4','.5'})
set(ax(1),'XMinorTick','on','YMinorTick','on','YTick',[0:10:50],'YTickLabel',{'0','10','20','30','40','50'})
xlabel('Time (s)')

plot(axes_handles(4),x(1:25:end),y(1:25:end),'.','Color',c(example_dose,:));
ff= fit(x(~isnan(x)),y(~isnan(x)),'poly1');
l = plot(axes_handles(4),sort(x),ff(sort(x)),'k');
L = {};
L{1} = strcat('Gain=',oval(ff.p1),'Hz/V, r^2=',oval(rsquare(y,x)));


xlabel(axes_handles(4),'Projected Stimulus (V)')
ylabel(axes_handles(4),'ORN Response (Hz)')
set(axes_handles(4),'YColor',c(example_dose,:),'XColor','r')

lh = legend(l,L,'Location','northwest');

movePlot(axes_handles(2),'up',.02)

plot(axes_handles(5),nanmean(PID(a:z,:)),nanstd(PID(a:z,:)),'k+')
plot(axes_handles(5),[0 2],[0 2],'k--')
xlabel(axes_handles(5),'\mu_{stimulus} (V)')
ylabel(axes_handles(5),'\sigma_{stimulus (V)}')
set(axes_handles(5),'XLim',[0 2],'YLim',[0 2])

prettyFig('plw',1.3,'lw',1.5,'fs',14)
set(lh,'Position',[0.42 0.32 0.2112 0.0275],'box','off')
labelFigure
xlabel(axes_handles(2),'Time (s)')

if being_published
	snapnow
	delete(gcf)
end




%    ##      ## ######## ########  ######## ########      ######      ###    #### ##    ## 
%    ##  ##  ## ##       ##     ## ##       ##     ##    ##    ##    ## ##    ##  ###   ## 
%    ##  ##  ## ##       ##     ## ##       ##     ##    ##         ##   ##   ##  ####  ## 
%    ##  ##  ## ######   ########  ######   ########     ##   #### ##     ##  ##  ## ## ## 
%    ##  ##  ## ##       ##     ## ##       ##   ##      ##    ##  #########  ##  ##  #### 
%    ##  ##  ## ##       ##     ## ##       ##    ##     ##    ##  ##     ##  ##  ##   ### 
%     ###  ###  ######## ########  ######## ##     ##     ######   ##     ## #### ##    ## 


%% Figure 2: ORN gain decreases with increasing stimulus intensity, similar to the Weber-Fechner Law
% Odorant stimuli drawn from distributions with similar variances but increasing means (A) elicit ORN responses with decreasing variances (B). After extracting linear filters for all stimulus paradigms, a plot of the ORN response vs. the linear prediction (C) shows a systematic decrease in slope. Plotting lines to each of these clouds of points determines the neuron gain in each case. Neuron gain decreases with the mean stimulus (D). This stimulus-dependent decrease in gain is well described by a power law with an exponent close to -1 (the Weber-Fechner Prediction). For comparison, a power law with the exponent fixed at -1 is also shown (dashed line). 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
clear ax
ax(1) = subplot(3,10,1:8);
ax(2) = subplot(3,10,9:10);
ax(3) = subplot(3,10,11:18);
ax(4) = subplot(3,10,19:20);

ax(5) = subplot(3,3,7); 
ax(6) = subplot(3,3,8);
ax(7) = subplot(3,3,9);
for i = 1:length(ax)
	hold(ax(i),'on');
end  

% plot all the stimuli
for i = 1:max(paradigm)
	plot(ax(1),time,nanmean(PID(a:z,paradigm==i),2),'Color',c(i,:))
end

% plot the stimulus distributions 
a = floor(35/dt);
z = floor(55/dt);

for i = 1:max(paradigm)
	plot_hist = PID(a:z,paradigm == i);
	[hy,hx]  = hist(plot_hist(:),50);
	hy = hy/sum(hy);
	plot(ax(2),hy,hx,'Color',c(i,:));
end

% plot the responses
for i = 1:max(paradigm)
	temp = fA(a:z,paradigm==i);
	plot(ax(3),time,nanmean(temp,2),'Color',c(i,:));
	[hx,hy] = hist(mean(temp,2),30);
	hx = hx/sum(hx);
	plot(ax(4),hx,hy,'Color',c(i,:));
end


% show gain changes for all paradigms -- average over neurons 
ss = 100;
all_x = 0:0.1:2;
axes(ax(5)), hold(ax(5),'on')
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	[~,orn_io_data(i)] = plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end

mean_stim = nanmean(PID(a:z,:));

% show gain changes -- gain vs. mean stimulus
for i = 1:max(paradigm)
	plot(ax(6),mean_stim(paradigm==i),fA_gain(paradigm==i),'+','Color',c(i,:));
end


% fit a power law with exponent -1
mean_stim = mean_stim(:);
fA_gain = fA_gain(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(mean_stim(~isnan(fA_gain)),fA_gain(~isnan(fA_gain)),'power1',options);
l = plot(ax(6),sort(mean_stim),cf(sort(mean_stim)),'r');
r2 = rsquare(nonnans(fA_gain),cf(nonnans(mean_stim)));
set(ax(6),'XScale','log','YScale','log','YLim',[10 300],'XLim',[.1 2.5])

% rescale by Weber law
ss = 50;
allx = [];
ally = [];
axes(ax(7)), hold(ax(7),'on')
for i = 1:8 % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);

	allx = [allx mean(y)+cf(mean(s))*(x(1:ss:end))];
	ally = [ally y(1:ss:end)];

	x = cf(nanmean(s))*(x);
	x = x - nanmean(x);
	x = x + nanmean(y);
	plotPieceWiseLinear(x(1e3:end),y(1e3:end),'nbins',40,'Color',c(i,:));
end

plot(ax(7),[0 80],[0 80],'r')


% cosmetics
set(ax(1),'XLim',[35 55],'YLim',[0 2])
set(ax(2),'YLim',[0 2],'YTick',[])

set(ax(3),'XLim',[35 55],'YLim',[0 70])
set(ax(4),'YLim',[0 70],'YTick',[])

set(ax(7),'XLim',[0 70],'YLim',[0 70])

ylabel(ax(1),'Stimulus (V)')
ylabel(ax(3),'ORN Response (Hz)')
xlabel(ax(3),'Time (s)')
xlabel(ax(4),'Probability')
xlabel(ax(5),'Projected Stimulus (V)')
ylabel(ax(5),'ORN Response (Hz)')
xlabel(ax(6),'Mean Stimulus (V)')
ylabel(ax(6),'ORN Gain (Hz/V)')

xlabel(ax(7),[' Projected Stimulus ' char(10) 'Rescaled by Weber Law'])
ylabel(ax(7),'ORN Response (Hz)')

prettyFig('plw',1.3,'lw',1.5,'fs',14)
labelFigure

if being_published
	snapnow
	delete(gcf)
end


%% Sanity Check
% In this section, we check if my filter extractions works (comparing it to Damon's FFT-based methods) if we still see the overall phenomenology if we only use a single filter to project all the data. 

% back out Damon's filters here
K_Damon = NaN*K;
for i = 1:width(PID)
	if (paradigm(i)==1)
		K_Damon(:,i) = backOutFilter(PID(a:z,i),fA(a:z,i),'offset',100,'filter_length',700);
	end
end

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(2,3,1), hold on
for i = 1:max(paradigm)
	temp = nanmean(K(:,paradigm==i),2);
	plot(filtertime,temp,'Color',c(i,:))
end
xlabel('Filter Lag (s)')
ylabel('Filter amplitude')
title('Filters extracted at different mean stimuli')


subplot(2,3,2), hold on
temp = nanmean(K(:,paradigm==1),2);
single_K = temp;
temp =temp/max(temp);
plot(filtertime,temp,'k')
temp = nanmean(K_Damon(:,paradigm==1),2); temp =temp/max(temp);
plot(filtertime,temp,'r')
legend({'Srinivas Filter','Damon Filter'})
xlabel('Filter Lag (s)')
ylabel('Filter (a.u.)')

% use this to make projections everywhere
fp = NaN*fA;
for i = 1:width(PID)
	fp(:,i) = convolve(1e-3*(1:length(PID)),PID(:,i),single_K,filtertime);
end

subplot(2,3,3), hold on
ss = 100;
all_x = 0:0.1:2;
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fp(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('ORN Response (Hz)')
title('Using a single filter')

% find gain in each trial
single_K_gain = NaN(width(PID),1);
single_K_r2 = NaN(width(PID),1);
r2 = NaN(width(PID),1);
for i = 1:width(PID)
	try
		temp = fit(fp(a:z,i),fA(a:z,i),'poly1');
		single_K_gain(i) = temp.p1;
		single_K_r2(i) = rsquare(fp(a:z,i),fA(a:z,i));
		r2(i) = rsquare(fA_pred(a:z,i),fA(a:z,i));
	catch
	end
end

subplot(2,3,4), hold on
plot(nanmean(PID(a:z,:)),single_K_gain,'k+')
set(gca,'XScale','log','YScale','log')

% fit a power law with exponent -1
mean_stim = nanmean(PID(a:z,:));
single_K_gain = single_K_gain(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(nonnans(mean_stim),nonnans(single_K_gain),'power1',options);
plot(sort(mean_stim),cf(sort(mean_stim)),'r');


% also fit a power law with no constraint
mean_stim = nanmean(PID(a:z,:));
single_K_gain = single_K_gain(:);
cf = fit(nonnans(mean_stim),nonnans(single_K_gain),'power1');
l = plot(sort(mean_stim),cf(sort(mean_stim)),'k--');
legend(l,['\alpha=' oval(cf.b)])

xlabel('Mean Stim (V)')
ylabel('Gain calc. using single filter (Hz/V)')

subplot(2,3,5), hold on
plot([0 1],[0 1],'k--')
plot(r2,single_K_r2,'k+')
xlabel('r^2, best filters')
ylabel('r^2 using single filter ')

% also estimate gain using variances of stimulus and response
frac_var = NaN(width(PID),1);
for i = 1:width(PID)
	try
		frac_var(i) = std(fA(a:z,i))/std(PID(a:z,i));
	catch
	end
end
subplot(2,3,6), hold on
for i = 1:width(PID)
	plot(mean_stim(i),frac_var(i),'+','Color',c(paradigm(i),:))
end
% fit a power law with exponent -1
mean_stim = nanmean(PID(a:z,:));
frac_var = frac_var(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(nonnans(mean_stim),nonnans(frac_var),'power1',options);
plot(sort(mean_stim),cf(sort(mean_stim)),'r');
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('\sigma_{Firing Rate}/\sigma_{Stimulus} (Hz/V)')

prettyFig('plw',1.3,'lw',1.5,'fs',14,'FixLogX',true)
labelFigure
if being_published
	snapnow
	delete(gcf)
end

%% Timescale of gain changes in the firing rate
% In this section, we attempt to determine the timescale of the gain changes by analyzing gain in small bins along the course of the experiment. The bin size is 50ms, and we compute the gain trial-wise, averaging them over paradigms. Only gains where the r2 of linear fit is > 0.8 are retained. Different colours indicate increasing mean stimulus, consistent with the rest of this document. We use a single filter (show in the previous figure) to compute all the gains. 

if ~exist('inst_gain','var')
	inst_gain = NaN*fA;
	inst_gain_r2 = NaN*fA;
	window_size = 5e2;
	for i = 1:width(PID)
		textbar(i,width(PID))
		for t = 5e3:10:7e3
			x = fp(t-window_size:t,i);
			y = fA(t-window_size:t,i);
			[temp,gof] = fit(x(:),y(:),'poly1');
			inst_gain(t,i) = temp.p1;
			inst_gain_r2(t,i) = gof.rsquare;
		end
	end
	inst_gain(inst_gain_r2<.8) = NaN;
end

timescale_firing = NaN(max(paradigm),1);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
ax(1) = subplot(1,3,1); hold on
ax(2) = subplot(1,3,2); hold on
for i = 1:max(paradigm)
	temp = nanmean(inst_gain(:,paradigm==i),2);
	time = 1e-3*(1:length(inst_gain));
	time = time(~isnan(temp));
	temp = temp(~isnan(temp));
	time = time - 5;

	% fit an exponential to get a timescale
	fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[mean(temp(1:5)) 1 1 1],'Lower',[min(temp(1:5)) 0 0 0 ],'Upper',[max(temp(1:5)) Inf Inf Inf]);
	ft = fittype('A*exp(-x/tau+c)+d','options',fo);
	ff = fit(time(:),temp(:),ft);
	timescale_firing(i) = ff.tau*1e3;

	plot(ax(1),time,temp,'+-','Color',c(i,:))
	% also plot the fit
	plot(ax(2),time,ff(time),'Color',c(i,:))
end
xlabel(ax(1),'Time since odour onset (s)')
ylabel(ax(1),'Gain (Hz/V)')
set(ax(1),'XScale','linear','yScale','log')
title(ax(1),'Gain as a function of time')

xlabel(ax(2),'Time since odour onset (s)')
ylabel(ax(2),'Gain (Hz/V)')
set(ax(2),'XScale','linear','yScale','log')
title(ax(2),'Exponential fits')

subplot(1,3,3), hold on
plot(1:max(paradigm),timescale_firing,'k+')
xlabel('Paradigm')
set(gca,'YLim',[0 200])
ylabel('Timescale of firing (ms)')
title('Timescale of gain control')
suptitle('Timescale of firing gain control')

prettyFig('plw',1.3,'lw',1.5,'fs',14,'FixLogX',true)
labelFigure

if being_published
	snapnow
	delete(gcf)
end

%% Timescales of change in LFP gain
% We now repeat this analysis, but for the LFP. 

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end

% extract filter from lowest dose
K_LFP = K*NaN;
K_LFP_Damon = K*NaN;
a_LFP = 15e3; % to avoid the long, 10s filter
z_LFP = 45e3; % ditto
for i = 1:width(PID)
	if (paradigm(i)==1)
		K_LFP_Damon(:,i) = backOutFilter(PID(a_LFP:z_LFP,i),filtered_LFP(a_LFP:z_LFP,i),'offset',100,'filter_length',700);
		temp = fitFilter2Data(PID(a_LFP:z_LFP,i),filtered_LFP(a_LFP:z_LFP,i),'offset',200,'filter_length',900,'reg',1);
		K_LFP(:,i) = temp(101:end-100);
	end
end

% find gains everywhere
[~,~,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a_LFP,'z',z_LFP);

% make linear predictions using the average filter
K_LFP_av = nanmean(K_LFP(:,paradigm==1),2);
LFP_pred = NaN*LFP;
for i = 1:width(PID)
	LFP_pred(:,i) = convolve(1e-3*(1:length(PID)),PID(:,i),K_LFP_av,filtertime);
end


if ~exist('inst_gain_LFP','var')
	inst_gain_LFP = NaN*fA;
	inst_gain_LFP_r2 = NaN*fA;
	window_size = 5e2;
	for i = 1:width(PID)
		textbar(i,width(PID))
		for t = 5e3:10:7e3
			x = LFP_pred(t-window_size:t,i);
			y = LFP(t-window_size:t,i);
			if ~any(isnan(x)) & ~any(isnan(y))
				[temp,gof] = fit(x(:),y(:),'poly1');
				inst_gain_LFP(t,i) = temp.p1;
				inst_gain_LFP_r2(t,i) = gof.rsquare;
			end
		end
	end
	inst_gain_LFP(inst_gain_LFP_r2 < .8) = NaN;
end


figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
temp = nanmean(K_LFP_Damon(:,paradigm==1),2);
temp = temp/-min(temp);
plot(filtertime,temp,'r')
temp = nanmean(K_LFP(:,paradigm==1),2);
temp = temp/-min(temp);
plot(filtertime,temp,'k')
legend({'Damon Filter','Srinivas Filter'})

timescale_LFP = NaN(max(paradigm),1);

ax(1) = subplot(2,2,2); hold on
ax(2) = subplot(2,2,3); hold on
for i = 1:max(paradigm)
	temp = nanmean(inst_gain_LFP(:,paradigm==i),2);
	time = 1e-3*(1:length(inst_gain_LFP));
	time = time(~isnan(temp));
	temp = temp(~isnan(temp));
	time = time - 5;

	% fit an exponential to get a timescale
	fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[mean(temp(1:5)) 1 1 1],'Lower',[min(temp(1:5)) 0 0 0 ],'Upper',[max(temp(1:5)) Inf Inf Inf]);
	ft = fittype('A*exp(-x/tau+c)+d','options',fo);
	ff = fit(time(:),temp(:),ft);
	timescale_LFP(i) = ff.tau*1e3;

	plot(ax(1),time,temp,'+-','Color',c(i,:))
	% also plot the fit
	plot(ax(2),time,ff(time),'Color',c(i,:))
end
xlabel(ax(1),'Time since odour onset (s)')
ylabel(ax(1),'Gain (Hz/V)')
set(ax(1),'XScale','linear','yScale','log')
title(ax(1),'Gain as a function of time')

xlabel(ax(2),'Time since odour onset (s)')
ylabel(ax(2),'Gain (Hz/V)')
set(ax(2),'XScale','linear','yScale','log')
title(ax(2),'Exponential fits')

subplot(2,2,4), hold on
plot(1:max(paradigm),timescale_firing,'k+')
xlabel('Paradigm')
set(gca,'YLim',[0 200])
ylabel('Timescale of firing (ms)')
title('Timescale of gain control')
suptitle('Timescale of firing gain control')

prettyFig('plw',1.3,'lw',1.5,'fs',14,'FixLogX',true)
labelFigure

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
% 
pFooter;
