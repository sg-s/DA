% Paper Figures
% makes all the figures for the paper
% 
% created by Srinivas Gorur-Shandilya at 12:57 , 21 January 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

% this determines which figures to do. 
fig1 = false; 	% how we determine gain
fig2 = false;	% weber-like gain control
fig3 = false;	% speed-gain tradeoff
fig4 = false;	% fast gain control
fig5 = false;	% fast gain control widely observed
fig6 = true; 	% natural stimuli
fig7 = false; 	% switching experiment
fig8 = false;	% LFP experiments
fig9 = false; 	% models to explain this

%    ######## ####  ######   ##     ## ########  ########       ##   
%    ##        ##  ##    ##  ##     ## ##     ## ##           ####   
%    ##        ##  ##        ##     ## ##     ## ##             ##   
%    ######    ##  ##   #### ##     ## ########  ######         ##   
%    ##        ##  ##    ##  ##     ## ##   ##   ##             ##   
%    ##        ##  ##    ##  ##     ## ##    ##  ##             ##   
%    ##       ####  ######    #######  ##     ## ########     ###### 

clearvars -except being_published fig*

if fig1

%% Figure 1: ORN gain can be estimated by measuring responses to Gaussian inputs
% Gaussian odorant inputs with short correlation times (A), elicit flickering responses in ORNs that track the odorant stimulus well (B). A linear filter K can be extracted from the odorant input and the firing rate output of the neuron (C). The slope of the residuals in a plot of the firing response vs. the linear prediction (D) is defined as the gain of the ORN in this stimulus paradigm. Here, we measure the ORN gain in the linear regime: the linear filter accounts for 96% of the variance in the ORN response (red line), and adding an output nonlinearity (dotted black line), only accounts for an additional 1% of the variance in the responses. The odorant used is ethyl acetate, stimulating the ab3A neuron. Shading in all plots shows the standard error of the mean. 

figure('outerposition',[0 0 800 700],'PaperUnits','points','PaperSize',[800 700]); hold on
axes_handles(1) = subplot(7,2,1:4); hold on
axes_handles(2) = subplot(7,2,5:8); hold on
axes_handles(3) = subplot(7,2,[9 11 13]); hold on
axes_handles(4) = subplot(7,2,[10 12 14]); hold on
% load cached data
load('MeanShiftedGaussians.mat')

% shorten paradigm names by throwing out 'MFC'
short_paradigm_names = paradigm_names;
for i = 1:length(paradigm_names)
	short_paradigm_names{i} = paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end);
end


% some global parameters
nbins = 50;
histx = [];
histy = [];
dt = 1e-3;
all_pid = [];


c = parula(length(paradigm_names)+1);

mean_pid = NaN(length(c),1);

% plot lowest dose stimulus
plot_these=find(strcmp(paradigm_names{1}, combined_data.paradigm));
plot_this = mean2(combined_data.PID(plot_these,:));
time = dt*(1:length(plot_this));
plot(axes_handles(1),time,plot_this,'Color',c(1,:))
ylabel(axes_handles(1),'Stimulus (V)')


% plot lowest dose response
plot_this = mean2(combined_data.fA(:,plot_these));
time = dt*(1:length(plot_this));
plot(axes_handles(2),time,plot_this,'Color',c(1,:))
ylabel(axes_handles(2),'ORN Response (Hz)')
set(axes_handles(1),'XLim',[15 55])
set(axes_handles(2),'XLim',[15 55])
xlabel(axes_handles(2),'Time (s)')

% load the data cut and processed
load('MSG_per_neuron.mat','MSG_data')
% if ~exist('MSG_data','var')
% 	load('MSG_per_neuron.mat','MSG_data')

% 	% back out all filters
% 	for i = 1:8
% 		for j = 1:13
% 			disp([i j])
% 			if width(MSG_data(i,j).stim) > 1
% 				this_stim = mean2(MSG_data(i,j).stim);
% 				this_stim = this_stim - mean(this_stim);
% 				this_stim = this_stim/std(this_stim);
% 				this_resp = mean2(MSG_data(i,j).resp);
% 				this_resp = this_resp - mean(this_resp);
% 				this_resp = this_resp/std(this_resp);
% 				[K,filtertime_full] = fitFilter2Data(this_stim,this_resp,'reg',1,'filter_length',1999,'offset',300);
% 				filtertime_full = filtertime_full;
% 				filtertime = (-200:800);
% 				K = interp1(filtertime_full,K,filtertime);
% 				MSG_data(i,j).K = K;
% 			end
% 		end
% 	end
% end
% save('MSG_per_neuron.mat','MSG_data')

% plot the filter for the lowest dose 
filtertime = (-200:800)*1e-3;
K = NaN*filtertime;
for i = 1:13 % get all the filters for the lowest dose
	K = [K; MSG_data(1,i).K];
end
K(1,:) = [];
err = std(K);
err = err/sqrt(width(K));
axes(axes_handles(3))
shadedErrorBar(filtertime,mean2(K),err,{'Color',c(1,:)})
set(axes_handles(3),'XLim',[min(filtertime) max(filtertime)])
xlabel(axes_handles(3),'Lag (s)')
ylabel(axes_handles(3),'Filter K (norm)')

% make linear predictions everywhere
% and also calculate the r2 of each -- this will be used as weights
for i = 1:8
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			this_stim = mean2(MSG_data(i,j).stim);
			this_resp = mean2(MSG_data(i,j).resp);
			MSG_data(i,j).fp = convolve(MSG_data(i,j).time,mean2(MSG_data(i,j).stim),[0 0 MSG_data(i,j).K(3:end)],filtertime) ;
			MSG_data(i,j).r2 = rsquare(MSG_data(i,j).fp,mean2(MSG_data(i,j).resp));
		end
	end
end


% plot linear prediction vs. data for the lowest dose. 
ss = 25;
y = mean2([MSG_data(1,:).resp]);
x = mean2([MSG_data(1,:).fp]);
plot(axes_handles(4),x(1:ss:end),y(1:ss:end),'.','Color',c(1,:));

ff= fit(x(~isnan(x)),y(~isnan(x)),'poly1');
l = plot(axes_handles(4),sort(x),ff(sort(x)),'r');
L = {};
L{1} = strcat('Gain=',oval(ff.p1),'Hz/V, r^2=',oval(rsquare(y,x)));

% also fit a hill function
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
minx = min(x);
x = x - minx;
clear p
p.     A= 42.9682;
p.     k= 1.0570;
p.     n= 2.4301;
p.offset= 2.9976;
[x,idx] = sort(x); y = y(idx);
l(2) = plot(axes_handles(4),x + minx,hill4(p,x),'k--');
L{2} = strcat('Hill fit, r^2=',oval(rsquare(hill4(p,x),y)));

xlabel(axes_handles(4),'K \otimes s')
ylabel(axes_handles(4),'Response (Hz)')


legend(l,L,'Location','northwest');


PrettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')

if being_published
	snapnow
	delete(gcf)
end


end


%        ######## ####  ######   ##     ## ########  ########     #######  
%        ##        ##  ##    ##  ##     ## ##     ## ##          ##     ## 
%        ##        ##  ##        ##     ## ##     ## ##                 ## 
%        ######    ##  ##   #### ##     ## ########  ######       #######  
%        ##        ##  ##    ##  ##     ## ##   ##   ##          ##        
%        ##        ##  ##    ##  ##     ## ##    ##  ##          ##        
%        ##       ####  ######    #######  ##     ## ########    ######### 

%% Figure 2: ORN gain decreases with increasing stimulus intensity, similar to the Weber-Fechner Law
% Odorant stimuli drawn from distributions with similar variances but increasing means (A) elicit ORN responses with decreasing variances (B). After extracting linear filters for all stimulus paradigms, a plot of the ORN response vs. the linear prediction (C) shows a systematic decrease in slope. Plotting lines to each of these clouds of points determines the neuron gain in each case. Neuron gain decreases with the mean stimulus (D). This decrease is well described by a power law with an exponent close to -1 (the Weber-Fechner Prediction). For comparison, a power law with the exponent fixed at -1 is also fitted (dashed line). 

if fig2

figure('outerposition',[0 0 800 700],'PaperUnits','points','PaperSize',[800 700]); hold on
axes_handles(5) = subplot(2,2,1); hold on;
axes_handles(6) = subplot(2,2,2); hold on;
axes_handles(7) = subplot(2,2,3); hold on;
axes_handles(8) = subplot(2,2,4); hold on;

% plot the stimulus distributions 
a = floor(15/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	plot_hist = (combined_data.PID(plot_these,a:z));
	[hy,hx]  = hist(plot_hist(:),50);
	hy = hy/length(plot_these);
	plot(axes_handles(5),hx,hy,'Color',c(i,:));
end

xlabel(axes_handles(5),'Stimulus (V)')
ylabel(axes_handles(5),'count')


% plot the response distributions 
a = floor(15/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	temp =  [MSG_data(i,:).resp];
	temp = mean2(temp);
	[hy,hx]  = hist(temp,50);
	plot(axes_handles(6),hx,hy,'Color',c(i,:));
end

xlabel(axes_handles(6),'Response (Hz)')
ylabel(axes_handles(6),'count')


% show gain changes for all paradigms -- average over neurons 
ss = 50;
for i = 1:8 % iterate over all paradigms 
	y = ([MSG_data(i,:).resp]);
	x = ([MSG_data(i,:).fp]);
	if ~isvector(x)
		x = mean2(x);
	end

	if ~isvector(y)
		y = mean2(y);
	end 

	plot(axes_handles(7),x(1:ss:end),y(1:ss:end),'.','Color',c(i,:))
end

xlabel(axes_handles(7),'K\otimes s')
ylabel(axes_handles(7),'Neuron Response (Hz)')


% compute gain changes on a per-neuron basis
gain = NaN(8,13);
mean_stim = NaN(8,13);
mean_resp = NaN(8,13);
gain_err = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			x = MSG_data(i,j).fp;
			if ~isvector(x)
				x = mean2(x);
			end
			if ~isvector(y)
				y = mean2(y);
			end 
			
			% trim NaNs again
			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));
			mean_resp(i,j) = mean(mean([MSG_data(i,j).resp]));

			% get the weights for the each
			temp = confint(temp);
			gain_err(i,j) = diff(temp(:,1))/2;
		end
	end	
end

% show gain changes -- gain vs. mean stimulus
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		errorbar(axes_handles(8),mean_stim(i,j),gain(i,j),gain_err(i,j),'+','Color',c(i,:));
	end
end

mean_stim = mean_stim(~isnan(mean_stim));
mean_resp = mean_resp(~isnan(mean_resp));
gain = gain(~isnan(gain));
gain_err = gain_err(~isnan(gain_err));


cf = fit(mean_stim(:),gain(:),'power1','Weights',1./gain_err);
set(axes_handles(8),'XScale','log','YScale','log','YLim',[1 200],'XLim',[.5 3.5])
% set(axes_handles(8),'XScale','linear','YScale','linear','YLim',[1 45],'XLim',[.5 3.5])
xlabel(axes_handles(8),'Mean Stimulus (V)')
ylabel(axes_handles(8),'Neuron Gain (Hz/V)')
l(1)=plot(axes_handles(8),sort(mean_stim),cf(sort(mean_stim)),'k');
r2 = rsquare(cf(mean_stim(:)),gain(:));
L = strcat('y = \alpha x^{\beta}, \beta= ',oval(cf.b), ',r^2 = ',oval(r2));

% fit a power law with exponent -1
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(mean_stim(:),gain(:),'power1',options);
l(2)=plot(axes_handles(8),sort(mean_stim),cf(sort(mean_stim)),'k--');
r2 = rsquare(cf(mean_stim(:)),gain(:));
legend(l,{L, strcat('y = \alpha x^{-1}, \alpha= ',oval(cf.a), ', r^2 = ',oval(r2))} )

PrettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')

if being_published
	snapnow
	delete(gcf)
end

end


%            ######## ####  ######   ##     ## ########  ########     #######  
%            ##        ##  ##    ##  ##     ## ##     ## ##          ##     ## 
%            ##        ##  ##        ##     ## ##     ## ##                 ## 
%            ######    ##  ##   #### ##     ## ########  ######       #######  
%            ##        ##  ##    ##  ##     ## ##   ##   ##                 ## 
%            ##        ##  ##    ##  ##     ## ##    ##  ##          ##     ## 
%            ##       ####  ######    #######  ##     ## ########     #######  


if fig3

%% Figure 3: ORNs speed up responses on increasing stimulus mean
% ORN gain may permit responses to occur with smaller delays, as ORNs need to integrate a smaller stimulus duration to respond. Speedup in responses can be estimated using the peak times of linear filters fit to increasing mean concentrations of odorant (colours, in A). ORN response speedups with increasing stimulus mean can also be estimated in a model-free manner using the spike-triggered average (STA, B). 

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on


% fit parametric filters to the raw, neuron-wise filters extracted earlier 
% for i = 1:8
% 	for j = 1:13
% 		if ~isempty(MSG_data(i,j).K)
% 			d.stimulus = MSG_data(i,j).K(200:end);
% 			d.response = MSG_data(i,j).K(200:end);
% 			for k = 1:5
% 				MSG_data(i,j).p = FitModel2Data(@FitFilter,d,MSG_data(i,j).p);
% 			end
% 		end
% 	end
% end

% compute peak locations of all these filters
clear l 
l = zeros(8,1);
peak_loc_K = NaN(8,13);
subplot(1,3,1), hold on
for i = 1:8
	for j = 1:13
		if ~isempty(MSG_data(i,j).K)
			K2 = FitFilter(MSG_data(i,j).K(200:end),MSG_data(i,j).p);
			filtertime = 1e-3*(1:length(K2));
			l(i)=plot(filtertime,K2,'Color',c(i,:));
			[~,loc] = max(K2);
			peak_loc_K(i,j) = filtertime(loc);
		end
	end
end


set(gca,'XLim',[-.01 .5])
xlabel('Lag (s)')
ylabel('Filter (norm)')
L = paradigm_names;
for i = 1:length(L)
	L{i} = L{i}(strfind(L{i},'-')+1:end);
end
legend(l,L)

before = 1e4;
after = 5e3;
% compute STA for all the data
% allfiles = dir('/local-data/DA-paper/fast-flicker/orn/*.mat');
% allSTA = [];
% mean_stim = []; 
% dil = [];
% for i= 1:length(allfiles)
% 	load(['/local-data/DA-paper/fast-flicker/orn/' allfiles(i).name])
% 	for j = 1:length(spikes)
% 		if length(spikes(j).A) > 1 && ~strcmp(ControlParadigm(j).Name,'end')
% 			these_spikes = spikes(j).A(:,35e4:55e4);
% 			this_stim = data(j).PID(:,35e4:55e4);

% 			if isfield(spikes,'discard')
% 				rm_this = find(spikes(j).discard);
% 				rm_these_spikes = [];
% 				if ~isempty(rm_this)
% 					this_stim(rm_this,:) = [];
% 					for k = 1:length(rm_this)
% 						if rm_this(k) > width(these_spikes)
% 						else
% 							rm_these_spikes = [rm_these_spikes k];
% 						end
% 					end
% 					these_spikes(rm_these_spikes,:) = [];
% 				end
% 			end

% 			this_STA = STA(these_spikes,this_stim,before,after);
% 			allSTA = [allSTA this_STA];
% 			mean_stim = [mean_stim; mean(this_stim,2)];
% 			this_dil = str2double(ControlParadigm(j).Name(strfind(ControlParadigm(j).Name,'-')+1:strfind(ControlParadigm(j).Name,'%')-1));
% 			dil = [dil;this_dil*ones(width(this_stim),1)  ];
% 		end
% 	end
% end
% % save this
% save('allSTA.mat','allSTA','dil','mean_stim');




% plot the STA
old_mean_stim = mean_stim;
load('allSTA.mat','allSTA','dil','mean_stim');
subplot(1,3,2), hold on
udil = unique(dil);
for i = 1:8
	temp = (allSTA(:,dil == udil(i)));
	for j = 1:width(temp)
		temp(:,j) = temp(:,j)/max(temp(:,j));
	end
	t = -after:before;
	t = t*1e-4; 
	errorShade(t,flipud(mean2(temp)),flipud(sem(temp)),'Color',c(i,:));
end
set(gca,'XLim',[-.2 .5])
xlabel('Lag (s)')
ylabel('STA (norm)')

% find the peak of each STA
peak_STA = NaN*mean_stim;
for i = 1:length(mean_stim)
	[~,loc] = max(allSTA(:,i));
	peak_STA(i) = (1e4 - loc)/10;
end

x = NaN(8,1); y = x; ex = x; ey = x;
peak_STA(peak_STA>900) = NaN; % throw out some bullshit values
for i = 1:8
	x(i) = mean2(mean_stim(dil == udil(i)));
	ex(i) = sem(mean_stim(dil == udil(i)));
	y(i) = nanmean(peak_STA(dil == udil(i)));
	ey(i) = sem(nonnans(peak_STA(dil == udil(i))));
end

subplot(1,3,3), hold on
clear l
l(1) = errorbar(x,y,ey,'k')
peak_loc_K = peak_loc_K(~isnan(peak_loc_K));
l(2) = plot(old_mean_stim(:),peak_loc_K(:)/dt,'r+');

% calculate Spearman's rho (http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
s2 = spear(old_mean_stim,peak_loc_K);
s1 = spear(mean_stim,peak_STA);

legend(l,{strcat('STA, \rho=',oval(s1)), strcat('Filter, \rho=',oval(s2))})
set(gca,'YLim',[-10 130])
ylabel('Peak time (ms)')
xlabel('Mean Stimulus (V)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

end

%           ######## ####  ######   ##     ## ########  ########    ##        
%           ##        ##  ##    ##  ##     ## ##     ## ##          ##    ##  
%           ##        ##  ##        ##     ## ##     ## ##          ##    ##  
%           ######    ##  ##   #### ##     ## ########  ######      ##    ##  
%           ##        ##  ##    ##  ##     ## ##   ##   ##          ######### 
%           ##        ##  ##    ##  ##     ## ##    ##  ##                ##  
%           ##       ####  ######    #######  ##     ## ########          ##  


if fig4

%% Figure 4: ORNs can change gain on a fast time scale 
% Sensors that control gain in response to changing stimulus statistics must do so relevant timescales: a gain control mechanism that is too slow risks saturation or insensitivity, while a gain control mechanism that is too quick to change leads to fluctuating responses. To mimic an odor signal that changes its statistics rapidly, we recorded ORN responses to a flickering stimulus with a large variance (A-B).  The linear filter extracted from this data (C) convolved with the stimulus yields a linear prediction (D). Comparing the responses (y-axis) to the linear prediction (x-axis) when the stimulus is locally low (green) vs. locally high (red) shows that neuron gain at these two times is significantly different (E). Within this single stimulus presentation, gain varies inversely with the recent stimulus (F). Repeating the analysis over a range of time scales identifies a region of sub-second timescales over which ORNs significantly change their gain in a stimulus-driven manner. The odour used is ethyl acetate, and recordings are from ab3A. The vertical line in G indicates the history length used in E. 

clearvars -except being_published fig*

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);
B_spikes = vertcat(B_spikes,spikes(4).B);

% A spikes --> firing rate
hash = DataHash(full(all_spikes));
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

fig_handle=figure('Units','pixels','outerposition',[100 100 1240 720],'Color','w','PaperUnits','points','PaperSize',[1240 720],'Toolbar','none');
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position',[70 490 402.5 105]);
axes_handles(2)=axes('Units','pixels','Position',[542.5 437.5 105 105]);
axes_handles(3)=axes('Units','pixels','Position',[735 560 472.5 105]);
axes_handles(4)=axes('Units','pixels','Position',[735 437.5 472.5 105]);
axes_handles(5)=axes('Units','pixels','Position',[70 70 262.5 262.5]);
axes_handles(6)=axes('Units','pixels','Position',[437.5 70 262.5 262.5]);
axes_handles(7)=axes('Units','pixels','Position',[805 70 350 262.5]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end

% set up a colour map
c = parula(8);


% plot stimulus
plot(axes_handles(1),tA,mean2(PID),'k')
set(axes_handles(1),'XLim',[10 60])
ylabel(axes_handles(1),'Stimulus (V)')

% plot response
plot(axes_handles(3),tA,mean2(fA),'k')
set(axes_handles(3),'XLim',[10 60],'XTickLabel',{});
ylabel(axes_handles(3),'Firing Rate (Hz)')

% extract Linear model for each trial 
K = [];
for i = [3:10 13:20]
	[this_K, ~, filtertime_full] = FindBestFilter(PID(:,i),fA(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	K = [K;this_K];
end

% plot linear filter
axes(axes_handles(2))
plot_this = mean2(K);
err = std(K)/sqrt(width(K));
shadedErrorBar(filtertime,plot_this,err,{'Color',c(3,:)})
xlabel(axes_handles(2),'Lag (s)')
ylabel(axes_handles(2),'Filter K')


% make a linear prediction using a filter fit to the mean data (this is almost exactly the same)
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);
K = K/max(K);
fp = convolve(tA,mean2(PID),K,filtertime);
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;


% plot prediction and prediction quality
l=plot(axes_handles(4),tA,fp,'Color',c(3,:));
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2)))
ylabel(axes_handles(4),'K\otimes stimulus (Hz)')


% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles([5 7]);

hl_min = .1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',mean2(fA),'prediction',fp,'stimulus',mean2(PID),'time',tA,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',1,'engine',@GainAnalysis5);
set(axes_handles(7),'XLim',[.09 11]) % to show .1 and 10 on the log scale

% show the p-value
axes(axes_handles(5))
text(10,60,'p < 0.01')

% plot gain vs preceding stimulus
[x,y] = MakeFig6G(mean2(PID),mean2(fA),fp,500);
gain_time = mean(diff(tA))*(1:length(x));
rm_this = (isnan(x) | isnan(y));
x(rm_this) = [];
y(rm_this) = [];
gain_time(rm_this) = [];
ss = 50;
plot(axes_handles(6),x(1:ss:end),y(1:ss:end),'k.')
xlabel(axes_handles(6),'Stimulus in preceding 500ms')
ylabel(axes_handles(6),'Relative gain')

% fit a Weber-scaling to these points
[x,idx] = sort(x);
y = y(idx);
% get equal bins in x
my = NaN(100,1); mx = my;
a = 1;
ss = floor(length(x)/length(my));
for i = 1:length(my)-1
	z = a+ss;
	my(i) = max(y(a:z));
	mx(i) = mean(x(a:z));
	a = z;
end
my(end) = []; mx(end) = [];
mx(1:20) = []; my(1:20) = [];

fo = fitoptions('rat01');
fo.StartPoint = [.4 -.08];
ff = fit(mx(:),my(:),'rat01',fo);
l = plot(axes_handles(6),.17:.01:max(x),ff(.17:.01:max(x)),'r');
legend(l,'Weber Envelope');

% plot gain
gain = y;
% plot(axes_handles(7),gain_time,gain,'k')
% ylabel(axes_handles(7),'Relative Gain')
% set(axes_handles(7),'XLim',[10 60],'YLim',[0 7])
% xlabel(axes_handles(7),'Time (s)')


% link some axes
linkaxes(axes_handles([3:4]))

% fix some labels
xlabel(axes_handles(4),'Time (s)')
xlabel(axes_handles(1),'Time (s)')
ylabel(axes_handles(7),'Relative Gain')
set(axes_handles(7),'YScale','log','YTick',[0.5 1 2],'YLim',[0.4 3.5])
set(axes_handles(4),'YLim',[0 100])
ylabel(axes_handles(5),'Firing Rate (Hz)')
xlabel(axes_handles(5),'K\otimes stimulus (Hz)')
title(axes_handles(5),'')

% Indicate regions of high and low stimulus on the stimulus
hl = round(history_lengths(9)*1e3);
shat = ComputeSmoothedStimulus(mean2(PID),hl);

n = floor(sum(~isnan(mean2(fA)))*.33);
shat(1:hl) = Inf; % the initial segment where we can't estimate shat is excluded
shat(isnan(shat)) = Inf;
shat(isnan(mean2(fA))) = Inf;
[~, t_low] = sort(shat,'ascend');
t_low = t_low(1:n); % this is an index
t_low = tA(t_low); % t_low is now a time. 
 
shat = ComputeSmoothedStimulus(mean2(PID),hl);
shat(1:hl) = -Inf;
shat(isinf(shat)) = -Inf;
shat(isnan(mean2(fA))) = -Inf;
[~, t_high] = sort(shat,'descend');
t_high = t_high(1:n);
t_high  = tA(t_high);

plot(axes_handles(1),t_low,1+0*t_low,'g.')
plot(axes_handles(1),t_high,1+0*t_low,'r.')


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end

end

%   ######## ####  ######   ##     ## ########  ########    ######## 
%   ##        ##  ##    ##  ##     ## ##     ## ##          ##       
%   ##        ##  ##        ##     ## ##     ## ##          ##       
%   ######    ##  ##   #### ##     ## ########  ######      #######  
%   ##        ##  ##    ##  ##     ## ##   ##   ##                ## 
%   ##        ##  ##    ##  ##     ## ##    ##  ##          ##    ## 
%   ##       ####  ######    #######  ##     ## ########     ######  


if fig5

%% Figure 4: Gain Control is widely observed
% The previous result showed fast gain control using ethyl acetate, a volatile odor, in ab3A, an antennal ORN. To determine if fast gain control is widely observed, and not a peculiarity of the odor-receptor combination used, I reanalyzed a published data set of ORN responses to flickering odor stimuli (from Martelli et al. J. Neuro. 2013). The data set consists of six different odors and two ORN types: one in the antenna and one in the maxillary palp. My analysis of fast gain control is self-consistent across experimental replicates (A-C), and shows that fast gain control is observed for different odors (D-F) and in two different neurons (G-I).  In each row, the first column shows the linear filter extracted from the stimulus and the response. The second column shows the principal components fit to the highest (lowest) 1/3 of stimulus filtered over some history length in red (green). The third column shows how relative gain varies at times when the filtered stimulus in in the top 1/3 (green) or bottom 1/3 (red), for various history lengths. Only significant (p<0.01) points are shown. In every dataset, significant gain control is observed on a sub-second time scale, and the gain control always acts to amplify responses to recently weak stimuli and suppress responses to recently large stimuli (red dots always below green dots in third column). 

clearvars -except being_published fig*
load('CMData_Gain.mat')
load('CM_Data_filters.mat')
combined_data_file = ('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat');
load(combined_data_file)
filtertime = -200:700;
filtertime = filtertime*1e-3;
history_lengths = logspace(log10(.200),1,30);


s = 860;
figure('outerposition',[0 0 s s],'PaperUnits','points','PaperSize',[s s]); hold on, clear s
axes_handles = [];
for i = 1:9
	axes_handles(i) = subplot(3,3,i); hold on
end

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  

% add a new field with correlation times
for i = 1:length(data)
	if any(strfind(data(i).original_name,'30ms'))
		data(i).corr_time = 30;
	end
	if any(strfind(data(i).original_name,'50ms'))
		data(i).corr_time = 50;
	end
	if any(strfind(data(i).original_name,'100ms'))
		data(i).corr_time = 100;
	end
end


% first row: experimental replicates
do_these = [2 7 9 13 15 16];

for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	plot(axes_handles(1),filtertime,mean2(K))
end
clear ph
ph(3:4) = axes_handles(2:3);


for i = do_these
	time = 1e-3*(1:length(mean2(data(i).PID)));
	stimulus = mean2(data(i).PID);
	prediction = mean2(data(i).LinearFit);
	response = mean2(data(i).fA);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	% fix the gain to be exactly 1
	x = prediction;
	y = response;
	rm_this = isnan(x) | isnan(y) ;
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	% ignore very low responses
	y(y<5) = NaN;

	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);
end


% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')


%    ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%    ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%    ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%    ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%    ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
   
%    ##    ## ######## ##     ## ########   #######  ##    ##  ######  
%    ###   ## ##       ##     ## ##     ## ##     ## ###   ## ##    ## 
%    ####  ## ##       ##     ## ##     ## ##     ## ####  ## ##       
%    ## ## ## ######   ##     ## ########  ##     ## ## ## ##  ######  
%    ##  #### ##       ##     ## ##   ##   ##     ## ##  ####       ## 
%    ##   ### ##       ##     ## ##    ##  ##     ## ##   ### ##    ## 
%    ##    ## ########  #######  ##     ##  #######  ##    ##  ######  


do_these = [11 17];
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(7),filtertime,mean2(K))];
end
legend(l,{'pb1A','ab3A'})


clear ph
ph(3:4) = axes_handles(8:9);
for i = do_these
	time = 1e-3*(1:length(mean2(data(i).PID)));
	stimulus = mean2(data(i).PID);
	prediction = mean2(data(i).LinearFit);
	response = mean2(data(i).fA);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	% fix the gain to be exactly 1
	x = prediction;
	y = response;
	rm_this = isnan(x) | isnan(y) ;
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	% ignore very low responses
	y(y<5) = NaN;

	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);
end



% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')



%     ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%     ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%     ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%     ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%     ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
    
%      #######  ########   #######  ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##    ## 
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##       
%     ##     ## ##     ## ##     ## ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##   ##         ## 
%     ##     ## ##     ## ##     ## ##     ## ##    ##  ##    ## 
%      #######  ########   #######   #######  ##     ##  ######  


do_these = [7 8 10 14 17 18];
odours = {'1but','1o3ol','dsucc','2ac','2but','5ol'};
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(4),filtertime,mean2(K))];
end
legend(l,odours)

clear ph
ph(3:4) = axes_handles(5:6);
for i = do_these
	time = 1e-3*(1:length(mean2(data(i).PID)));
	stimulus = mean2(data(i).PID);
	prediction = mean2(data(i).LinearFit);
	response = mean2(data(i).fA);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	% fix the gain to be exactly 1
	x = prediction;
	y = response;
	rm_this = isnan(x) | isnan(y) ;
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	% ignore very low responses
	y(y<5) = NaN;

	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);
end


% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')

% clean up -- remove the scatter points
for j = [2 5 8]
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		if strcmp(get(h(i),'Marker'),'.')
			rm_this = [rm_this i];
		end
	end
	delete(h(rm_this))
end

% remove the line indicating the example history plot
for j = [3 6 9]
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		try
			if  strcmp(get(h(i),'LineStyle'),'-.')
				rm_this = [rm_this i];
			end
		catch
		end
	end
	delete(h(rm_this))

	% remove the lines where the data isn't significant
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		try
			if  strcmp(get(h(i),'Type'),'line')
				rm_this = [rm_this i];
			end
		catch
		end
	end
	delete(h(rm_this))

	% make all the scatter plots smaller
	h=get(axes_handles(j),'Children');
	for i = 1:length(h)
		 set(h(i),'SizeData',256)
	end

	set(axes_handles(j),'XLim',[.1 10],'YLim',[.6 1.45])
end


% cosmetics
xlabel(axes_handles(1),'Lag (s)')
xlabel(axes_handles(4),'Lag (s)')
title(axes_handles(1),'Filters')
title(axes_handles(4),'Filters')
title(axes_handles(7),'Filters')
ylabel(axes_handles(3),'Relative Gain')
ylabel(axes_handles(6),'Relative Gain')
ylabel(axes_handles(9),'Relative Gain')
title(axes_handles(2),'')
title(axes_handles(5),'')
title(axes_handles(8),'')



PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

ylabel(axes_handles(1),'Exp. Replicates','FontSize',20)
ylabel(axes_handles(7),'Diff. ORNs','FontSize',20)
ylabel(axes_handles(4),'Diff. odors','FontSize',20)

if being_published
	snapnow
	delete(gcf)
end

end

%          ######## ####  ######   ##     ## ########  ########     #######  
%          ##        ##  ##    ##  ##     ## ##     ## ##          ##     ## 
%          ##        ##  ##        ##     ## ##     ## ##          ##        
%          ######    ##  ##   #### ##     ## ########  ######      ########  
%          ##        ##  ##    ##  ##     ## ##   ##   ##          ##     ## 
%          ##        ##  ##    ##  ##     ## ##    ##  ##          ##     ## 
%          ##       ####  ######    #######  ##     ## ########     #######  



%% Figure 6: ORNs change gain rapidly to naturalistic stimuli. 
% 

clearvars -except being_published fig*

if fig6


load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;
B_spikes = spikes(2).B;


% A spikes --> firing rate
hash = DataHash(full(all_spikes));
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

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,3,1:3), hold on
errorShade(tA(1:10:end),mean2(PID(1:10:end,:)),sem(PID(1:10:end,:)),'Color',[0 0 0]);
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(2,3,4), hold on
for i = 1:width(PID)
	[y,x] = histcounts(PID(:,i),300);x(1) = [];
	plot(x,y,'Color',[.5 .5 .5])
end
set(gca,'XScale','log','YScale','log','XLim',[.1 10])
xlabel('Stimulus (V)')
ylabel('count')

% make a linear filter
[K, filtertime_full] = fitFilter2Data(mean2(PID),mean2(fA),'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean2(PID),K,filtertime);

% correct for trivial scaling
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;

shat = ComputeSmoothedStimulus(mean2(PID),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

% make the output analysis plot
subplot(2,3,5), hold on
ss = 1;
cc = parula(100);
c= cc(shat,:);
scatter(fp(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')
xlabel('Linear Prediction (Hz)')
ylabel('Actual response (Hz)')

% plot gain vs stimulus for all these whiffs
subplot(2,3,6), hold on
fp = convolve(tA,mean2(PID),K,filtertime);
shat = ComputeSmoothedStimulus(mean2(PID),500);

% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = ComputeOnsOffs(R>10);
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

errorbar(mean_stim,gain,gain_err,'k+')

ff = fit(mean_stim(:),gain(:),'power1','Weights',1./gain_err);

l = plot(sort(mean_stim),ff(sort(mean_stim)),'r');
L = strcat('y=\alpha x^{\beta}, \beta=',oval(ff.b),'. r^2=',oval(rsquare(ff(mean_stim),gain)));
legend(l,L)
xlabel('Mean Stimulus in preceding 500ms (V)')
ylabel('Gain (Hz/V)')

PrettyFig('fs=12;');

if being_published
	snapnow
	delete(gcf)
end

end


%      ##     ##  #######  ########  ######## ##        ######  
%      ###   ### ##     ## ##     ## ##       ##       ##    ## 
%      #### #### ##     ## ##     ## ##       ##       ##       
%      ## ### ## ##     ## ##     ## ######   ##        ######  
%      ##     ## ##     ## ##     ## ##       ##             ## 
%      ##     ## ##     ## ##     ## ##       ##       ##    ## 
%      ##     ##  #######  ########  ######## ########  ######  

if fig9


%% Figure 5: Models to explain observed phenomena 

figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
clear axes_handles
for i = 1:10
	axes_handles(i) = subplot(2,5,i); hold on
end

% first column shows fit quality
for i = 1:5:10
	plot(axes_handles(i),[0 5],[0 5],'k--')
	xlabel(axes_handles(i),'(P_{S}/P_{N})^{1/2}','interpreter','tex')
	ylabel(axes_handles(i),'(P_{S}/P_{R})^{1/2}','interpreter','tex')
end

% Label the columns
title(axes_handles(1),'Fit Quality')
title(axes_handles(2),'Weber Scaling')
title(axes_handles(3),'Response Speedup')
title(axes_handles(4),'Pulse Speedup')
title(axes_handles(5),'Fast Gain Control')



%       ##       ##    ##    ##     ##  #######  ########  ######## ##       
%       ##       ###   ##    ###   ### ##     ## ##     ## ##       ##       
%       ##       ####  ##    #### #### ##     ## ##     ## ##       ##       
%       ##       ## ## ##    ## ### ## ##     ## ##     ## ######   ##       
%       ##       ##  ####    ##     ## ##     ## ##     ## ##       ##       
%       ##       ##   ###    ##     ## ##     ## ##     ## ##       ##       
%       ######## ##    ##    ##     ##  #######  ########  ######## ######## 

load('LVF_data.mat')
tA = 1e-3*(1:length(data(1).stim));
% extract filters for each neuron
for i = 1:length(data)
	[this_K, ~, filtertime_full] = FindBestFilter(mean2(data(i).stim),mean2(data(i).resp),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*1e-3;
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	data(i).K = this_K;
end

% make linear predictions and extract non-linearities 
for i = 1:length(data)
	fp = convolve(tA,mean2(data(i).stim),data(i).K,filtertime);
	R = mean2(data(i).resp);
	temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
	fp = fp*temp.p1;
	fp = fp+temp.p2;
	data(i).LinearFit = fp;
	data(i).LNFit = hill(p_LN(i),data(i).LinearFit);

	% show fit quality of DA model for Large Variance Flicker
	[qx, qy] = GeffenMeister(data(i).resp,data(i).LNFit);
	plot(axes_handles(1),qx,qy,'ko')
end

% show gain analysis -- da model
ph = []; ph(3) = axes_handles(5);
history_lengths = 0.4890;
time = 1e-3*(1:length(data(1).LinearFit));
p=GainAnalysisWrapper2('response',mean2([data.resp]),'prediction',mean2([data.LNFit]),'stimulus',mean2([data.stim]),'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(1),'engine',@GainAnalysis5);


% show the p-value
axes(axes_handles(5))
if p(1) == 0
	text(10,60,'p<0.01')
else
	text(10,60,strkat('p = ',oval(p(1))))
end

%    ##       ##    ##    ########  ##     ## ##        ######  ########  ######  
%    ##       ###   ##    ##     ## ##     ## ##       ##    ## ##       ##    ## 
%    ##       ####  ##    ##     ## ##     ## ##       ##       ##       ##       
%    ##       ## ## ##    ########  ##     ## ##        ######  ######    ######  
%    ##       ##  ####    ##        ##     ## ##             ## ##             ## 
%    ##       ##   ###    ##        ##     ## ##       ##    ## ##       ##    ## 
%    ######## ##    ##    ##         #######  ########  ######  ########  ######  



load('PulseData.mat')

% average filters
K = mean2(vertcat(data.K));
K = K/max(K);

% create the data to be fit
% clear d
% time = 1e-3*(1:length(PulseData(1).stim(:,1)));
% for i = 1:length(PulseData)
% 	d(i).stimulus = convolve(time,mean2(PulseData(i).stim),K,filtertime);
% 	d(i).response = mean2(PulseData(i).resp);
% end

clear p
p.A =  225.1576;
p.k =  382.2206;
p.n =  0.3952;

% this stores quality metrics
qx = NaN(length(PulseData),1);
qy = NaN(length(PulseData),1);

% make predictions everywhere
for i = 1:length(PulseData)
	PulseData(i).pred = NaN*PulseData(i).stim;
	for j = 1:width(PulseData(i).stim)
		PulseData(i).pred(:,j) = convolve(time,PulseData(i).stim(:,j),K,filtertime);
		PulseData(i).pred(:,j) = hill(PulseData(i).pred(:,j),p);
	end
	[qx(i), qy(i)] = GeffenMeister(PulseData(i).resp,PulseData(i).pred);
end


% show this on the plot
plot(axes_handles(1),qx,qy,'k+')

% show the timing data
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).resp(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
resp_time_data = resp_half_time - stim_half_time;
foreground_stim(c:end) = [];

% clean up
rm_this = foreground_stim < 1e-2 | resp_time_data < 0 | resp_time_data > 100;
foreground_stim(rm_this) = [];
resp_time_data(rm_this) = [];

plot(axes_handles(4),foreground_stim,resp_time_data,'k+')
set(axes_handles(4),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
xlabel(axes_handles(4),'Mean Stimulus (V)')
ylabel(axes_handles(4),'\tau_{ORN}-\tau_{PID} (ms)','interpreter','tex')

% show the timing data for the LN model
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).pred(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
foreground_stim(c:end) = [];

resp_time = resp_half_time - stim_half_time;

% clean up exactly like we did the data
foreground_stim(rm_this) = [];
resp_time(rm_this) = [];

l=plot(axes_handles(4),foreground_stim,resp_time,'r+');

% show the spearman rho with the data
s=rsquare(resp_time,resp_time_data);
set(axes_handles(4),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
legend(l,strcat('r^2=',oval(s)))

% ##       ##    ##    ##      ## ######## ########  ######## ########  
% ##       ###   ##    ##  ##  ## ##       ##     ## ##       ##     ## 
% ##       ####  ##    ##  ##  ## ##       ##     ## ##       ##     ## 
% ##       ## ## ##    ##  ##  ## ######   ########  ######   ########  
% ##       ##  ####    ##  ##  ## ##       ##     ## ##       ##   ##   
% ##       ##   ###    ##  ##  ## ##       ##     ## ##       ##    ##  
% ######## ##    ##     ###  ###  ######## ########  ######## ##     ## 

load('MSG_per_neuron.mat','MSG_data')
load('LN_Fit_to_MSG.mat','p')
% fit a LN model for each neuron
% clear p
% for i = 1:13
% 	clear d
% 	c = 1;
% 	for j = 1:8
% 		if width(MSG_data(j,i).stim)>1
% 			d(c).stimulus = mean2(MSG_data(j,i).stim);
% 			d(c).response = mean2(MSG_data(j,i).resp);
% 			c = c+1;
% 		end
% 	end
% 	if exist('d','var')
% 		p(i) = FitModel2Data(@pLNModel,d,'p0',p(i),'nsteps',2);
% 	end
% end

% compute all the model predictions and the fit qualities 
qx = NaN(8,13);
qy = NaN(8,13);
for i = 1:8 % there are 8 paradimgs in the MSG data
	for j = 1:13 % there are 13 neurons
		if ~isempty(p(j).A) && width(MSG_data(i,j).stim) > 1
			for k = 1:width(MSG_data(i,j).stim)
				MSG_data(i,j).fp_LN(:,k) = pLNModel(MSG_data(i,j).stim(:,k),p(j));
			end
			% compute the geffen-Meister error metric
			[qx(i,j) qy(i,j)] = GeffenMeister(MSG_data(i,j).resp,MSG_data(i,j).fp_LN);
		end
	end
end

% plot on the error plot, and colour code by stimulus
c = parula(9);
for i = 1:8
	this_qx = qx(i,~isnan(qx(i,:)));
	this_qy = qy(i,~isnan(qy(i,:)));
	plot(axes_handles(1),this_qx,this_qy,'x','Color',c(i,:))
end


% get the gain scaling factor -- this is to correct for the fact that we normalised the filter heights everywhere
x = mean2(horzcat(MSG_data(1,:).stim));
y = mean2(horzcat(MSG_data(1,:).resp));
gsf = std(y)/std(x);


% show the gain plot again, and also the gain for the da model 
gain = NaN(8,13);
gain_DA = NaN(8,13);
mean_stim = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			x = MSG_data(i,j).fp;
			x_da = MSG_data(i,j).fp_LN;
			if ~isvector(x)
				x = mean2(x);
			end
			if ~isvector(y)
				y = mean2(y);
			end 
			if ~isvector(x_da)
				x_da = mean2(x_da);
			end 


			% trim NaNs again
			rm_this = isnan(x) | isnan(y) | isnan(x_da);
			x(rm_this) = [];
			x_da(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

			temp=fit(x(:),x_da(:),'poly1');
			gain_DA(i,j) = temp.p1;

			% get the units of gain right
			gain(i,j) = gain(i,j)*gsf;
			gain_DA(i,j) = gain_DA(i,j)*gsf;

		end
	end	
end
plot(axes_handles(2),mean_stim(:),gain(:),'kx')
l=plot(axes_handles(2),mean_stim(:),gain_DA(:),'rx');
xlabel(axes_handles(2),'Mean Stimulus (V)')
ylabel(axes_handles(2),'Neuron Gain (Hz/V)')
set(axes_handles(2),'XScale','log','YScale','log')

% show the spearman rho 
a = gain(:);
b = gain_DA(:);
rm_this = isnan(a) | isnan(b);
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));


% ##       ##    ##     ######  ########  ######## ######## ########  ##     ## ########  
% ##       ###   ##    ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## 
% ##       ####  ##    ##       ##     ## ##       ##       ##     ## ##     ## ##     ## 
% ##       ## ## ##     ######  ########  ######   ######   ##     ## ##     ## ########  
% ##       ##  ####          ## ##        ##       ##       ##     ## ##     ## ##        
% ##       ##   ###    ##    ## ##        ##       ##       ##     ## ##     ## ##        
% ######## ##    ##     ######  ##        ######## ######## ########   #######  ##        
   


peak_loc_data = NaN(8,13);
peak_loc_DA = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			this_resp = MSG_data(i,j).resp;
			this_stim = MSG_data(i,j).stim;
			this_pred = MSG_data(i,j).fp_LN;
			if ~isvector(this_resp)
				this_resp = mean2(this_resp);
			end
			if ~isvector(this_stim)
				this_stim = mean2(this_stim);
			end
			if ~isvector(this_pred)
				this_pred = mean2(this_pred);
			end

			a = this_resp - mean(this_resp);
			b = this_stim - mean(this_stim);
			a = a/std(a);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			x = x/max(x);
			[~,loc] = max(x);
			peak_loc_data(i,j) = t(loc);

			% throw out the first 5 seconds
			this_pred(1:5e3) = [];
			b(1:5e3) = [];
			a = this_pred - mean(this_pred);
			a = a/std(a);
			x = xcorr(a,b); % positive peak means a lags b
			x = x/max(x);
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			[~,loc] = max(x);
			peak_loc_DA(i,j) = t(loc);
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

		end
	end
end 

plot(axes_handles(3),mean_stim(:),1e3*peak_loc_data(:),'kx')
l=plot(axes_handles(3),mean_stim(:),1e3*peak_loc_DA(:),'rx');
set(axes_handles(3),'XLim',[-.1 4],'YLim',[5 120])
xlabel(axes_handles(3),'Mean Stimulus (V)')
ylabel(axes_handles(3),'Peak xcorr. (ms)')

a = peak_loc_data(:);
b = peak_loc_DA(:);
rm_this = isnan(a) | isnan(b) | a < 0 | b < 0;
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));

return



%       ########     ###               ##     ##  #######  ########  ######## ##       
%       ##     ##   ## ##              ###   ### ##     ## ##     ## ##       ##       
%       ##     ##  ##   ##             #### #### ##     ## ##     ## ##       ##       
%       ##     ## ##     ##            ## ### ## ##     ## ##     ## ######   ##       
%       ##     ## #########            ##     ## ##     ## ##     ## ##       ##       
%       ##     ## ##     ##            ##     ## ##     ## ##     ## ##       ##       
%       ########  ##     ##            ##     ##  #######  ########  ######## ######## 


% DA Model
load('LVF_data.mat')
% generate predictions
for i = 1:length(data)
	data(i).fp = DAModelv2(mean2(data(i).stim),p(i));
	% censor the first 10 seconds
	resp = data(i).resp(1e4:end,:);
	fp = data(i).fp(1e4:end);
	% show fit quality of DA model for Large Variance Flicker
	[qx, qy] = GeffenMeister(resp,fp);
	plot(axes_handles(6),qx,qy,'ko')
end


% show gain analysis -- da model
ph = []; ph(3) = axes_handles(10);
history_lengths = 0.4890;
time = 1e-3*(1:length(data(1).fp));
[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',mean2([data.resp]),'prediction',mean2([data.fp]),'stimulus',mean2([data.stim]),'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(1),'engine',@GainAnalysis5);


% show the p-value
axes(axes_handles(10))
text(10,60,strkat('p = ',oval(p(1))))



%     ########     ###             ########  ##     ## ##        ######  ########  ######  
%     ##     ##   ## ##            ##     ## ##     ## ##       ##    ## ##       ##    ## 
%     ##     ##  ##   ##           ##     ## ##     ## ##       ##       ##       ##       
%     ##     ## ##     ##          ########  ##     ## ##        ######  ######    ######  
%     ##     ## #########          ##        ##     ## ##             ## ##             ## 
%     ##     ## ##     ##          ##        ##     ## ##       ##    ## ##       ##    ## 
%     ########  ##     ##          ##         #######  ########  ######  ########  ######  


% pulses

load('PulseData.mat')
% make all the model predictions
qx = NaN(length(PulseData),1);
qy = NaN(length(PulseData),1);
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		PulseData(i).fp(:,j) = DAModelv2(PulseData(i).stim(:,j),p);
	end
	[qx(i), qy(i)] = GeffenMeister(PulseData(i).resp,PulseData(i).fp);
end

% show this on the plot
plot(axes_handles(6),qx,qy,'k+')


% show the timing data
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).resp(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
resp_time = resp_half_time - stim_half_time;
foreground_stim(c:end) = [];

% clean up
rm_this = foreground_stim < 1e-2 | resp_time < 0 | resp_time > 100;
foreground_stim(rm_this) = [];
resp_time(rm_this) = [];
resp_time_data = resp_time;

plot(axes_handles(9),foreground_stim,resp_time,'k+')
set(axes_handles(9),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
xlabel(axes_handles(9),'Mean Stimulus (V)')
ylabel(axes_handles(9),'\tau_{ORN}-\tau_{PID} (ms)','interpreter','tex')

% show the timing data for the DA model
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).fp(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
foreground_stim(c:end) = [];

resp_time = resp_half_time - stim_half_time;

% clean up exactly like we did the data
foreground_stim(rm_this) = [];
resp_time(rm_this) = [];

l=plot(axes_handles(9),foreground_stim,resp_time,'r+');

s=rsquare(resp_time,resp_time_data);
set(axes_handles(4),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
legend(l,strcat('r^2=',oval(s)))


%      ########     ###              ##      ## ######## ########  ######## ########  
%      ##     ##   ## ##             ##  ##  ## ##       ##     ## ##       ##     ## 
%      ##     ##  ##   ##            ##  ##  ## ##       ##     ## ##       ##     ## 
%      ##     ## ##     ##           ##  ##  ## ######   ########  ######   ########  
%      ##     ## #########           ##  ##  ## ##       ##     ## ##       ##   ##   
%      ##     ## ##     ##           ##  ##  ## ##       ##     ## ##       ##    ##  
%      ########  ##     ##            ###  ###  ######## ########  ######## ##     ## 

% show weber scaling 
load('MSG_per_neuron.mat','MSG_data')
load('DA_Fit_to_MSG.mat')

% compute all the model predictions and the fit qualities 
qx = NaN(8,13);
qy = NaN(8,13);
for i = 1:8 % there are 8 paradimgs in the MSG data
	for j = 1:13 % there are 13 neurons
		if ~isempty(p(j).A) && width(MSG_data(i,j).stim) > 1
			for k = 1:width(MSG_data(i,j).stim)
				MSG_data(i,j).fp_DA(:,k) = DAModelv2(MSG_data(i,j).stim(:,k),p(j));
			end
			% compute the geffen-Meister error metric
			[qx(i,j) qy(i,j)] = GeffenMeister(MSG_data(i,j).resp,MSG_data(i,j).fp_DA);
		end
	end
end

% plot on the error plot, and colour code by stimulus
c = parula(9);
for i = 1:8
	this_qx = qx(i,~isnan(qx(i,:)));
	this_qy = qy(i,~isnan(qy(i,:)));
	plot(axes_handles(6),this_qx,this_qy,'x','Color',c(i,:))
end


% get the gain scaling factor -- this is to correct for the fact that we normalised the filter heights everywhere
x = mean2(horzcat(MSG_data(1,:).stim));
y = mean2(horzcat(MSG_data(1,:).resp));
gsf = std(y)/std(x);


% show the gain plot again, and also the gain for the da model 
gain = NaN(8,13);
gain_DA = NaN(8,13);
mean_stim = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			x = MSG_data(i,j).fp;
			x_da = MSG_data(i,j).fp_DA;
			if ~isvector(x)
				x = mean2(x);
			end
			if ~isvector(y)
				y = mean2(y);
			end 
			if ~isvector(x_da)
				x_da = mean2(x_da);
			end 


			% trim NaNs again
			rm_this = isnan(x) | isnan(y) | isnan(x_da);
			x(rm_this) = [];
			x_da(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

			temp=fit(x(:),x_da(:),'poly1');
			gain_DA(i,j) = temp.p1;

			% get the units of gain right
			gain(i,j) = gain(i,j)*gsf;
			gain_DA(i,j) = gain_DA(i,j)*gsf;

		end
	end	
end
plot(axes_handles(7),mean_stim(:),gain(:),'kx')
l=plot(axes_handles(7),mean_stim(:),gain_DA(:),'rx');
xlabel(axes_handles(7),'Mean Stimulus (V)')
ylabel(axes_handles(7),'Neuron Gain (Hz/V)')
set(axes_handles(7),'XScale','log','YScale','log')

% show the spearman rho 
a = gain(:);
b = gain_DA(:);
rm_this = isnan(a) | isnan(b);
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));




%    ########     ###        ######  ########  ######## ######## ########  ##     ## ########  
%    ##     ##   ## ##      ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## 
%    ##     ##  ##   ##     ##       ##     ## ##       ##       ##     ## ##     ## ##     ## 
%    ##     ## ##     ##     ######  ########  ######   ######   ##     ## ##     ## ########  
%    ##     ## #########          ## ##        ##       ##       ##     ## ##     ## ##        
%    ##     ## ##     ##    ##    ## ##        ##       ##       ##     ## ##     ## ##        
%    ########  ##     ##     ######  ##        ######## ######## ########   #######  ##        


% xcorr for mean shifted gaussians
peak_loc_data = NaN(8,13);
peak_loc_DA = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			this_resp = MSG_data(i,j).resp;
			this_stim = MSG_data(i,j).stim;
			this_pred = MSG_data(i,j).fp_DA;
			if ~isvector(this_resp)
				this_resp = mean2(this_resp);
			end
			if ~isvector(this_stim)
				this_stim = mean2(this_stim);
			end
			if ~isvector(this_pred)
				this_pred = mean2(this_pred);
			end

			a = this_resp - mean(this_resp);
			b = this_stim - mean(this_stim);
			a = a/std(a);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			x = x/max(x);
			[~,loc] = max(x);
			peak_loc_data(i,j) = t(loc);

			% throw out the first 5 seconds
			this_pred(1:5e3) = [];
			b(1:5e3) = [];
			a = this_pred - mean(this_pred);
			a = a/std(a);
			x = xcorr(a,b); % positive peak means a lags b
			x = x/max(x);
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			[~,loc] = max(x);
			peak_loc_DA(i,j) = t(loc);
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

		end
	end
end 

plot(axes_handles(8),mean_stim(:),1e3*peak_loc_data(:),'kx')
l=plot(axes_handles(8),mean_stim(:),1e3*peak_loc_DA(:),'rx');
set(axes_handles(8),'XLim',[-.1 4],'YLim',[5 120])
xlabel(axes_handles(8),'Mean Stimulus (V)')
ylabel(axes_handles(8),'Peak xcorr. (ms)')

a = peak_loc_data(:);
b = peak_loc_DA(:);
rm_this = isnan(a) | isnan(b) | a < 0 | b < 0;
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end

end

%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

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



