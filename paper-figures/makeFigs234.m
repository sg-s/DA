% GainChangesWithMeanStimulus.m
% 
% created by Srinivas Gorur-Shandilya at 5:49 , 23 September 2015. Contact me at http://srinivas.gs/contact/
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


clearvars -except being_published fig*


%% Figure 1: ORN gain can be estimated by measuring responses to Gaussian inputs
% Gaussian odorant inputs with short correlation times (A), elicit flickering responses in ORNs that track the odorant stimulus well (B). A linear filter K can be extracted from the odorant input and the firing rate output of the neuron (C). The slope of the residuals in a plot of the firing response vs. the linear prediction (D) is defined as the gain of the ORN. Here, we measure the ORN gain in the linear regime: the linear filter accounts for 96% of the variance in the ORN response (red line), and adding an output nonlinearity (dotted black line), only accounts for an additional 1%. The odorant used is ethyl acetate, stimulating the ab3A neuron. Shading in all plots shows the standard error of the mean. 

figure('outerposition',[0 0 800 700],'PaperUnits','points','PaperSize',[800 700]); hold on
axes_handles(1) = subplot(7,2,1:4); hold on
axes_handles(2) = subplot(7,2,5:8); hold on
axes_handles(3) = subplot(7,2,[9 11 13]); hold on
axes_handles(4) = subplot(7,2,[10 12 14]); hold on
% load cached data
load('../data/MeanShiftedGaussians.mat')

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

% plot some stimulus
example_dose = 3;
plot_these=find(strcmp(paradigm_names{example_dose}, combined_data.paradigm));
plot_this = mean2(combined_data.PID(plot_these,:));
time = dt*(1:length(plot_this));
axes(axes_handles(1))
plot(time,plot_this,'Color',c(example_dose,:));
ylabel(axes_handles(1),'Stimulus (V)')
set(axes_handles(1),'XLim',[45 55])

% load the data cut and processed
load('../data/MSG_per_neuron.mat','MSG_data')
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
K = reshape([MSG_data(example_dose,:).K],1001,length([MSG_data(example_dose,:).K])/1001);
axes(axes_handles(3))
shadedErrorBar(filtertime,mean2(K),sem(K),{'Color',c(example_dose,:)})
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


% plot linear prediction vs. data for the example dose. 
ss = 25;
y = mean2([MSG_data(example_dose,:).resp]);
x = mean2([MSG_data(example_dose,:).fp]);
time = 35+1e-3*(1:length([MSG_data(example_dose,:).fp]));


[ax,plot1,plot2] = plotyy(axes_handles(2),time,y,time,x);
set(ax(1),'XLim',[45 55])
set(ax(2),'XLim',[45 55],'YLim',[-.7 .8 ])
set(plot1,'Color',c(example_dose,:))
set(plot2,'Color','r')
ylabel(ax(1),'ORN Response (Hz)')
ylabel(ax(2),'Projected Stimulus')
set(axes_handles(2),'box','off')

plot(axes_handles(4),x(1:ss:end),y(1:ss:end),'.','Color',c(example_dose,:));

ff= fit(x(~isnan(x)),y(~isnan(x)),'poly1');
l = plot(axes_handles(4),sort(x),ff(sort(x)),'r');
L = {};
L{1} = strcat('Gain=',oval(ff.p1),'Hz/V, r^2=',oval(rsquare(y,x)));

% % also fit a hill function
% rm_this = isnan(x) | isnan(y);
% x(rm_this) = []; y(rm_this) = [];
% minx = min(x);
% x = x - minx;
% clear p
% p.     A= 42.9682;
% p.     k= 1.0570;
% p.     n= 2.4301;
% p.y_offset= 2.9976;
% [x,idx] = sort(x); y = y(idx);
% l(2) = plot(axes_handles(4),x + minx,hill4(x,p),'k--');
% L{2} = strcat('Hill fit, r^2=',oval(rsquare(hill4(x,p),y)));


xlabel(axes_handles(4),'Projected Stimulus')
ylabel(axes_handles(4),'Response (Hz)')

legend(l,L,'Location','northwest');

prettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')

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
for i = 1:8
	plot(ax(1),MSG_data(1,1).time,mean2([MSG_data(i,:).stim]),'Color',c(i,:))
end

% plot the stimulus distributions 
a = floor(35/dt);
z = floor(55/dt);

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	plot_hist = (combined_data.PID(plot_these,a:z));
	[hy,hx]  = hist(plot_hist(:),50);
	hy = hy/length(plot_these);
	hy = hy/sum(hy);
	plot(ax(2),hy,hx,'Color',c(i,:));
end

% plot the responses
for i = 1:8
	temp = [MSG_data(i,:).resp];
	plot(ax(3),MSG_data(1,1).time,mean2(temp),'Color',c(i,:));
	[hx,hy] = hist(mean2(temp),30);
	hx = hx/sum(hx);
	plot(ax(4),hx,hy,'Color',c(i,:));
end


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

	plot(ax(5),x(1:ss:end),y(1:ss:end),'.','Color',c(i,:))
end


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

		end
	end	
end

% calculate the gain error as the standard deviation over each paradigm
gain_err = nanstd(gain');
gain_err(gain_err == 0) = Inf;

% show gain changes -- gain vs. mean stimulus
for i = 1:8
	plot(ax(6),(mean_stim(i,:)),(gain(i,:)),'+','Color',c(i,:));
end


cf = fit(nanmean(mean_stim')',nanmean(gain')','power1','Weights',1./gain_err(:));


l(1)=plot(ax(6),nanmean(mean_stim')',cf(nanmean(mean_stim')'),'k');
r2 = rsquare(nonnans(gain),cf(nonnans(mean_stim)));
L = strcat('y = \alpha x^{\beta}, \beta= ',oval(cf.b), ',r^2 = ',oval(r2));

% fit a power law with exponent -1
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
options.Weights = 1./gain_err(:);
cf = fit(nanmean(mean_stim')',nanmean(gain')','power1',options);
l(2)=plot(ax(6),nanmean(mean_stim')',cf(nanmean(mean_stim')'),'k--');
r2 = rsquare(nonnans(gain),cf(nonnans(mean_stim)));
legend(l,{L, strcat('y = \alpha x^{-1}, \alpha= ',oval(cf.a), ', r^2 = ',oval(r2))},'Location','northoutside')

% rescale by Weber law
ss = 50;
allx = [];
ally = [];
for i = 1:8 % iterate over all paradigms 
	y = ([MSG_data(i,:).resp]);
	x = ([MSG_data(i,:).fp]);
	s = ([MSG_data(i,:).stim]);
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(s)
		s = mean2(s);
	end
	if ~isvector(y)
		y = mean2(y);
	end 

	allx = [allx mean(y)+cf(mean2(s))*(x(1:ss:end))];
	ally = [ally y(1:ss:end)];
	plot(ax(7),mean(y)+cf(mean2(s))*(x(1:ss:end)),y(1:ss:end),'.','Color',c(i,:))
end
allx = allx(:); ally = ally(:); 
rm_this = isnan(allx) | isnan(ally);
allx(rm_this) = [];
ally(rm_this) = [];
ff2 = fit(allx(:),ally(:),'poly1');
clear l
l = plot(ax(7),-5:0.1:45,ff2(-5:0.1:45),'r');
legend(l,strcat('r^2=',oval(rsquare(ally,ff2(allx)))),'Location','northwest');


% cosmetics
set(ax(1),'XLim',[45 55],'YLim',[0 4.5])
set(ax(2),'YLim',[0 4.5],'YTick',[])

set(ax(3),'XLim',[45 55],'YLim',[0 40])
set(ax(4),'YLim',[0 40],'YTick',[])
set(ax(6),'XScale','log','YScale','log','YLim',[1 50],'XLim',[.5 3.5])
set(ax(7),'XLim',[0 50],'YLim',[0 50])

ylabel(ax(1),'Stimulus (V)')
ylabel(ax(3),'ORN Response (Hz)')
xlabel(ax(3),'Time (s)')
xlabel(ax(4),'Probability')
xlabel(ax(5),'Projected Stimulus')
ylabel(ax(5),'ORN Response (Hz)')
xlabel(ax(6),'Mean Stimulus (V)')
ylabel(ax(6),'Neuron Gain (Hz/V)')

xlabel(ax(7),[' Projected Stimulus ' char(10) 'Rescaled by Weber Law'])
ylabel(ax(7),'Neuron Response (Hz)')



prettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')

if being_published
	snapnow
	delete(gcf)
end


%         ######  ##     ## ########  ########     ######## ####  ######         ##   
%        ##    ## ##     ## ##     ## ##     ##    ##        ##  ##    ##      ####   
%        ##       ##     ## ##     ## ##     ##    ##        ##  ##              ##   
%         ######  ##     ## ########  ########     ######    ##  ##   ####       ##   
%              ## ##     ## ##        ##           ##        ##  ##    ##        ##   
%        ##    ## ##     ## ##        ##           ##        ##  ##    ##        ##   
%         ######   #######  ##        ##           ##       ####  ######       ###### 

%% Supplementary Figure 1

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
ss = 50;
all_x = []; all_y = [];
for i = 1:8 % iterate over all paradigms 
	y = ([MSG_data(i,:).resp]);
	s = ([MSG_data(i,:).stim]);
	x = ([MSG_data(i,:).fp]);
	if ~isvector(x)
		x = mean2(x);
	end

	if ~isvector(s)
		s = mean2(s);
	end

	if ~isvector(y)
		y = mean2(y);
	end 

	all_x = [all_x x+mean(s)];
	all_y = [all_y y];
	plot(x(1:ss:end)+mean(s),y(1:ss:end),'.','Color',c(i,:))
end

p.A= 34.4229;
p.k= 0.7455;
p.n= 1.5538;
L = ['Hill fit, r^2=' oval(rsquare(all_y(:),hill(all_x(:),p)))];
all_x = nonnans(sort(all_x(:)));
all_x = linspace(all_x(1),all_x(end),100);
l = plot(all_x,hill(all_x,p),'r');
legend(l,L,'Location','southeast')
xlabel(['Projected Stimulus +' char(10) 'mean stimulus'])
ylabel('Neuron Response (Hz)')
set(gca,'XLim',[0 5],'YLim',[0 45])


% show non-linearities to each cloud
clear p
subplot(1,2,2), hold on
for i = 1:8 % iterate over all paradigms 
	y = ([MSG_data(i,:).resp]);
	s = ([MSG_data(i,:).stim]);
	x = ([MSG_data(i,:).fp]);
	if ~isvector(x)
		x = mean2(x);
	end

	if ~isvector(s)
		s = mean2(s);
	end

	if ~isvector(y)
		y = mean2(y);
	end 

	clear d
	d.stimulus = x + mean(s);
	d.response = y;
	p(i) = fitModel2Data(@hill,d,'nsteps',0);
	plot(all_x(all_x<max(d.stimulus)),hill(all_x(all_x<max(d.stimulus)),p(i)),'Color',c(i,:))
end
set(gca,'XLim',[0 5],'YLim',[0 45])
xlabel(['Projected Stimulus +' char(10) 'mean stimulus'])
ylabel('Neuron Response (Hz)')

prettyFig('fs=20;');

if being_published
	snapnow
	delete(gcf)
end


%        ######  ########  ######## ######## ########  ##     ## ########   ######  
%       ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## ##    ## 
%       ##       ##     ## ##       ##       ##     ## ##     ## ##     ## ##       
%        ######  ########  ######   ######   ##     ## ##     ## ########   ######  
%             ## ##        ##       ##       ##     ## ##     ## ##              ## 
%       ##    ## ##        ##       ##       ##     ## ##     ## ##        ##    ## 
%        ######  ##        ######## ######## ########   #######  ##         ######  

%% Figure 3: ORNs speed up responses on increasing stimulus mean
% ORN gain may permit responses to occur with smaller delays, as ORNs need to integrate the stimulus over a smaller duration to respond. Speedup in responses can be estimated using the peak times of linear filters fit to increasing mean concentrations of odorant (colours, in A). ORN response speedups with increasing stimulus mean can also be estimated in a model-free manner using the spike-triggered average (STA, B). Both methods suggest that response delays decrease with increasing mean stimulus (C). 

figure('outerposition',[0 0 1100 800],'PaperUnits','points','PaperSize',[1100 800]); hold on
subplot(2,2,1), hold on

ac_mean = zeros(20e3+1,8);
ac_std = zeros(20e3+1,8);
a = [];
for i = 1:8
	this_ac = [];
	for j = 1:13
		
		stim = (MSG_data(i,j).stim);
		if ~isempty(stim)
			if ~isvector(stim)
				stim = mean2(stim);
			end
			[~,~,temp]=findCorrelationTime(stim);
			this_ac = [this_ac temp];
		end
	end
	if isvector(this_ac)
		ac_mean(:,i) = (this_ac);
	else
		ac_mean(:,i) = mean2(this_ac);
	end
	ac_std(:,i) = sem(this_ac);
end

time = 1e-3*(1:length(ac_mean));
for i = 1:8
	[~,si(i)] = errorShade(time,ac_mean(:,i),ac_std(:,i),'Color',c(i,:));
end
for i = 1:8
	uistack(si(i), 'bottom')
end
set(gca,'XScale','log','XLim',[5e-3 2],'YLim',[-.4 1])
xlabel('Lag (s)')
ylabel('Stimulus Autocorrelation')

ac_mean = zeros(20e3+1,8);
ac_std = zeros(20e3+1,8);
a = [];
for i = 1:8
	this_ac = [];
	for j = 1:13
		
		stim = (MSG_data(i,j).resp);
		if ~isempty(stim)
			if ~isvector(stim)
				stim = mean2(stim);
			end
			[~,~,temp]=findCorrelationTime(stim);
			this_ac = [this_ac temp];
		end
	end
	if isvector(this_ac)
		ac_mean(:,i) = (this_ac);
	else
		ac_mean(:,i) = mean2(this_ac);
	end
	ac_std(:,i) = sem(this_ac);
end

subplot(2,2,2), hold on
time = 1e-3*(1:length(ac_mean));
for i = 1:8
	[~,si(i)] = errorShade(time,ac_mean(:,i),ac_std(:,i),'Color',c(i,:));
end
for i = 1:8
	uistack(si(i), 'bottom')
end
set(gca,'XScale','log','XLim',[5e-3 2],'YLim',[-.4 1])
xlabel('Lag (s)')
ylabel('Response Autocorrelation')


% fit parametric filters to the raw, neuron-wise filters extracted earlier 
% for i = 1:8
% 	for j = 1:13
% 		if ~isempty(MSG_data(i,j).K)
% 			d.stimulus = MSG_data(i,j).K(200:end);
% 			d.response = MSG_data(i,j).K(200:end);
% 			for k = 1:5
% 				MSG_data(i,j).p = fitModel2Data(@FitFilter,d,MSG_data(i,j).p);
% 			end
% 		end
% 	end
% end

% show filter speedups
subplot(2,3,4), hold on

% compute peak locations of all these filters
clear l 
l = zeros(8,1);
peak_loc_K = NaN(8,13);
for i = 1:8
	for j = 1:13
		if ~isempty(MSG_data(i,j).K)
			K2 = pFilter(MSG_data(i,j).K(200:end),MSG_data(i,j).p);
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
legend(l,L);

before = 1e4;
after = 5e3;

% compute STA for all the data
% allfiles = dir('/local-data/DA-paper/fast-flicker/orn/*.mat');
% allSTA = [];
% mean_stim = []; 
% dil = [];
% for i= 1:length(allfiles)
% 	load(['/local-data/DA-paper/fast-flicker/orn/' allfiles(i).name])
% 	for j = 1:length(spi`kes)
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
% 						`
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
% save this
% save('allSTA.mat','allSTA','dil','mean_stim');


% plot the STA
subplot(2,3,5), hold on

old_mean_stim = mean_stim;
load('../data/allSTA.mat','allSTA','dil','mean_stim');

udil = unique(dil);
for i = 1:8
	temp = (allSTA(:,dil == udil(i)));
	for j = 1:width(temp)
		temp(:,j) = temp(:,j)/max(temp(:,j));
	end
	t = -after:before;
	t = t*1e-3; 
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

subplot(2,3,6), hold on
clear l
l(1) = errorbar(x,y,ey,'k');
peak_loc_K = peak_loc_K(~isnan(peak_loc_K));
l(2) = plot(nonnans(old_mean_stim(:)),peak_loc_K(:)/dt,'r+');

% calculate Spearman's rho (http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
s2 = spear(nonnans(old_mean_stim),peak_loc_K);
s1 = spear(mean_stim,peak_STA);

legend(l,{strcat('STA, \rho=',oval(s1)), strcat('Filter, \rho=',oval(s2))})
set(gca,'YLim',[-10 130])
ylabel('Peak time (ms)')
xlabel('Mean Stimulus (V)')

prettyFig('fs=20;');

if being_published
	snapnow
	delete(gcf)
end

prettyFig()

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