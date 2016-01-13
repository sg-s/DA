% GainChangesWithMeanStimulus.m
% 
% created by Srinivas Gorur-Shandilya at 5:49 , 23 September 2015. Contact me at http://srinivas.gs/contact/
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

% load data 

[PID, ~, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

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
ylabel(ax(2),'Projected Stimulus')
set(axes_handles(2),'box','off')
set(ax(1),'XMinorTick','on','YMinorTick','on')
set(ax(2),'XMinorTick','on','YMinorTick','on')
xlabel('Time (s)')

plot(axes_handles(4),x(1:25:end),y(1:25:end),'.','Color',c(example_dose,:));
ff= fit(x(~isnan(x)),y(~isnan(x)),'poly1');
l = plot(axes_handles(4),sort(x),ff(sort(x)),'k');
L = {};
L{1} = strcat('Gain=',oval(ff.p1),'Hz/V, r^2=',oval(rsquare(y,x)));


xlabel(axes_handles(4),'Projected Stimulus (V)')
ylabel(axes_handles(4),'ORN Response (Hz)')
set(axes_handles(4),'YColor',c(example_dose,:),'XColor','r')

legend(l,L,'Location','northwest');

prettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')



movePlot(axes_handles(3),'down',.03)
movePlot(axes_handles(4),'down',.03)
movePlot(axes_handles(1),'up',.06)
movePlot(axes_handles(2),'up',.03)

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



prettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=0;','FixLogY=0;')

if being_published
	snapnow
	delete(gcf)
end



% test figure: optimal coding in the mean shifted gaussian case
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(3,3,1), hold on
x = 0:0.05:2;
y = zeros(length(x),max(paradigm));
for i = 1:max(paradigm)
	stim = PID(a:z,paradigm == i);
	[hy,hx]  = hist(stim(:),50);
	temp = histcounts(stim(:),x);
	temp = temp/sum(temp);
	temp = cumsum(temp);
	y(2:end,i) = temp;
	hy = hy/sum(hy);
	plot(hx,hy,'Color',c(i,:));
end
ylabel('p(stimulus)')
set(gca,'XLim',[0 2])

% normalise ORN data
m = max(max([orn_io_data.y]));
for i = 1:max(paradigm)
	orn_io_data(i).y = orn_io_data(i).y/m;
end

% now plot the integrals -- these are the theoretical predictions
% also plot the data
subplot(3,3,4), hold on
for i = 1:max(paradigm)
	plot(x,y(:,i),'Color',c(i,:));
	plot(orn_io_data(i).x,orn_io_data(i).y,'+','Color',c(i,:))
end
xlabel('Stimulus (V)')
ylabel('Response (norm)')
set(gca,'XLim',[0 2])

% now show the correlation between predicted gain and actual gain
subplot(3,3,7), hold on
orn_gain = NaN(max(paradigm),1);
predicted_gain = NaN*orn_gain;
for i = 1:max(paradigm)
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',temp(~isnan(temp))','poly1');
	orn_gain(i) = temp.p1;
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',y(~isnan(temp),i),'poly1');
	predicted_gain(i) = temp.p1;
end
plot(predicted_gain,orn_gain,'k+')
plot([1e-3 10],[1e-3 10],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.1 3])
xlabel('Predicted Gain')
ylabel('ORN Gain')


% now plot the stimulus + best fit assuming Kd constraint
clear p
p = cache('Kd_constrained_MSG_fit2');
if isempty(p)
	clear d
	d.stimulus = x;
	d.response = y(:,1)';
	p = fitModel2Data(@hill,d,'make_plot',false,'nsteps',1000,'display_type','iter');
	ub = p;
	ub.k = Inf;
	lb = p;
	lb.k = 0;
	for i = 2:max(paradigm)
		clear d
		d.stimulus = x;
		d.response = y(:,i)';
		p(i) = fitModel2Data(@hill,d,'p0',p(1),'ub',ub,'lb',lb,'make_plot',false,'nsteps',300,'display_type','iter');
	end
	cache('Kd_constrained_MSG_fit2',p);
end

subplot(3,3,2), hold on

for i = 1:max(paradigm)
	plot(x,y(:,i),'Color',c(i,:));
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
end
title('Only K_D can vary')
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

subplot(3,3,5), hold on
for i = 1:max(paradigm)
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
	plot(orn_io_data(i).x,orn_io_data(i).y,'+','Color',c(i,:))
end
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

% now show the correlation between predicted gain and actual gain
subplot(3,3,8), hold on
orn_gain = NaN(max(paradigm),1);
predicted_gain = NaN*orn_gain;
for i = 1:max(paradigm)
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',temp(~isnan(temp))','poly1');
	orn_gain(i) = temp.p1;
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	pred = hill(x,p(i));
	temp = fit(x(~isnan(temp))',pred(~isnan(temp))','poly1');
	predicted_gain(i) = temp.p1;
end
plot(predicted_gain,orn_gain,'k+')
plot([1e-3 10],[1e-3 10],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.1 3])
xlabel('Predicted Gain')

% now plot the stimulus + best fit assuming n constraint
clear p
p = cache('n_constrained_MSG_fit2');
if isempty(p)
	clear d
	d.stimulus = x;
	d.response = y(:,1)';
	p = fitModel2Data(@hill,d,'make_plot',false,'nsteps',1000,'display_type','iter');
	ub = p;
	ub.n = Inf;
	lb = p;
	lb.n = 0;
	for i = 2:max(paradigm)
		clear d
		d.stimulus = x;
		d.response = y(:,i)';
		p(i) = fitModel2Data(@hill,d,'p0',p(1),'ub',ub,'lb',lb,'make_plot',false,'nsteps',300,'display_type','iter');
	end
	cache('n_constrained_MSG_fit2',p);
end

subplot(3,3,3), hold on
for i = 1:max(paradigm)
	plot(x,y(:,i),'Color',c(i,:));
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
end
title('Only n can vary')
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

subplot(3,3,6), hold on
for i = 1:max(paradigm)
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
	plot(orn_io_data(i).x,orn_io_data(i).y,'+','Color',c(i,:))
end
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

% now show the correlation between predicted gain and actual gain
subplot(3,3,9), hold on
orn_gain = NaN(max(paradigm),1);
predicted_gain = NaN*orn_gain;
for i = 1:max(paradigm)
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',temp(~isnan(temp))','poly1');
	orn_gain(i) = temp.p1;
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	pred = hill(x,p(i));
	temp = fit(x(~isnan(temp))',pred(~isnan(temp))','poly1');
	predicted_gain(i) = temp.p1;
end
plot(predicted_gain,orn_gain,'k+')
plot([1e-3 10],[1e-3 10],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.1 3],'XTick',[1e-10 1e-6 1e-3 1e1])
xlabel('Predicted Gain')

prettyFig('plw=1.3;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end


% another supplementary figure showing that variance changes are not important here. 
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(nanmean(PID(a:z,:)),nanstd(PID(a:z,:)),'k+')
plot([0 2],[0 2],'k--')
xlabel('\mu_{stimulus} (V)')
ylabel('\sigma_{stimulus (V)}')

prettyFig('plw=1.3;','lw=1.5;','fs=20;')

if being_published
	snapnow
	delete(gcf)
end

% subplot(1,3,2), hold on
% for i = 1:max(paradigm)
% 	plot(mean_stim(paradigm==i),fA_gain(paradigm==i),'+','Color',c(i,:));
% end
% l = plot(NaN,NaN,'k+');
% legend(l,['\rho = ' oval(spear(mean_stim,fA_gain))])
% set(gca,'XScale','log','YScale','log')
% xlabel('Mean Stimulus (V)')
% ylabel('ORN Gain (Hz/V)')

% subplot(1,3,3), hold on
% std_stim = nanstd(PID(a:z,:));
% for i = 1:max(paradigm)
% 	plot(std_stim(paradigm==i),fA_gain(paradigm==i),'+','Color',c(i,:));
% end
% l = plot(NaN,NaN,'k+');
% legend(l,['\rho = ' oval(spear(std_stim(:),fA_gain))])
% set(gca,'XScale','log','YScale','log','XLim',[min(std_stim) max(std_stim)])
% xlabel('Std Stimulus (V)')
% ylabel('ORN Gain (Hz/V)')


% % OK, now we show that both the mean and the variance can account for gain changes. now we show that in this stimulus, the mean and the variance co-vary

% all_block_sizes = factor2(20e3);
% all_block_sizes = all_block_sizes(6:24);
% clear l r2
% r2 = NaN(length(all_block_sizes),8);
% for j = 1:8
% 	for i = 1:length(all_block_sizes)
% 		temp = [MSG_data(j,:).stim];
% 		temp = temp(2:end,:);
% 		temp = temp(:);
% 		temp = reshape(temp,all_block_sizes(i),length(temp)/all_block_sizes(i));
% 		r2(i,j) = rsquare(mean(temp),std(temp));
% 	end
% 	plot(all_block_sizes,r2(:,j),'Color',c(j,:))
% end
% set(gca,'XScale','log','YLim',[0 1])
% ylabel('r^2(\sigma,\mu)')
% xlabel('Window Length (ms)')


% % compute gain changes on a per-neuron basis
% gain = NaN(8,13);
% mean_stim = NaN(8,13);
% std_stim = NaN(8,13);
% for i = 1:8 % iterate over all paradigms 
% 	for j = 1:13
% 		if width(MSG_data(i,j).stim) > 1
% 			y = MSG_data(i,j).resp; % average over all neurons 
% 			x = MSG_data(i,j).fp;
% 			if ~isvector(x)
% 				x = mean(x,2);
% 			end
% 			if ~isvector(y)
% 				y = mean(y,2);
% 			end 
			
% 			% trim NaNs again
% 			rm_this = isnan(x) | isnan(y);
% 			x(rm_this) = [];
% 			y(rm_this) = [];

% 			temp=fit(x(:),y(:),'poly1');
% 			gain(i,j) = temp.p1;
% 			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));
% 			std_stim(i,j) = mean(std([MSG_data(i,j).stim]));
% 		end
% 	end	
% end

% subplot(1,3,2), hold on
% plot(mean_stim(:),gain(:),'k+')

% subplot(1,3,3), hold on
% plot(std_stim(:),gain(:),'k+')


%        ######  ########  ######## ######## ########  ##     ## ########   ######  
%       ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## ##    ## 
%       ##       ##     ## ##       ##       ##     ## ##     ## ##     ## ##       
%        ######  ########  ######   ######   ##     ## ##     ## ########   ######  
%             ## ##        ##       ##       ##     ## ##     ## ##              ## 
%       ##    ## ##        ##       ##       ##     ## ##     ## ##        ##    ## 
%        ######  ##        ######## ######## ########   #######  ##         ######  

% no longer doing this figure

% Figure 3: ORNs speed up responses on increasing stimulus mean
% ORN gain may permit responses to occur with smaller delays, as ORNs need to integrate the stimulus over a smaller duration to respond. Speedup in responses can be estimated using the peak times of linear filters fit to increasing mean concentrations of odorant (colours, in A). ORN response speedups with increasing stimulus mean can also be estimated in a model-free manner using the spike-triggered average (STA, B). Both methods suggest that response delays decrease with increasing mean stimulus (C). 

% figure('outerposition',[0 0 1100 800],'PaperUnits','points','PaperSize',[1100 800]); hold on
% clear axes_handles
% axes_handles(1) = subplot(2,2,1); hold on

% ac_mean = zeros(20e3+1,8);
% ac_std = zeros(20e3+1,8);
% clear a
% for i = 1:8
% 	this_ac = [];
% 	for j = 1:13
		
% 		stim = (MSG_data(i,j).stim);
% 		if ~isempty(stim)
% 			if ~isvector(stim)
% 				stim = mean(stim,2);
% 			end
% 			[~,~,temp]=findCorrelationTime(stim);
% 			this_ac = [this_ac temp];
% 		end
% 	end
% 	if isvector(this_ac)
% 		ac_mean(:,i) = (this_ac);
% 	else
% 		ac_mean(:,i) = mean(this_ac,2);
% 	end
% 	ac_std(:,i) = sem(this_ac);
% end

% time = 1e-3*(1:length(ac_mean));
% for i = 1:8
% 	[~,si(i)] = errorShade(time,ac_mean(:,i),ac_std(:,i),'Color',c(i,:));
% end
% for i = 1:8
% 	uistack(si(i), 'bottom')
% end
% set(gca,'XScale','log','XLim',[5e-3 2],'YLim',[-.4 1])
% xlabel('Lag (s)')
% ylabel('Stimulus Autocorrelation')

% ac_mean = zeros(20e3+1,8);
% ac_std = zeros(20e3+1,8);
% a = [];
% for i = 1:8
% 	this_ac = [];
% 	for j = 1:13
		
% 		stim = (MSG_data(i,j).resp);
% 		if ~isempty(stim)
% 			if ~isvector(stim)
% 				stim = mean(stim,2);
% 			end
% 			[~,~,temp]=findCorrelationTime(stim);
% 			this_ac = [this_ac temp];
% 		end
% 	end
% 	if isvector(this_ac)
% 		ac_mean(:,i) = (this_ac);
% 	else
% 		ac_mean(:,i) = mean(this_ac,2);
% 	end
% 	ac_std(:,i) = sem(this_ac);
% end

% axes_handles(2) = subplot(2,2,2); hold on
% time = 1e-3*(1:length(ac_mean));
% for i = 1:8
% 	[~,si(i)] = errorShade(time,ac_mean(:,i),ac_std(:,i),'Color',c(i,:));
% end
% for i = 1:8
% 	uistack(si(i), 'bottom')
% end
% set(gca,'XScale','log','XLim',[5e-3 2],'YLim',[-.4 1])
% xlabel('Lag (s)')
% ylabel('Response Autocorrelation')


% % fit parametric filters to the raw, neuron-wise filters extracted earlier 
% % for i = 1:8
% % 	for j = 1:13
% % 		if ~isempty(MSG_data(i,j).K)
% % 			d.stimulus = MSG_data(i,j).K(200:end);
% % 			d.response = MSG_data(i,j).K(200:end);
% % 			for k = 1:5
% % 				MSG_data(i,j).p = fitModel2Data(@FitFilter,d,MSG_data(i,j).p);
% % 			end
% % 		end
% % 	end
% % end

% % show filter speedups
% axes_handles(3) = subplot(2,2,3); hold on

% % compute peak locations of all these filters
% clear l 
% l = zeros(8,1);
% peak_loc_K = NaN(8,13);
% mean_stim_K = NaN(8,13);
% for i = 1:8
% 	for j = 1:13
% 		if ~isempty(MSG_data(i,j).K)
% 			K2 = pFilter(MSG_data(i,j).K(200:end),MSG_data(i,j).p);
% 			filtertime = 1e-3*(1:length(K2));
% 			l(i)=plot(filtertime,K2,'Color',c(i,:));
% 			[~,loc] = max(K2);
% 			peak_loc_K(i,j) = filtertime(loc);

% 			mean_stim_K(i,j) = mean(mean(MSG_data(i,j).stim));

% 		end
% 	end
% end


% set(gca,'XLim',[-.01 .5])
% xlabel('Lag (s)')
% ylabel('Filter (norm)')
% L = paradigm_names;
% for i = 1:length(L)
% 	L{i} = L{i}(strfind(L{i},'-')+1:end);
% end
% legend(l,L);


% % compute STA for all the data
% before = 1e4;
% after = 5e3;

% allfiles = dir('/local-data/DA-paper/fast-flicker/orn/*.mat');
% allSTA = [];
% mean_stim = []; 
% dil = [];
% for i= 1:length(allfiles)
% 	load(['/local-data/DA-paper/fast-flicker/orn/' allfiles(i).name])
% 	disp(['/local-data/DA-paper/fast-flicker/orn/' allfiles(i).name])
% 	for j = 1:length(spikes)
% 		textbar(j,length(spikes))
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

% 			mean_stim = [mean_stim; mean(this_stim,2)];
% 			this_dil = str2double(ControlParadigm(j).Name(strfind(ControlParadigm(j).Name,'-')+1:strfind(ControlParadigm(j).Name,'%')-1));
% 			dil = [dil;this_dil*ones(width(this_stim),1)  ];

% 			this_STA = STA(these_spikes,this_stim,'normalise',true,'regulariseParameter',1,'before',before,'after',after);
% 			allSTA = [allSTA this_STA];
			
% 		end
% 	end
% end
% % normalise all of them
% for i = 1:width(allSTA)
% 	allSTA(:,i) = allSTA(:,i)/max(allSTA(:,i));
% end

% % save this
% save('.cache/allSTA.mat','allSTA','dil','mean_stim');


% % plot the STA
% subplot(2,3,5), hold on

% old_mean_stim = mean_stim;
% load('.cache/allSTA.mat','allSTA','dil','mean_stim');

% udil = unique(dil);
% for i = 1:8
% 	temp = (allSTA(:,dil == udil(i)));
% 	for j = 1:width(temp)
% 		temp(:,j) = temp(:,j)/max(temp(:,j));
% 	end
% 	t = 1e-3*(1:length(temp)) - .44;
% 	errorShade(t,flipud(mean2(temp)),flipud(sem(temp)),'Color',c(i,:));
% end
% set(gca,'XLim',[-.2 .5])
% xlabel('Lag (s)')
% ylabel('STA (norm)')

% % find the peak of each STA
% peak_STA = NaN*mean_stim;
% for i = 1:length(mean_stim)
% 	[~,loc] = max(allSTA(:,i));
% 	peak_STA(i) = (1e3 - loc) + 60;
% end

% x = NaN(8,1); y = x; ex = x; ey = x;
% peak_STA(peak_STA>900) = NaN; % throw out some bullshit values
% for i = 1:8
% 	x(i) = mean2(mean_stim(dil == udil(i)));
% 	ex(i) = sem(mean_stim(dil == udil(i)));
% 	y(i) = nanmean(peak_STA(dil == udil(i)));
% 	ey(i) = sem(nonnans(peak_STA(dil == udil(i))));
% end

% axes_handles(4) = subplot(2,2,4); hold on


% l = errorbar(nanmean(mean_stim'),1e3*nanmean(peak_loc_K'),1e3*nanstd(peak_loc_K')./sqrt(sum(~isnan(peak_loc_K)')),'k.');


% % calculate Spearman's rho (http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
% s2 = spear(nonnans(mean_stim_K(:)),nonnans(peak_loc_K(:)));

% legend(l,{strcat('\rho=',oval(s2))})
% set(gca,'YLim',[20 130])
% ylabel('Peak time (ms)')
% xlabel('Mean Stimulus (V)')

% prettyFig('fs=20;');

% set(axes_handles(1),'XTick',[1e-2 1e-1 1e0 ])
% set(axes_handles(2),'XTick',[1e-2 1e-1 1e0 ])


% if being_published
% 	snapnow
% 	delete(gcf)
% end

% prettyFig()

% if being_published
% 	snapnow
% 	delete(gcf)
% end


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
% This file has the following external dependencies:
showDependencyHash(mfilename);

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
