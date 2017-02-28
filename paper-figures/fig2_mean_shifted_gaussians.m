% GainChangesWithMeanStimulus.m
% 
% created by Srinivas Gorur-Shandilya at 5:49 , 23 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

% use dataManager 
dm = dataManager;

% assemble data

clear cdata
cdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
A_spikes = cdata.A_spikes;
B_spikes = cdata.B_spikes;
cdata = cleanMSGdata(cdata);

K = (mean(cdata.K2(:,cdata.paradigm==1),2));

v2struct(cdata)

clearvars LFP_pred LFP_gain fB cdata

% remove some bad trials
bad_trials =  max(fA) == 0;
fA(:,bad_trials) = [];
fA_pred(:,bad_trials) = [];
PID(:,bad_trials) = [];
K2(:,bad_trials) = [];
orn(bad_trials) = [];
fly(bad_trials) = [];
fA_gain(bad_trials) = [];
paradigm(bad_trials) = [];
LFP(:,bad_trials) = [];
A_spikes(bad_trials,:) = [];


% some core variables
dt = 1e-3;
c = parula(max(paradigm)+1); % colour scheme




;;     ;;    ;;;    ;;;; ;;    ;;    ;;;;;;;; ;;;;  ;;;;;;   
;;;   ;;;   ;; ;;    ;;  ;;;   ;;    ;;        ;;  ;;    ;;  
;;;; ;;;;  ;;   ;;   ;;  ;;;;  ;;    ;;        ;;  ;;        
;; ;;; ;; ;;     ;;  ;;  ;; ;; ;;    ;;;;;;    ;;  ;;   ;;;; 
;;     ;; ;;;;;;;;;  ;;  ;;  ;;;;    ;;        ;;  ;;    ;;  
;;     ;; ;;     ;;  ;;  ;;   ;;;    ;;        ;;  ;;    ;;  
;;     ;; ;;     ;; ;;;; ;;    ;;    ;;       ;;;;  ;;;;;;   


%% Figure 2: ORN gain decreases with increasing stimulus intensity, similar to the Weber-Fechner Law
%

figure('PaperUnits','centimeters','PaperSize',[20 20],'Position',[100 100 800 800]); hold on
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
time = 1e-3*(1:length(PID));
for i = 1:max(paradigm)
	plot(ax(1),time(a:z),nanmean(PID(a:z,paradigm==i),2),'Color',c(i,:))
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
	temp(:,sum(temp) == 0) = [];
	plot(ax(3),time(a:z),nanmean(temp,2),'Color',c(i,:));
	[hx,hy] = hist(mean(temp,2),30);
	hx = hx/sum(hx);
	plot(ax(4),hx,hy,'Color',c(i,:));
end

% show gain changes for all paradigms -- average over neurons 
ss = 100;
all_x = 0:0.1:2;
axes(ax(5)), hold(ax(5),'on')


clear th
for i = 1:max(paradigm) % iterate over all paradigms 
	y = fA(a:z,paradigm == i);
	x = fA_pred(a:z,paradigm == i);
	s = PID(a:z,paradigm == i);
	rm_this = sum(y)==0;
	y(:,rm_this) = [];
	x(:,rm_this) = [];
	s(:,rm_this) = [];
	y = nanmean(y,2);
	x = nanmean(x,2);
	s = nanmean(s,2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	
	[~,orn_io_data(i)] = plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:),'show_error',false,'LineWidth',3);
	t = oval(rsquare(x,y));
	t = t(2:end);
	th(i) = text(.05+max(orn_io_data(i).x),5+max(orn_io_data(i).y),t,'Color',c(i,:));

end

% fix a few positions
P = [
    0.3238   61.0108         0    0.5068   61.6890         0    0.5693  56.8031         0    0.5647   51.2603         0    0.6709   46.8964         0    0.6753   42.5207         0    0.8600   42.1649         0 1.0810   40.1730         0    1.3000   34.1413         0    1.6061 35.0734         0];
P = reshape(P,3,10);
for i = 1:10
	th(i).Position = P(:,i);
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
plot(ax(6),sort(mean_stim),cf(sort(mean_stim)),'r');
set(ax(6),'XScale','log','YScale','log','YLim',[10 300],'XLim',[.1 2.5])

% rescale by Weber law
ss = 50;

axes(ax(7)), hold(ax(7),'on')
for i = 1:8 % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);

	% x = x - nanmean(x);
	% x = x + nanmean(nanmean(s));

	x = cf(nanmean(s))*(x);
	x = x - nanmean(x);
	x = x + nanmean(y);
	[~,temp]= plotPieceWiseLinear(x(1e3:end),y(1e3:end),'nbins',40,'Color',c(i,:),'show_error',false,'LineWidth',3,'make_plot',false);
	plot(ax(7),temp.x,temp.y,'Color',c(i,:),'LineWidth',3)
end

plot(ax(7),[0 80],[0 80],'r')


% cosmetics
set(ax(1),'XLim',[35 55],'YLim',[0 2])
set(ax(2),'YLim',[0 2])

set(ax(3),'XLim',[35 55],'YLim',[0 70])
set(ax(4),'YLim',[0 70])

set(ax(7),'XLim',[0 70],'YLim',[0 70])

ylabel(ax(1),'Stimulus (V)')
ylabel(ax(3),'ab3A firing rate (Hz)')
xlabel(ax(3),'Time (s)')
xlabel(ax(4),'Probability')
xlabel(ax(5),'Projected stimulus (V)')
ylabel(ax(5),'ab3A firing rate (Hz)')
xlabel(ax(6),'Mean stimulus (V)')
ylabel(ax(6),'ab3A ORN gain (Hz/V)')

xlabel(ax(7),[' Projected stimulus ' char(10) 'rescaled by Weber''s law'])
ylabel(ax(7),'ab3A firing rate (Hz)')

prettyFig('plw',1.3,'lw',1.5,'fs',.5,'font_units','centimeters')

axes(ax(6))
text(.2, 30,'$\sim 1/s$','interpreter','latex','Color','r','FontSize',20)

% move some axes to the right
ax(2).Position(1) = .8;
ax(4).Position(1) = .8;

% shrink data for smaller file sizes
shrinkDataInPlot(ax(1),5)
shrinkDataInPlot(ax(3),5)

% move top two rows up a bit
ax(1).Position(2) = .75;
ax(2).Position(2) = .75;
ax(3).Position(2) = .44;
ax(4).Position(2) = .44;

labelFigure;

% deinteresect some axes
ax(1).XLim(2) = 55.05;
deintersectAxes(ax(1))
deintersectAxes(ax(2))
ax(3).XLim(2) = 55.05;
deintersectAxes(ax(3))
deintersectAxes(ax(4))
deintersectAxes(ax(6))

% fix some origins
ax(5).YLim(1) = 0;
ax(5).XLim(1) = 0;
ax(5).Position(1) = .1;

if being_published
	snapnow
	delete(gcf)
end


 ;;;;;;  ;;     ;; ;;;;;;;;  ;;;;;;;;     ;;;;;;;; ;;;;  ;;;;;;   
;;    ;; ;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;    ;;  
;;       ;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;        
 ;;;;;;  ;;     ;; ;;;;;;;;  ;;;;;;;;     ;;;;;;    ;;  ;;   ;;;; 
      ;; ;;     ;; ;;        ;;           ;;        ;;  ;;    ;;  
;;    ;; ;;     ;; ;;        ;;           ;;        ;;  ;;    ;;  
 ;;;;;;   ;;;;;;;  ;;        ;;           ;;       ;;;;  ;;;;;;   

%% Supp figure showing that NL models can't reproduce this, but if we vary the k_D, it can

figure('outerposition',[0 0 1300 901],'PaperUnits','points','PaperSize',[1300 901]); hold on

% show the input nonlinearity 
kD = mean(mean(PID(a:z,paradigm==1)));

x = logspace(-2,2,100);

A = 1./(1+ (kD./x));

subplot(2,3,1); hold on
plot(x,A,'k');
set(gca,'XScale','log')
xlabel('Stimulus (V)')
ylabel('a')

% show the filter
subplot(2,3,2); hold on
filtertime = (1:length(K))*1e-3 - .1;
plot(filtertime,K/norm(K),'k')
set(gca,'XLim',[-.1 .6])
xlabel('Lag (ms)')
ylabel('Filter (norm)')

% generate responses using this model
NL_pred = NaN*fA;
for i = 1:size(fA,2)
	S = PID(:,i);
	S = 1./(1 + kD./S);
	NL_pred(:,i) = convolve(time,S,K,filtertime);
end

% make the i/o plot


subplot(2,3,3); hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = 100*NL_pred(a:z,paradigm == i);
	x = convolve(time,PID(:,i),K,filtertime); x = x(a:z);
	s = PID(a:z,paradigm == i);
	rm_this = sum(y)==0;
	y(:,rm_this) = [];
	x(:,rm_this) = [];
	s(:,rm_this) = [];
	y = nanmean(y,2);
	x = nanmean(x,2);
	s = nanmean(s,2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:),'show_error',false,'LineWidth',3);
end
xlabel('Projected Stimulus (V)')
ylabel('Prediction (Hz)')

% now vary the kD with the mean stimulus

% show the input nonlinearity 

x = logspace(-2,2,100);

subplot(2,3,4); hold on
for i = 1:max(paradigm)
	kD = mean(mean(PID(a:z,paradigm==i)));
	A = 1./(1+ (kD./x));
	plot(x,A,'Color',c(i,:));
end
set(gca,'XScale','log')
xlabel('Stimulus (V)')
ylabel('a')

% show the filter
subplot(2,3,5); hold on
plot(filtertime,K/norm(K),'k')
set(gca,'XLim',[-.1 .6])
xlabel('Lag (ms)')
ylabel('Filter (norm)')

% generate responses using this model
NL_pred = NaN*fA;
for i = 1:size(fA,2)
	kD = mean(mean(PID(a:z,i)));
	S = PID(:,i);
	S = 1./(1 + kD./S);
	NL_pred(:,i) = convolve(time,S,K,filtertime);
end

% make the i/o plot
subplot(2,3,6); hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = 100*NL_pred(a:z,paradigm == i);
	x = convolve(time,PID(:,i),K,filtertime); x = x(a:z);
	s = PID(a:z,paradigm == i);
	rm_this = sum(y)==0;
	y(:,rm_this) = [];
	x(:,rm_this) = [];
	s(:,rm_this) = [];
	y = nanmean(y,2);
	x = nanmean(x,2);
	s = nanmean(s,2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:),'show_error',false,'LineWidth',3);
end
xlabel('Projected Stimulus (V)')
ylabel('Prediction (Hz)')

prettyFig();
labelFigure('x_offset',0,'font_size',28)


if being_published
	snapnow
	delete(gcf)
end


%% Supp Fig for eLIFE

figure('outerposition',[0 0 802 801],'PaperUnits','points','PaperSize',[802 801]); hold on
clear ax

p = find(paradigm == 1); p(5) = []; 
c = [.7 .7 .7];

ax(1) = subplot(4,1,1); hold on
plot(time,PID(:,p),'Color',c)
plot(time,mean(PID(:,p),2),'k','LineWidth',2)
set(gca,'XLim',[35 45.03],'YLim',[0 0.55])
ylabel('Stimulus (V)')

ax(2) = subplot(4,1,2); hold on
plot(time,LFP(:,p),'Color',c)
plot(time,mean(LFP(:,p),2),'k','LineWidth',2)
set(gca,'XLim',[35 45.03],'YLim',[-4 3])
ylabel('ab3 LFP (mV)')

ax(3) = subplot(4,1,3); hold on
y = bandPass(LFP(:,p(end)),10,Inf);
plot(time,y,'k')
set(gca,'XLim',[35 45.03],'YLim',[-.2 .3])
ylabel('Filtered LFP')

ax(4) = subplot(4,1,4); hold on
raster2(A_spikes(p,:),[],0.5,'k')
set(gca,'XLim',[35 45.03],'YLim',[.5 7.5])
ylabel('ORN #')
xlabel('Time (s)')

prettyFig('plw',1);

labelFigure('x_offset',-.1)

for i = 1:4
	deintersectAxes(ax(i))
end

if being_published
	snapnow
	delete(gcf)
end

%  ######  ##     ## ########  ########     ######## ####  ######   
% ##    ## ##     ## ##     ## ##     ##    ##        ##  ##    ##  
% ##       ##     ## ##     ## ##     ##    ##        ##  ##        
%  ######  ##     ## ########  ########     ######    ##  ##   #### 
%       ## ##     ## ##        ##           ##        ##  ##    ##  
% ##    ## ##     ## ##        ##           ##        ##  ##    ##  
%  ######   #######  ##        ##           ##       ####  ######   

clear ax
figure('outerposition',[0 0 1100 800],'PaperUnits','points','PaperSize',[1100 800]); hold on 
for i = 6:-1:1
	ax(i) = subplot(2,3,i); hold on
end

plot(ax(1),nanmean(PID(a:z,:)),nanstd(PID(a:z,:)),'k+')
plot(ax(1),[0 2],[0 2],'k--')
xlabel(ax(1),'\mu_{stimulus} (V)')
ylabel(ax(1),'\sigma_{stimulus (V)}')
set(ax(1),'XLim',[0 2],'YLim',[0 2])
title(ax(1),'ethyl acetate')


% also estimate gain using variances of stimulus and response
mean_stim = nanmean(PID(a:z,:));
frac_var = NaN(width(PID),1);
for i = 1:width(PID)
	try
		frac_var(i) = std(fA(a:z,i))/std(PID(a:z,i));
	catch
	end
end
% for i = 1:width(PID)
% 	plot(ax(4),mean_stim(i),frac_var(i),'+','Color',c(paradigm(i),:))
% end
% fit a power law with exponent -1

% frac_var = frac_var(:);
% options = fitoptions(fittype('power1'));
% options.Lower = [-Inf -1];
% options.Upper = [Inf -1];
% cf = fit(nonnans(mean_stim),nonnans(frac_var),'power1',options);
% plot(ax(4),sort(mean_stim),cf(sort(mean_stim)),'r');
% set(ax(4),'XScale','log','YScale','log','XLim',[.1 3])
% xlabel(ax(4),'Mean stimulus (V)')
% ylabel(ax(4),'\sigma_{Firing rate}/\sigma_{Stimulus} (Hz/V)')
% title(ax(4),['ab3A' char(10) 'ethyl acetate'])

% we're not going to use ratio of std. devs
cla(ax(4))

% instead, plot single-filter gain estimation 

% compute a single filter
temp1 = PID(:,paradigm==1);
temp2 = fA(:,paradigm==1);
temp1(:,5) = [];
temp2(:,5) = [];
K = fitFilter2Data(nanmean(temp1(30e3:55e3,:),2),nanmean(temp2(30e3:55e3,:),2),'offset',200,'reg',1);
K = K(100:end-101);
filtertime = 1e-3*((1:length(K)) - 100);

% re-project and re-calculate gain
single_K_pred = fA_pred*NaN;
single_K_gain = fA_gain*NaN;

for i = 1:length(orn)
	single_K_pred(:,i) = convolve(time,PID(:,i),K,filtertime);
	[ff,gof] = fit(single_K_pred(35e3:55e3,i),fA(35e3:55e3,i),'poly1');
	single_K_gain(i) = ff.p1;
end
% for i = 1:length(orn)
% 	temp = K2(:,i);
% 	temp = temp/norm(temp);
% 	filtertime = 1e-3*((1:length(temp)) - 100);
% 	single_K_pred(:,i) = convolve(time,PID(:,i),temp,filtertime);
% 	[ff,gof] = fit(single_K_pred(45e3:55e3,i),fA(45e3:55e3,i),'poly1');
% 	single_K_gain(i) = ff.p1;
% end


mean_stim = mean(PID);
rm_this = isnan(single_K_gain) | single_K_gain == 0;
mean_stim(rm_this) = NaN;
single_K_gain(rm_this) = NaN;

c = parula(10);
for i = 1:width(PID)
	plot(ax(4),mean_stim(i),single_K_gain(i),'+','Color',c(paradigm(i),:))
end
set(ax(4),'XScale','log','YScale','log')

cf = fit(vectorise(mean_stim),single_K_gain,'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
plot(ax(4),mean_stim,cf(mean_stim),'r')
ax(4).XLim = [.1 10];
ax(4).YLim = [10 1e3];
xlabel(ax(4),'Mean Stimulus (V)')
title(ax(4),'single filter gain estimation')
ylabel(ax(4),'ab3A ORN gain (Hz/V)')


% ##      ## ######## ########  ######## ########   ######  
% ##  ##  ## ##       ##     ## ##       ##     ## ##    ## 
% ##  ##  ## ##       ##     ## ##       ##     ## ##       
% ##  ##  ## ######   ########  ######   ########   ######  
% ##  ##  ## ##       ##     ## ##       ##   ##         ## 
% ##  ##  ## ##       ##     ## ##       ##    ##  ##    ## 
%  ###  ###  ######## ########  ######## ##     ##  ######  

%  ######   ######## ##    ## ######## ########     ###    ##       ##       ##    ## 
% ##    ##  ##       ###   ## ##       ##     ##   ## ##   ##       ##        ##  ##  
% ##        ##       ####  ## ##       ##     ##  ##   ##  ##       ##         ####   
% ##   #### ######   ## ## ## ######   ########  ##     ## ##       ##          ##    
% ##    ##  ##       ##  #### ##       ##   ##   ######### ##       ##          ##    
% ##    ##  ##       ##   ### ##       ##    ##  ##     ## ##       ##          ##    
%  ######   ######## ##    ## ######## ##     ## ##     ## ######## ########    ##    

%  #######  ########   ######  ######## ########  ##     ## ######## ########  
% ##     ## ##     ## ##    ## ##       ##     ## ##     ## ##       ##     ## 
% ##     ## ##     ## ##       ##       ##     ## ##     ## ##       ##     ## 
% ##     ## ########   ######  ######   ########  ##     ## ######   ##     ## 
% ##     ## ##     ##       ## ##       ##   ##    ##   ##  ##       ##     ## 
% ##     ## ##     ## ##    ## ##       ##    ##    ## ##   ##       ##     ## 
%  #######  ########   ######  ######## ##     ##    ###    ######## ########  


% define what we want to work on
data_hashes = {'bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1','f11c4a5792d0c9fec7c40fd6aa2fce40'};
% for i = 1:length(data_hashes)
% 	cdata = consolidateData2(dm.getPath(data_hashes{i}));
% 	[oval(length(cdata.orn)) ' ' oval(length(unique(cdata.orn))) ' ' oval(length(unique(cdata.fly)))]
% end
odour_names = {'1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab2A','ab2A','pb1A'};

plot_here = ax([2 3 5 6]);
% core loop
for i = length(data_hashes):-1:1
	clear cdata
	cdata = consolidateData2(dm.getPath(data_hashes{i}));
	cdata = cleanMSGdata(cdata);

	% plot gain as we normally calculate it
	clear ph
	ph(2) = plot_here(i);
	plotMSGGain(cdata,ph);

	title(ph(2),odour_names{i});

	t = [orn_names{i} ' gain (Hz/V)'];
	ylabel(ph(2),t);
end

% fix some axes
plot_here(4).YLim = [.7e2 .7e4];
for i = 1:length(plot_here)
	plot_here(i).XLim(2) = 1.06*plot_here(i).XLim(2);
end

prettyFig('plw',2,'lw',2,'fs',.6,'font_units','centimeters','FixLogX',true)

labelFigure('font_size',25,'column_first',true);


for i = 2:length(ax)
	try
		deintersectAxes(ax(i));
	catch
	end
end

if being_published
	snapnow
	delete(gcf)
end



return

%% Sanity Check
% In this section, we check if my filter extractions works (comparing it to Damon's FFT-based methods) if we still see the overall phenomenology if we only use a single filter to project all the data. 

% back out Damon's filters here
K_Damon = NaN*K2;
for i = 1:width(PID)
	if (paradigm(i)==1)
		K_Damon(:,i) = backOutFilter(PID(a:z,i),fA(a:z,i),'offset',100,'filter_length',700);
	end
end

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(2,3,1), hold on
for i = 1:max(paradigm)
	temp = nanmean(K2(:,paradigm==i),2);
	plot(filtertime,temp,'Color',c(i,:))
end
xlabel('Filter Lag (s)')
ylabel('Filter amplitude')
title('Filters extracted at different mean stimuli')


subplot(2,3,2), hold on
temp = nanmean(K2(:,paradigm==1),2);
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
	[~,temp] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);
	plot(temp.x,temp.y,'Color',c(i,:));
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
ylabel('\sigma_{Firing rate}/\sigma_{Stimulus} (Hz/V)')

prettyFig('plw',1.3,'lw',1.5,'fs',14,'FixLogX',true)
labelFigure
if being_published
	snapnow
	delete(gcf)
end




%% Version Info
% 
pFooter;
