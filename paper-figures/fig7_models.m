%% fig7_models.m
% 


dm = dataManager;
pHeader;

% ######## ####  ######          #### ##    ## #### ######## 
% ##        ##  ##    ##          ##  ###   ##  ##     ##    
% ##        ##  ##                ##  ####  ##  ##     ##    
% ######    ##  ##   ####         ##  ## ## ##  ##     ##    
% ##        ##  ##    ##          ##  ##  ####  ##     ##    
% ##        ##  ##    ##  ###     ##  ##   ###  ##     ##    
% ##       ####  ######   ###    #### ##    ## ####    ##    


% plots to make
main_fig = figure('outerposition',[0 0 901 902],'PaperUnits','points','PaperSize',[901 902]); hold on
supp_fig = figure('outerposition',[0 0 901 603],'PaperUnits','points','PaperSize',[901 603]); hold on

clear ax
ax.nat_stim = subplot(3,3,1:3,'Parent',main_fig); 
ax.DA_MSG = subplot(3,3,4,'Parent',main_fig); 
ax.r2_DA_LN = subplot(3,3,5,'Parent',main_fig); 

ax.LFP_model_stim = subplot(3,3,7,'Parent',main_fig); 
ax.LFP_model_gain = subplot(3,3,8,'Parent',main_fig); 
ax.LFP_model_resp = subplot(3,3,9,'Parent',main_fig); 

% make insets in the first plot
ax.inset1 = axes('Parent',main_fig); 
ax.inset1.Position = [0.2442    0.8625    0.1154    0.0543];

ax.inset2 = axes('Parent',main_fig); 
ax.inset2.Position = [0.58 0.8625 0.17 0.0543];

% supp plots
ax.DA_lag1 = subplot(2,3,4,'Parent',supp_fig); 
ax.DA_lag2 = subplot(2,3,5,'Parent',supp_fig); 
ax.DA_vs_lin_pred = subplot(2,3,2,'Parent',supp_fig); 
ax.DA_vs_ORN = subplot(2,3,3,'Parent',supp_fig); 
ax.LN_vs_lin_pred = subplot(2,3,1,'Parent',supp_fig); 

% hold all plots
fn = fieldnames(ax);
for i = 1:length(fn)
	ax.(fn{i}).NextPlot = 'add';
end

% define colors, etc. 
c = lines(10);
LFP_color = c(4,:);
firing_color = c(5,:);
model_color = [1 0 0];

%{

 #######  ##       ########     ##    ##    ###    ########     ######  ######## #### ##     ## 
##     ## ##       ##     ##    ###   ##   ## ##      ##       ##    ##    ##     ##  ###   ### 
##     ## ##       ##     ##    ####  ##  ##   ##     ##       ##          ##     ##  #### #### 
##     ## ##       ##     ##    ## ## ## ##     ##    ##        ######     ##     ##  ## ### ## 
##     ## ##       ##     ##    ##  #### #########    ##             ##    ##     ##  ##     ## 
##     ## ##       ##     ##    ##   ### ##     ##    ##       ##    ##    ##     ##  ##     ## 
 #######  ######## ########     ##    ## ##     ##    ##        ######     ##    #### ##     ## 

%}

% first, show the natural stimulus trace together with LN and DA models

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

% fit a LN model to this
ft = fittype('hill4(x,A,K_D,n,x_offset)');
y = mean(fA,2);
x = fp;
rm_this = isnan(fp) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
ff = fit(x(:),y(:),ft,'StartPoint',[50 1 2 0],'Lower',[1 0 1 -10],'Upper',[200 5 10 10],'MaxIter',1e4);

% now plot the nonlinearity vs. the linear prediction, and colour by mean stimulus in the preceding 500ms
LN_pred = ff(fp);
cc = colourByMeanStimulus(mean(PID,2),true(length(PID),1));
scatter(ax.LN_vs_lin_pred,fp,LN_pred,20,cc,'filled')


% fit a DA model to this
clear p
p.   s0 = 0;
p.  n_z = 2;
p.tau_z = 147.3750;
p.  n_y = 2;
p.tau_y = 27.2500;
p.    C = 0.5000;
p.    A = 170.4375;
p.    B = 2.7656;

DA_R = DAModelv2(mean(PID,2),p);

% plot the DA response vs the lin pred and the ORN response
scatter(ax.DA_vs_lin_pred,fp,DA_R,20,cc,'filled')
scatter(ax.DA_vs_ORN,mean(fA,2),DA_R,20,cc,'filled')


% make the first plot
clear l
plot(ax.nat_stim,time(1:10:end),mean(fA,2),'k');
plot(ax.nat_stim,time(1:10:end),DA_R,'r');
plot(ax.nat_stim,time(1:10:end),ff(fp),'b');
l(1) = plot(ax.nat_stim,NaN,NaN,'k','LineWidth',2);
l(2) = plot(ax.nat_stim,NaN,NaN,'r','LineWidth',2);
l(3) = plot(ax.nat_stim,NaN,NaN,'b','LineWidth',2);
legend1 = legend(l,{'Data','DA Model','LN Model'});

plot(ax.inset1,time(1:10:end),mean(fA,2),'k')
plot(ax.inset1,time(1:10:end),DAModelv2(mean(PID,2),p),'r')
plot(ax.inset1,time(1:10:end),ff(fp),'b')

plot(ax.inset2,time(1:10:end),mean(fA,2),'k')
plot(ax.inset2,time(1:10:end),DAModelv2(mean(PID,2),p),'r')
plot(ax.inset2,time(1:10:end),ff(fp),'b')

set(ax.inset1,'Xlim',[8.5 9.5])
set(ax.inset2,'Xlim',[26 29])

ax.r2_DA_LN.NextPlot = 'add';
plot(ax.r2_DA_LN,rsquare(mean(fA,2),ff(fp)),rsquare(mean(fA,2),DAModelv2(mean(PID,2),p)),'k+')

%{

##    ## ######## ##      ##    ##    ##    ###    ########     ######  ######## #### ##     ## 
###   ## ##       ##  ##  ##    ###   ##   ## ##      ##       ##    ##    ##     ##  ###   ### 
####  ## ##       ##  ##  ##    ####  ##  ##   ##     ##       ##          ##     ##  #### #### 
## ## ## ######   ##  ##  ##    ## ## ## ##     ##    ##        ######     ##     ##  ## ### ## 
##  #### ##       ##  ##  ##    ##  #### #########    ##             ##    ##     ##  ##     ## 
##   ### ##       ##  ##  ##    ##   ### ##     ##    ##       ##    ##    ##     ##  ##     ## 
##    ## ########  ###  ###     ##    ## ##     ##    ##        ######     ##    #### ##     ## 

%}

% analyse kinetics of LFP and firing rate during the naturalistic stimulus presentation
load(dm.getPath('aeb361c027b71938021c12a6a12a85cd'),'-mat');

% first, fit the DA model to the data. (this was pre-fit, we're simply loading the fit here)
load(dm.getPath('b1b883899ab8c5ce1aed465819e75fce'));

% generate DA model responses
dd = ORNData;
for i = 1:length(od)
	dd(i).stimulus = od(i).stimulus;
	for j = 1:od(i).n_trials
		dd(i).firing_rate(:,j) = DAModelv2(dd(i).stimulus(:,j),p(i));
	end
end

% fit LN models to this data
LN_od = ORNData;
for i = 1:length(od)
	LN_od(i).stimulus = nanmean(od(i).stimulus,2);
	LN_od(i).firing_rate = nanmean(od(i).firing_rate,2);
	LN_od(i).K_firing = nanmean(od(i).K_firing,2);
end
LN_od = fitNL(LN_od,'firing_rate');


% how well does this do, compared to a LN model?
r2_DA = NaN(length(od),1);
r2_LN = NaN(length(od),1);
for i = 1:length(od)
	resp = od(i).firing_rate;
	rm_this = sum(resp) == 0 | isnan(sum(resp));
	resp(:,rm_this) = [];
	pred = LN_od(i).firing_projected;
	rm_this = sum(pred) == 0;
	pred(:,rm_this) = [];
	r2_LN(i) = rsquare(nanmean(resp,2),nanmean(pred,2));

	pred = dd(i).firing_rate;
	rm_this = sum(pred) == 0;
	pred(:,rm_this) = [];
	r2_DA(i) = rsquare(nanmean(resp,2),nanmean(pred,2));
end

ax.r2_DA_LN.NextPlot = 'add';
plot(ax.r2_DA_LN,[0 1],[0 1],'k--')
plot(ax.r2_DA_LN,r2_LN,r2_DA,'k+')
ylabel(ax.r2_DA_LN,'r^2 (DA Model)')
xlabel(ax.r2_DA_LN,'r^2 (LN Model)')
set(ax.r2_DA_LN,'XLim',[0 1],'YLim',[0 1])
r2_DA

% now show that the lags are constant with the DA model

min_acceptable_corr = .5;
min_acceptable_lag = 2;
clear l
for i = 1:length(od)
	S = nanmean(od(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);
	DA_R = nanmean(dd(i).firing_rate,2);


	mean_x = vectorise(computeSmoothedStimulus(S,200));
	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,1e3,25);

	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% plot
	axes(ax.DA_lag1)
	l(2) = plotPieceWiseLinear(mean_x,lag,'Color',firing_color,'nbins',19);

	axes(ax.DA_lag2)
	plotPieceWiseLinear(time_since_thresh_crossing,lag,'Color',firing_color,'nbins',19);

	% now do the same with the DA model prediction
	mean_x = vectorise(computeSmoothedStimulus(S,200));
	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,DA_R,1e3,25);

	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	axes(ax.DA_lag1)
	l(3) = plotPieceWiseLinear(mean_x,lag,'Color',model_color,'nbins',19);

	axes(ax.DA_lag2)
	plotPieceWiseLinear(time_since_thresh_crossing,lag,'Color',model_color,'nbins',19);

	mean_x = vectorise(computeSmoothedStimulus(S,200));
	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,X,1e3,25);
	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% plot
	axes(ax.DA_lag1)
	l(1) = plotPieceWiseLinear(mean_x,lag,'Color',LFP_color,'nbins',19);

	axes(ax.DA_lag2)
	plotPieceWiseLinear(time_since_thresh_crossing,lag,'Color',LFP_color,'nbins',19);

end

% labels -- main figure
xlabel(ax.DA_lag1,'\mu_{Stimulus} in preceding 200ms (V)')
ylabel(ax.DA_lag1,'Lag (ms)')
set(ax.DA_lag1,'YLim',[0 140],'XLim',[0 0.45])
L = legend(l,{'LFP','Firing Rate','DA Model'},'Location','southeast');
L.Position(1) = .7;

set(ax.DA_lag2,'YLim',[0 140],'XLim',[10 1e4],'XScale','log','XTick',[10 1e2 1e3 1e4])
xlabel(ax.DA_lag2,'Time since odor encounter (ms)')
ylabel(ax.DA_lag2,'Lag (ms)')

% ##     ##  ######   ######      ########     ###    ########    ###    
% ###   ### ##    ## ##    ##     ##     ##   ## ##      ##      ## ##   
% #### #### ##       ##           ##     ##  ##   ##     ##     ##   ##  
% ## ### ##  ######  ##   ####    ##     ## ##     ##    ##    ##     ## 
% ##     ##       ## ##    ##     ##     ## #########    ##    ######### 
% ##     ## ##    ## ##    ##     ##     ## ##     ##    ##    ##     ## 
% ##     ##  ######   ######      ########  ##     ##    ##    ##     ## 


%% fit DA model to data from fig 2
clear cdata
cdata = consolidateData2(dm.getPath('93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);

% % average across paradigms
% data = [];
% for i = 1:max(cdata.paradigm)
% 	s = (cdata.PID(30e3:55e3,cdata.paradigm == i));
% 	f = (cdata.fA(30e3:55e3,cdata.paradigm == i));
% 	rm_this = sum(f) == 0 | isnan(sum(f));
% 	s(:,rm_this) = []; f(:,rm_this) = [];
% 	data(i).stimulus = mean(s,2);
% 	data(i).response = mean(f,2);
% 	data(i).response(1:1e3) = NaN;
% end

clear p
p.   s0 = -0.1503;
p.  n_z = 2;
p.tau_z = 255.3750;
p.  n_y = 2;
p.tau_y = 20.0748;
p.    C = 0.2345;
p.    A = 2.8700e+03;
p.    B = 100;

% now generate responses using the DA model
cdata.DA_R = NaN*cdata.fA;
for i = 1:width(cdata.PID)
	cdata.DA_R(:,i) = DAModelv2(cdata.PID(:,i),p);
end

% fit a nonlinearity to the linear predictions
a = 25e3;
z = 55e3;
cdata.LN_pred = NaN*cdata.fA_pred;
for i = 1:length(cdata.paradigm)
	y = cdata.fA(a:z,i);
	x = cdata.fA_pred(a:z,i);
	s = cdata.PID(a:z,i);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	cdata.LN_pred(a:z,i) = x;
end
x = vectorise(cdata.LN_pred(a:z,:));
y = vectorise(cdata.fA(a:z,:));
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];

ft = fittype('hill4(x,A,K_D,n,x_offset)');
ff = fit(x(1:100:end),y(1:100:end),ft,'StartPoint',[50 1 2 0],'Lower',[1 0 1 -10],'Upper',[100 5 10 10],'MaxIter',1e4);

for i = 1:length(cdata.paradigm)
	cdata.LN_pred(:,i) = ff(cdata.LN_pred(:,i));
end


clear temp
r2_DA = NaN(max(cdata.orn),1);
r2_LN = NaN(max(cdata.orn),1);
for i = 1:max(cdata.orn)
	clear temp
	for j = 1:max(cdata.paradigm)
		fp = (cdata.DA_R(30e3:55e3,cdata.orn == i & cdata.paradigm == j));
		f = (cdata.fA(30e3:55e3,cdata.orn == i & cdata.paradigm == j));
		if size(f,2) > 0
			rm_this = isnan(sum(f)) | sum(f) == 0;
			fp(:,rm_this) = []; f(:,rm_this) = [];
			fp(fp<0) = NaN;
			[~,temp(j)] = plotPieceWiseLinear(mean(fp,2),mean(f,2),'Color',c(i,:),'nbins',50,'make_plot',false);
		end
	end
	r2_DA(i) = rsquare(vectorise([temp.x]),vectorise([temp.y]));

	clear temp
	for j = 1:max(cdata.paradigm)
		fp = (cdata.LN_pred(30e3:55e3,cdata.orn == i & cdata.paradigm == j));
		f = (cdata.fA(30e3:55e3,cdata.orn == i & cdata.paradigm == j));
		if size(f,2) > 0
			rm_this = isnan(sum(f)) | sum(f) == 0;
			fp(:,rm_this) = []; f(:,rm_this) = [];
			fp(fp<0) = NaN;
			[~,temp(j)] = plotPieceWiseLinear(mean(fp,2),mean(f,2),'Color',c(i,:),'nbins',50,'make_plot',false);
		end
	end
	r2_LN(i) = rsquare(vectorise([temp.x]),vectorise([temp.y]));
end


% make plot of DA model predictions vs. data, grouped by paradigm
c = parula(10);
for j = 1:max(cdata.paradigm)
	fp = (cdata.DA_R(30e3:55e3, cdata.paradigm == j));
	f = (cdata.fA(30e3:55e3, cdata.paradigm == j));
	if size(f,2) > 0
		rm_this = isnan(sum(f)) | sum(f) == 0;
		fp(:,rm_this) = []; f(:,rm_this) = [];
		fp(fp<0) = NaN;
		axes(ax.DA_MSG)
		[plot_handles,temp(j)] = plotPieceWiseLinear(mean(fp,2),mean(f,2),'Color',c(j,:),'nbins',50,'make_plot',true);
		delete(plot_handles.line(2:3))
		delete(plot_handles.shade)
	end
end

% now compare r2 for the DA and LN model
c = lines(5);
ax.r2_DA_LN.NextPlot = 'add';
plot(ax.r2_DA_LN,r2_LN,r2_DA,'o','Color',c(2,:))
r2_DA


% ##       ######## ########     ##     ##  #######  ########  ######## ##       
% ##       ##       ##     ##    ###   ### ##     ## ##     ## ##       ##       
% ##       ##       ##     ##    #### #### ##     ## ##     ## ##       ##       
% ##       ######   ########     ## ### ## ##     ## ##     ## ######   ##       
% ##       ##       ##           ##     ## ##     ## ##     ## ##       ##       
% ##       ##       ##           ##     ## ##     ## ##     ## ##       ##       
% ######## ##       ##           ##     ##  #######  ########  ######## ######## 


% now the LFP model

% core parameters 
T = 20e3; % total length
pulse_on = 15e3;
pulse_off = 16e3;

clear p
p.    A = 1e3;
p.    B = 100;
p.tau_y = 20;
p.tau_z = 400;
p.tau_r = 10;
p.tau_A = 100;


background_levels = logspace(-2,2,10);
c = parula(length(background_levels)+1);

% make sure we don't get negative responses
p.ko = -log(background_levels(1))*1.1;

all_gain = NaN*background_levels;
all_mu = NaN*background_levels;
for i = 1:length(background_levels)
	S = background_levels(i) + zeros(T,1);
	S(pulse_on:pulse_off) = 2*background_levels(i);
	time = 1e-3*(1:length(S)) - 10;
	plot(ax.LFP_model_stim,time,S,'Color',c(i,:));

	R = LFPmodel2(S,p);
	% plot(ax(2),R,'Color',c(i,:))

	temp = R - mean(R(pulse_on-1e3:pulse_on-1));

	plot(ax.LFP_model_resp,time,temp/max(temp(pulse_on:pulse_off)),'Color',c(i,:))


	% compute gain 
	g = max(temp(pulse_on:pulse_off))/(max(S) - min(S));
	plot(ax.LFP_model_gain,min(S),g,'+','Color',c(i,:));
	all_gain(i) = g;
	all_mu(i) = max(S) - min(S);
end


set(ax.LFP_model_stim,'XLim',[0 10])
set(ax.LFP_model_resp,'XLim',[0 10])

% draw a weber's fit to this
x = all_mu;
y = all_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax.LFP_model_gain,sort(x),cf(sort(x)),'r')

set(ax.LFP_model_stim,'YScale','log')
ylabel(ax.LFP_model_stim,'Stimulus (a.u.)')
ylabel(ax.LFP_model_resp,'r(t) (rescaled)')
set(ax.LFP_model_resp,'YLim',[-.1 1])
set(ax.LFP_model_gain,'XScale','log','YScale','log','XTick',logspace(-2,2,5),'YTick',logspace(-2,2,5))
xlabel(ax.LFP_model_gain,'Background S')
ylabel(ax.LFP_model_gain,'Gain (\DeltaR / \DeltaS)')



prettyFig(main_fig,'fs',12,'plw',1,'lw',1.5)


% cosmetic fixes
ax.nat_stim.XLim = [0 70];
ax.nat_stim.YLim(1) = -5;
ax.inset1.YLim = [-5 70];
ax.inset2.YLim = [-5 100];
xlabel(ax.nat_stim,'Time (s)');
ylabel(ax.nat_stim,'Firing rate (Hz)');

ax.r2_DA_LN.YTick = [0 .25 .5 .75 1];
ax.r2_DA_LN.XTick = [0 .25 .5 .75 1];


ylabel(ax.DA_MSG,'Firing rate (Hz)');
xlabel(ax.DA_MSG,'DA Model prediction (Hz)');
ax.DA_MSG.XLim(1) = 0;
ax.DA_MSG.YLim(1) = 0;

ax.LFP_model_stim.YMinorTick = 'off';
ax.LFP_model_gain.YMinorTick = 'off';
ax.LFP_model_gain.XMinorTick = 'off';
xlabel(ax.LFP_model_stim,'Time (s)')
xlabel(ax.LFP_model_resp,'Time (s)')

legend1.Box = 'off';

% now supp figures
xlabel(ax.LN_vs_lin_pred,'Projected Stimulus (V)')
ylabel(ax.LN_vs_lin_pred,'LN model prediction (Hz)')
ax.LN_vs_lin_pred.YLim = [0 100];
ax.LN_vs_lin_pred.XLim = [0 6];

xlabel(ax.DA_vs_lin_pred,'Projected Stimulus (V)')
ylabel(ax.DA_vs_lin_pred,'DA model prediction (Hz)')
ax.DA_vs_lin_pred.YLim = [0 120];
ax.DA_vs_lin_pred.XLim = [0 6];

xlabel(ax.DA_vs_ORN,'ORN firing rate (Hz)')
ylabel(ax.DA_vs_ORN,'DA model prediction (Hz)')
ax.DA_vs_ORN.YLim = [0 120];
ax.DA_vs_ORN.XLim = [0 120];

prettyFig(supp_fig,'fs',12,'plw',1,'lw',1.5)


if being_published	
	snapnow	
	delete(gcf)
end

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


