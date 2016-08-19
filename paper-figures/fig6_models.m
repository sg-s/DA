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


% make the main figure
main_fig = figure('outerposition',[0 0 901 602],'PaperUnits','points','PaperSize',[901 602]); hold on

clear ax
ax.nat_stim = subplot(2,3,1:3,'Parent',main_fig); 
ax.DA_MSG = subplot(2,3,5,'Parent',main_fig); 
ax.DA_kinetics = subplot(2,3,4); hold on
ax.r2_DA_LN = subplot(2,3,6,'Parent',main_fig); 

% make insets in the first plot
ax.inset1 = axes('Parent',main_fig); 
ax.inset1.Position = [0.22    0.8625    0.15    0.1];

ax.inset2 = axes('Parent',main_fig); 
ax.inset2.Position = [0.58 0.8625 0.17 0.1];

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

 ##    ##    ###    ########     ######  ######## #### ##     ## 
 ###   ##   ## ##      ##       ##    ##    ##     ##  ###   ### 
 ####  ##  ##   ##     ##       ##          ##     ##  #### #### 
 ## ## ## ##     ##    ##        ######     ##     ##  ## ### ## 
 ##  #### #########    ##             ##    ##     ##  ##     ## 
 ##   ### ##     ##    ##       ##    ##    ##     ##  ##     ## 
 ##    ## ##     ##    ##        ######     ##    #### ##     ## 

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
% scatter(ax.LN_vs_lin_pred,fp,LN_pred,20,cc,'filled')


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
% scatter(ax.DA_vs_lin_pred,fp,DA_R,20,cc,'filled')
% scatter(ax.DA_vs_ORN,mean(fA,2),DA_R,20,cc,'filled')


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

set(ax.inset1,'Xlim',[7.5 10.5])
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

cdata = consolidateData2('/local-data/DA-paper/data-for-paper/nat-stim/ab3-2ac');
v2struct(cdata)

% for c = 1:5
% 	for i = 1:max(orn)
% 		clear data
% 		R = fA(:,orn==i);
% 		X = LFP(:,orn==i);
% 		S = PID(:,orn==i);
% 		rm_this = sum(R) == 0 | isnan(sum(X));
% 		R(:,rm_this) = [];
% 		S(:,rm_this) = [];
% 		X(:,rm_this) = [];
% 		data.response = nanmean(R,2);
% 		data.response(1:5e3) = NaN;
% 		data.stimulus = nanmean(S,2);
% 		data.stimulus = data.stimulus - mean(data.stimulus(1:5e3));
% 		p(i) = fitModel2Data(@DAModelv2,data,'make_plot',true,'nsteps',10);
% 	end
% end
% save('/local-data/DA-paper/data-for-paper/fig7/DA_Model_fit_to_naturalistic_data.mat','p')


% first, fit the DA model to the data. (this was pre-fit, we're simply loading the fit here)
load('/local-data/DA-paper/data-for-paper/fig7/DA_Model_fit_to_naturalistic_data.mat','p')


% generate DA model responses
dd = ORNData;
for i = 1:max(orn)
	S = PID(:,orn==i);
	dd(i).stimulus = nanmean(S,2);
	dd(i).stimulus = dd(i).stimulus - mean(dd(i).stimulus(1:5e3));
	dd(i).firing_rate = DAModelv2(dd(i).stimulus,p(i));
end

% fit LN models to this data
load('/local-data/DA-paper/data-for-paper/nat-stim/v3_nat_stim.ORNData','-mat')
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

% ########     ###       ##     ##  #######  ########  ######## ##       
% ##     ##   ## ##      ###   ### ##     ## ##     ## ##       ##       
% ##     ##  ##   ##     #### #### ##     ## ##     ## ##       ##       
% ##     ## ##     ##    ## ### ## ##     ## ##     ## ######   ##       
% ##     ## #########    ##     ## ##     ## ##     ## ##       ##       
% ##     ## ##     ##    ##     ## ##     ## ##     ## ##       ##       
% ########  ##     ##    ##     ##  #######  ########  ######## ######## 

% ######## #### ########  #### ##    ##  ######      ##          ###     ######   
% ##        ##  ##     ##  ##  ###   ## ##    ##     ##         ## ##   ##    ##  
% ##        ##  ##     ##  ##  ####  ## ##           ##        ##   ##  ##        
% ######    ##  ########   ##  ## ## ## ##   ####    ##       ##     ## ##   #### 
% ##        ##  ##   ##    ##  ##  #### ##    ##     ##       ######### ##    ##  
% ##        ##  ##    ##   ##  ##   ### ##    ##     ##       ##     ## ##    ##  
% ##       #### ##     ## #### ##    ##  ######      ######## ##     ##  ######   


min_acceptable_corr = .8;
min_acceptable_lag = 20;
window_size = 1e3;
stim_thresh = .035;
clear l

axes(ax.DA_kinetics)
for i = 2:length(od) % first one has stimulus which is too low, different paradigm 
	textbar(i,length(od))
	S = nanmean(dd(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(dd(i).firing_rate,2);

	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,window_size,25);

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag < min_acceptable_lag | max_corr < min_acceptable_corr | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];

	% plot

	l(2) = plotPieceWiseLinear(mean_x,lag,'Color',firing_color,'nbins',19);
end

set(ax.DA_kinetics,'XScale','log','YLim',[0 160],'XLim',[1e-3 5])



% ##       ######## ########     ##     ##  #######  ########  ######## ##       
% ##       ##       ##     ##    ###   ### ##     ## ##     ## ##       ##       
% ##       ##       ##     ##    #### #### ##     ## ##     ## ##       ##       
% ##       ######   ########     ## ### ## ##     ## ##     ## ######   ##       
% ##       ##       ##           ##     ## ##     ## ##     ## ##       ##       
% ##       ##       ##           ##     ## ##     ## ##     ## ##       ##       
% ######## ##       ##           ##     ##  #######  ########  ######## ######## 

% fit modified DA model to LFP

% for c = 1:5
% 	for i = 1:max(orn)
% 		clear data
% 		R = fA(:,orn==i);
% 		X = LFP(:,orn==i);
% 		S = PID(:,orn==i);
% 		rm_this = sum(R) == 0 | isnan(sum(X));
% 		R(:,rm_this) = [];
% 		S(:,rm_this) = [];
% 		X(:,rm_this) = [];
% 		data.response = -nanmean(X,2);

% 		data.response(1:5e3) = NaN;
% 		data.stimulus = nanmean(S,2);
% 		data.stimulus = data.stimulus - mean(data.stimulus(1:5e3));
% 		q(i) = orderfields(fitModel2Data(@LFPmodel2_Euler,data,'make_plot',false,'nsteps',100,'p0',q(i)));
% 	end
% end
% save('/local-data/DA-paper/data-for-paper/fig7/LFP_Model_fit_to_naturalistic_data.mat','q')

% generate responses using LFP model
clear q
for i = 1:7
	q(i).tau_A = 1e4;
	q(i).A = 110;
	q(i).B = 10;
	q(i).tau_r = 2;
	q(i).s0 = 0;
	q(i).tau_y = 20;
	q(i).tau_z = 200;
end


dd = ORNData;
for i = 1:max(orn)
	S = PID(:,orn==i);
	dd(i).stimulus = nanmean(S,2);
	dd(i).stimulus = dd(i).stimulus - mean(dd(i).stimulus(1:5e3));
	dd(i).LFP = LFPmodel2_Euler(dd(i).stimulus,q(i));
end

min_acceptable_corr = .7;
min_acceptable_lag = 20;
window_size = 1e3;
stim_thresh = .035;


axes(ax.DA_kinetics)
for i = 2:length(od) % first one has stimulus which is too low, different paradigm 
	textbar(i,length(od))
	S = nanmean(dd(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(dd(i).LFP,2);

	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,window_size,25);

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag < min_acceptable_lag | max_corr < min_acceptable_corr | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];

	% plot

	l(1) = plotPieceWiseLinear(mean_x,lag,'Color',LFP_color,'nbins',19);
end

set(ax.DA_kinetics,'XScale','log','YLim',[0 160],'XLim',[1e-3 5],'XTick',[1e-3 1e-2 1e-1 1])
xlabel(ax.DA_kinetics,'\mu_{Stimulus} in preceding 1 s (V)')
ylabel(ax.DA_kinetics,'Lag (ms)')
legend2 = legend(l,{'LFP DA model','Firing DA model'},'Location','northwest');


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


% make plot of DA model predictions vs. linear prediction, grouped by paradigm
axes(ax.DA_MSG)
c = parula(10);
for j = 1:max(cdata.paradigm)
	fp = cdata.fA_pred(30e3:55e3, cdata.paradigm == j);
	f = cdata.fA(30e3:55e3, cdata.paradigm == j);
	s = cdata.PID(30e3:55e3, cdata.paradigm == j);
	if size(f,2) > 0
		rm_this = isnan(sum(f)) | sum(f) == 0;
		fp(:,rm_this) = []; f(:,rm_this) = [];
		fp(fp<0) = NaN;
		fp = nanmean(fp,2); 
		s = nanmean(s,2);
		fp = fp - nanmean(fp);
		fp = fp + nanmean(s);
		[plot_handles,temp(j)] = plotPieceWiseLinear(fp,nanmean(f,2),'Color',c(j,:),'nbins',50,'make_plot',true);
		delete(plot_handles.line(2:3))
		delete(plot_handles.shade)
	end
end

% now compare r2 for the DA and LN model
c = lines(5);
ax.r2_DA_LN.NextPlot = 'add';
plot(ax.r2_DA_LN,r2_LN,r2_DA,'o','Color',c(2,:))
r2_DA



% ##     ##    ###    ########  ####    ###    ##    ##  ######  ######## 
% ##     ##   ## ##   ##     ##  ##    ## ##   ###   ## ##    ## ##       
% ##     ##  ##   ##  ##     ##  ##   ##   ##  ####  ## ##       ##       
% ##     ## ##     ## ########   ##  ##     ## ## ## ## ##       ######   
%  ##   ##  ######### ##   ##    ##  ######### ##  #### ##       ##       
%   ## ##   ##     ## ##    ##   ##  ##     ## ##   ### ##    ## ##       
%    ###    ##     ## ##     ## #### ##     ## ##    ##  ######  ######## 


[PID, LFP, fA, paradigm, orn, fly] = consolidateData(dm.getPath('e30707e8e8ef6c0d832eee31eaa585aa'),1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
% for i = 1:width(LFP)
% 	a = find(~isnan(LFP(:,i)),1,'first');
% 	z = find(~isnan(LFP(:,i)),1,'last');
% 	LFP(a:z,i) = bandPass(LFP(a:z,i),1000,10)*10; % now in mV
% end

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% reshape the firing rate signals
reshaped_fA = fA(global_start:end-1e4-1,1:width(PID));
reshaped_fA = reshape(reshaped_fA,block_length,width(reshaped_fA)*length(reshaped_fA)/block_length);


% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];


% we are going to calculate only one filter/epoch
sr = 1e3; % sampling rate, Hz
if exist('.cache/VSA_K2.mat','file') == 2
	load('.cache/VSA_K2.mat','K2')
else
	filter_length = 1000;
	offset = 200;
	K2 = NaN(2,filter_length-offset,width(reshaped_fA));
	for i = 1:width(reshaped_fA)
		textbar(i,width(reshaped_PID))

		% calculate filter for large variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_fA(:,i);

		resp(1:1e3) = NaN;
		resp(5e3:end)= NaN;

		try
			this_K2 = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K2(1,:,i) = this_K2(100:end-101);
		catch 
		end

		% calculate filter for low variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_fA(:,i);

		resp(1:6e3) = NaN;

		try
			this_K2 = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K2(2,:,i) = this_K2(100:end-101);
		catch 
		end
	end
	mkdir('.cache')
	save('.cache/VSA_K2.mat','K2')
end


% make linear predictions on the de-trended data using a mean filter averaged over all cases
K2_mean = nanmean(squeeze(nanmean(K2,1)),2);
ft = -99:700;
fA_pred = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	fA_pred(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K2_mean,ft);
end

% fit a DA Model to this data
% clear data c p
% for i = 1:max(orn)
% 	c = find(orn==i,1,'last');
% 	data(i).stimulus = PID(2e4:10e4,c);
% 	data(i).response = fA(2e4:10e4,c);
% 	data(i).response(1:5e3) = NaN;
% end

% clear p
% for j = 1:5
% 	for i = [1 3 4 5]
% 		p(i) = fitModel2Data(@DAModelv2,data(i),'make_plot',true,'nsteps',1);
% 	end
% end

% we've already fit the model, just load the data
clear p
load('.cache/DA_model_fit_to_variance_data.mat','p')

% make DA model prediction on a per-neuron basis
clear data 
for i = [1 3 4 5]
	c = find(orn==i,1,'last');
	data(i).stimulus = PID(2e4:10e4,c);
	data(i).response = fA(2e4:10e4,c);
	
	data(i).DA_pred = DAModelv2(data(i).stimulus,p(i));
end

% fit a LN model to this
for i = [1 3 4 5]
	ft = -99:700;
	data(i).linear_pred = convolve(1e-3*(1:length(data(i).stimulus)),data(i).stimulus,K2_mean,ft);
	ft = fittype('hill4(x,A,K_D,n,x_offset)');
	x = data(i).linear_pred;
	y = data(i).response;
	rm_this = isnan(x)| isnan(y);
	ff = fit(x(~rm_this),y(~rm_this),ft,'StartPoint',[max(data(i).response) 2 2 0],'Lower',[0 0 1 -1],'Upper',[200 2 10 1],'MaxIter',1e4);
	data(i).LN_pred = ff(data(i).linear_pred);
end

% now plot r2 b/w DA and LN models for the contrast
for i = [1 3 4 5]
	r2_DA = rsquare(data(i).DA_pred,data(i).response);
	r2_LN = rsquare(data(i).LN_pred,data(i).response);
	plot(ax.r2_DA_LN,r2_LN,r2_DA,'bd');
end

% now fake some plots for a nice legend in the r2 plot
clear l 
l(1) = plot(ax.r2_DA_LN,NaN,NaN,'k+');
c = lines(5);
l(2) = plot(ax.r2_DA_LN,NaN,NaN,'o','Color',c(2,:));
l(3) = plot(ax.r2_DA_LN,NaN,NaN,'bd');
legend3 = legend(l,{'Nat. Stim.','changing mean','changing variance'});
legend3.Position = [0.7913    0.1512    0.1465    0.0775];

prettyFig(main_fig,'fs',12,'plw',1.5,'lw',1.5)

labelFigure('ignore_these',[ax.inset1 ax.inset2])

% cosmetic fixes
ax.nat_stim.XLim = [0 70];
ax.nat_stim.YLim(1) = -5;
ax.inset1.YLim = [-5 70];
ax.inset2.YLim = [-5 100];
xlabel(ax.nat_stim,'Time (s)');
ylabel(ax.nat_stim,'Firing rate (Hz)');

ax.r2_DA_LN.YTick = [0 .25 .5 .75 1];
ax.r2_DA_LN.XTick = [0 .25 .5 .75 1];


xlabel(ax.DA_MSG,'Projected Stimulus (V)');
ylabel(ax.DA_MSG,'DA Model prediction (Hz)');
ax.DA_MSG.XLim(1) = 0;
ax.DA_MSG.YLim(1) = 0;



if being_published	
	snapnow	
	delete(gcf)
end

% if being_published	
% 	snapnow	
% 	delete(gcf)
% end


%% Version Info
%
pFooter;


