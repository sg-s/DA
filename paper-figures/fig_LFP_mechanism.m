% fig_LFP_mechanism.m
% 
% created by Srinivas Gorur-Shandilya at 1:58 , 23 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

%% Mechanism of Gain Control: Local Field Potential 

figure('outerposition',[0 0 800 850],'PaperUnits','points','PaperSize',[800 850]); hold on

% ##     ## ########    ###    ##    ## 
% ###   ### ##         ## ##   ###   ## 
% #### #### ##        ##   ##  ####  ## 
% ## ### ## ######   ##     ## ## ## ## 
% ##     ## ##       ######### ##  #### 
% ##     ## ##       ##     ## ##   ### 
% ##     ## ######## ##     ## ##    ## 


%  ######  ##     ## #### ######## ######## ######## ########  
% ##    ## ##     ##  ##  ##          ##    ##       ##     ## 
% ##       ##     ##  ##  ##          ##    ##       ##     ## 
%  ######  #########  ##  ######      ##    ######   ##     ## 
%       ## ##     ##  ##  ##          ##    ##       ##     ## 
% ##    ## ##     ##  ##  ##          ##    ##       ##     ## 
%  ######  ##     ## #### ##          ##    ######## ########  


%  ######      ###    ##     ##  ######   ######  ####    ###    ##    ##  ######  
% ##    ##    ## ##   ##     ## ##    ## ##    ##  ##    ## ##   ###   ## ##    ## 
% ##         ##   ##  ##     ## ##       ##        ##   ##   ##  ####  ## ##       
% ##   #### ##     ## ##     ##  ######   ######   ##  ##     ## ## ## ##  ######  
% ##    ##  ######### ##     ##       ##       ##  ##  ######### ##  ####       ## 
% ##    ##  ##     ## ##     ## ##    ## ##    ##  ##  ##     ## ##   ### ##    ## 
%  ######   ##     ##  #######   ######   ######  #### ##     ## ##    ##  ######  

clear axes_handles


[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


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

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];

% band pass all the LFP
try 
	load('/local-data/DA-paper/LFP-MSG/september/filtered_LFP.mat','filtered_LFP')
catch
	filtered_LFP = LFP;
	for i = 1:width(LFP)
		filtered_LFP(:,i) = 10*bandPass(LFP(:,i),1e4,Inf);
	end
end

% define limits on data
a = 10e3; z = 50e3;

% show the stimulus paradigm-wise
c = parula(max(paradigm)+1);
subplot(3,3,1), hold on
time = 1e-3*(1:length(PID));
time = time(a:z);
for i = 1:max(paradigm)
	plot(time,nanmean(PID(a:z,paradigm==i),2),'Color',c(i,:));
end
set(gca,'XLim',[30 35])
ylabel('Stimulus (V)')
xlabel('Time (s)')

% show the LFP paradigm-wise
c = parula(max(paradigm)+1);
subplot(3,3,2), hold on
time = 1e-3*(1:length(PID));
time = time(a:z);
for i = 1:max(paradigm)
	plot(time,nanmean(filtered_LFP(a:z,paradigm==i),2),'Color',c(i,:));
end
set(gca,'XLim',[30 35])
xlabel('Time (s)')
ylabel('\DeltaLFP (mV)')

% extract filters and compute gains
[~,~,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% plot gain vs. mean stimulus and colour correctly
ms = mean(PID(a:z,:)); ms = ms(:);
subplot(3,3,3); hold on
for i = 1:max(paradigm)
	plot(ms(paradigm==i),LFP_gain(paradigm==i),'+','Color',c(i,:));
end


% plot power law fits for the LFP
x = mean(PID(a:z,:)); x = x(:);
y = LFP_gain(:);

xx = NaN(length(unique(paradigm)),1);
yy = NaN(length(unique(paradigm)),1);
ww = NaN(length(unique(paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(paradigm==i));
	yy(i) = mean(y(paradigm==i));
	ww(i) = 1./sem(y(paradigm==i));
end

fo = fitoptions('power1');
fo.Upper = [Inf -1];
fo.Lower = [-Inf -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
plot([.2 1.6],ff([.2 1.6]),'r');

set(gca,'XScale','log','YScale','log','XLim',[.09 2.5],'YLim',[.9 25])
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V')


%      ##     ##    ###    ########  ####    ###    ##    ##  ######  ######## 
%      ##     ##   ## ##   ##     ##  ##    ## ##   ###   ## ##    ## ##       
%      ##     ##  ##   ##  ##     ##  ##   ##   ##  ####  ## ##       ##       
%      ##     ## ##     ## ########   ##  ##     ## ## ## ## ##       ######   
%       ##   ##  ######### ##   ##    ##  ######### ##  #### ##       ##       
%        ## ##   ##     ## ##    ##   ##  ##     ## ##   ### ##    ## ##       
%         ###    ##     ## ##     ## #### ##     ## ##    ##  ######  ######## 
     
%      ##       ######## ########      ######      ###    #### ##    ## 
%      ##       ##       ##     ##    ##    ##    ## ##    ##  ###   ## 
%      ##       ##       ##     ##    ##         ##   ##   ##  ####  ## 
%      ##       ######   ########     ##   #### ##     ##  ##  ## ## ## 
%      ##       ##       ##           ##    ##  #########  ##  ##  #### 
%      ##       ##       ##           ##    ##  ##     ##  ##  ##   ### 
%      ######## ##       ##            ######   ##     ## #### ##    ## 


%      ########  ########  ######  ##     ##    ###    ########  ######## 
%      ##     ## ##       ##    ## ##     ##   ## ##   ##     ## ##       
%      ##     ## ##       ##       ##     ##  ##   ##  ##     ## ##       
%      ########  ######    ######  ######### ##     ## ########  ######   
%      ##   ##   ##             ## ##     ## ######### ##        ##       
%      ##    ##  ##       ##    ## ##     ## ##     ## ##        ##       
%      ##     ## ########  ######  ##     ## ##     ## ##        ######## 
     
%      ########     ###    ########    ###    
%      ##     ##   ## ##      ##      ## ##   
%      ##     ##  ##   ##     ##     ##   ##  
%      ##     ## ##     ##    ##    ##     ## 
%      ##     ## #########    ##    ######### 
%      ##     ## ##     ##    ##    ##     ## 
%      ########  ##     ##    ##    ##     ## 

clearvars -except being_published axes_handles
path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, ~, ~, orn] = consolidateData(path_name,1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	LFP(a:z,i) = LFP(a:z,i) - fastFiltFilt(ones(1e4,1),1e4,LFP(a:z,i));
	LFP(a:z,i) = LFP(a:z,i)*10; % to get the units right, now in mV
end

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_orn(rm_this) = [];

% we are going to calculate only one filter/epoch
sr = 1e3; % sampling rate, Hz
if exist('.cache/VSA_K.mat','file') == 2
	load('.cache/VSA_K.mat','K1')
else
	filter_length = 1000;
	offset = 200;
	K1 = NaN(800,width(reshaped_PID));
	for i = 1:width(reshaped_LFP)
		textbar(i,width(reshaped_PID))
		stim = reshaped_PID(:,i);
		resp = reshaped_LFP(:,i);

		resp(1:1e3) = NaN;
		resp(5e3:6e3)= NaN;

		try
			this_K1 = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K1(:,i) = this_K1(100:end-101);
		catch 
		end
	end
	mkdir('.cache')
	save('.cache/VSA_K.mat','K1')
end

% make linear predictions on the de-trended data using a mean filter averaged over all cases
K1_mean = nanmean(K1,2);
ft = -99:700;
LFP_pred = NaN*reshaped_LFP;
for i = 1:width(reshaped_LFP)
	LFP_pred(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K1_mean,ft);
end

% plot the stimulus
subplot(3,2,3); hold on
time = 1e-3*(1:length(reshaped_PID));
plot(time(1:10:end),reshaped_PID(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel('Time since switch (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 10],'YLim',[0 1.1])
plot([1 5],[1 1],'r','LineWidth',3)
plot([6 10],[1 1],'b','LineWidth',3)

% plot the LFP response
subplot(3,2,5); hold on
time = 1e-3*(1:length(reshaped_LFP));
plot(time(1:10:end),reshaped_LFP(1:10:end,1:2:end),'Color',[.5 .5 .5 .5]);
xlabel('Time since switch (s)')
ylabel('\DeltaLFP (mV)')
set(gca,'XLim',[0 10],'YLim',[-4 4.3])
plot([1 5],[3 3],'r','LineWidth',3)
plot([6 10],[3 3],'b','LineWidth',3)

% plot the distributions of the projected stimulus
ax(1) = subplot(3,2,4); hold on
ax(2) = subplot(3,2,6); hold on
temp = LFP_pred(1e3:5e3,:); temp = nonnans(temp(:));
x = min(min(LFP_pred)):0.02:max(max(LFP_pred));
y = histcounts(temp,x);
y = y/sum(y); 
plot(ax(1),x(2:end),y,'r')

% integrate it and re-plot as a prediction
plot(ax(2),x(2:end),cumsum(y),'r--')

temp = LFP_pred(6e3:end,:); temp = nonnans(temp(:));
y = histcounts(temp,x);
y = y/sum(y);
plot(ax(1),x(2:end),y,'b')

% integrate it and re-plot as a prediction
plot(ax(2),x(2:end),cumsum(y),'b--')


% now plot the actual i/o curve: high contrast
x = LFP_pred(1e3:5e3,:);
y = reshaped_LFP(1e3:5e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; cy(:,rm_this)= [];
[~,data_hi] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);

% low contrast
x = LFP_pred(6e3:9e3,:);
y = reshaped_LFP(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; cy(:,rm_this)= [];
[~,data_lo] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);

% normalise globally
m = min([data_lo.y data_hi.y]);
data_lo.y = data_lo.y - m; data_hi.y = data_hi.y - m;

M = max([data_lo.y data_hi.y]);
data_lo.y = data_lo.y/M; data_hi.y = data_hi.y/M;

% plot data
axes(ax(2))
[line_handle2, shade_handle2] = errorShade(data_hi.x,data_hi.y,data_hi.ye,'Color',[1 0 0],'Shading',.5);
uistack(shade_handle2,'bottom')

% plot data
[line_handle2, shade_handle2] = errorShade(data_lo.x,data_lo.y,data_lo.ye,'Color',[0 0 1],'Shading',.5);
uistack(shade_handle2,'bottom')
xlabel('Projected Stimulus (V)')
ylabel('\Delta LFP (norm)')

ylabel(ax(1),'Probability')
set(ax,'XLim',[-1.2 -.4])


% ########  ##    ## ##    ##    ###    ##     ## ####  ######   ######      #######  ######## 
% ##     ##  ##  ##  ###   ##   ## ##   ###   ###  ##  ##    ## ##    ##    ##     ## ##       
% ##     ##   ####   ####  ##  ##   ##  #### ####  ##  ##       ##          ##     ## ##       
% ##     ##    ##    ## ## ## ##     ## ## ### ##  ##  ##        ######     ##     ## ######   
% ##     ##    ##    ##  #### ######### ##     ##  ##  ##             ##    ##     ## ##       
% ##     ##    ##    ##   ### ##     ## ##     ##  ##  ##    ## ##    ##    ##     ## ##       
% ########     ##    ##    ## ##     ## ##     ## ####  ######   ######      #######  ##  

%  ######      ###    #### ##    ## 
% ##    ##    ## ##    ##  ###   ## 
% ##         ##   ##   ##  ####  ## 
% ##   #### ##     ##  ##  ## ## ## 
% ##    ##  #########  ##  ##  #### 
% ##    ##  ##     ##  ##  ##   ### 
%  ######   ##     ## #### ##    ## 

%  ######   #######  ##    ## ######## ########   #######  ##       
% ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
% ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
% ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
% ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
% ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%  ######   #######  ##    ##    ##    ##     ##  #######  ######## 

% now do the supplementary figure showing dynamics of gain change


% calculate instantaneous gain everywhere (for a subset of the data for now)
x = reshaped_PID(:); x = x(1:1e5);
y = reshaped_LFP(:); y = y(1:1e5);
z = LFP_pred(:); z = z(1:1e5);
[mean_stim,inst_gain,e] = findInstGain(x,y,z,50);

% reshape back
inst_gain = reshape(inst_gain,size(reshaped_PID,1),10);
e = reshape(e,size(reshaped_PID,1),size(reshaped_PID,2));
mean_stim = reshape(mean_stim,size(reshaped_PID,1),size(reshaped_PID,2));
inst_gain(e<.8) = NaN;
inst_gain(inst_gain == 0) = NaN;
inst_gain(1:1e3,:) = NaN;


figure('outerposition',[0 0 500 800],'PaperUnits','points','PaperSize',[500 800]); hold on

% show how gain changes with time
subplot(3,1,1), hold on
errorShade(time,nanmean(inst_gain,2),nanstd(inst_gain'),'Color',[0 0 0]);
xlabel('Time since switch (s)')
ylabel('Inst. Gain (Hz/V)')

% put a timescale on this change by finding the time to half asymptote 
tau_fA = NaN(width(inst_gain),1);
for i = 1:width(inst_gain)
	if mean(isnan(inst_gain(:,i))) < .3
		a = nanmean(inst_gain(1e3:5e3,i));
		z = nanmean(inst_gain(6e3:end,i));
		try
			tau_fA(i) = find(inst_gain(5e3:6e3,i) > a+ (z-a)/2,1,'first');
		catch
		end
	end
end
tau_fA(tau_fA==1) = NaN;


%    ########    ###     ######  ########     ######      ###    #### ##    ## 
%    ##         ## ##   ##    ##    ##       ##    ##    ## ##    ##  ###   ## 
%    ##        ##   ##  ##          ##       ##         ##   ##   ##  ####  ## 
%    ######   ##     ##  ######     ##       ##   #### ##     ##  ##  ## ## ## 
%    ##       #########       ##    ##       ##    ##  #########  ##  ##  #### 
%    ##       ##     ## ##    ##    ##       ##    ##  ##     ##  ##  ##   ### 
%    ##       ##     ##  ######     ##        ######   ##     ## #### ##    ## 

%     ######   #######  ##    ## ######## ########   #######  ##       
%    ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%    ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
%    ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%    ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%    ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%     ######   #######  ##    ##    ##    ##     ##  #######  ######## 


% clearvars -except axes_handles being_published h

% p = '/local-data/DA-paper/large-variance-flicker/LFP/';
% [PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes]  = consolidateData(p,1);


% % set to NaN firing rates that are 0
% fA(:,max(fA) == 0) = NaN;

% % throw out trials where we didn't record the LFP, for whatever reason
% not_LFP = 0*orn;
% for i = 1:width(LFP)
% 	not_LFP(i) = abs(mean2(LFP(:,i)));
% end
% LFP(:,not_LFP< 0.5) = NaN;


% % filter the LFP
% filtered_LFP = LFP;
% for i = 1:width(LFP)
% 	a = find(~isnan(LFP(:,i)),1,'first');
% 	z = find(~isnan(LFP(:,i)),1,'last');
% 	if isempty(a)
% 		a = 1;
% 	end
% 	if isempty(z)
% 		z = length(LFP);
% 	end
% 	try
% 		this_LFP = this_LFP - fastFiltFilt(ones(1e4,1),1e4,LFP(a:z,i));
% 		filtered_LFP(a:z,i) = this_LFP*10; % to get the units right, now in mV
% 	catch
% 	end
% end

% % extract LFP filters trialwise
% a = 10e3; z = 50e3;
% [K1,LFP_pred] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% % also extract firing rate filters
% [K2,fA_pred] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% % fit non-linearities by ORN for the LFP
% for i = 1:max(orn)
% 	do_these = orn == i;
% 	x = nanmean(LFP_pred(a:z,do_these),2);
% 	y = nanmean(filtered_LFP(a:z,do_these),2);
% 	ff = fit(x(:),y(:),'poly4');
% 	do_these = find(do_these);
% 	for j = 1:length(do_these)
% 		do_this = do_these(j);
% 		LFP_pred(:,do_this) = ff(LFP_pred(:,do_this));
% 	end
% end


% % find inst gain for the LFP and the firing rate
% inst_gain_LFP = NaN(length(PID),width(PID));
% inst_gain_fA = NaN(length(PID),width(PID));
% gain_err_LFP = NaN(length(PID),width(PID));
% gain_err_fA = NaN(length(PID),width(PID));


% for i = 1:max(orn)
% 	s = nanmean(PID(:,do_these),2);
% 	p = nanmean(LFP_pred(:,do_these),2);
% 	r = nanmean(filtered_LFP(:,do_these),2);
% 	[~,inst_gain_LFP(:,i),gain_err_LFP(:,i)] = makeFig6G(s,r,p,50);

% 	p = nanmean(fA_pred(:,do_these),2);
% 	r = nanmean(fA(:,do_these),2);
% 	[~,inst_gain_fA(:,i),gain_err_fA(:,i)] = makeFig6G(s,r,p,50);

% end

% % clean up
% for i = 1:max(orn)
% 	inst_gain_fA(gain_err_fA(:,i)<.8,i) = NaN;
% 	inst_gain_fA(inst_gain_fA(:,i)<0,i) = NaN;
% 	inst_gain_fA(1:20e3,i) = NaN;

% 	inst_gain_LFP(gain_err_LFP(:,i)<.8,i) = NaN;
% 	inst_gain_LFP(inst_gain_LFP(:,i)<0,i) = NaN;
% 	inst_gain_LFP(1:20e3,i) = NaN;
% end


% % now vary the gain filter and show it is the best possible
% history_lengths = logspace(-2,1,30);
% rho_fA = NaN(length(history_lengths),max(orn));
% rho_LFP = NaN(length(history_lengths),max(orn));
% for i = 1:max(orn)
% 	do_these = orn == i;
% 	this_stim = nanmean(PID(:,do_these),2);
% 	this_inst_gain = inst_gain_fA(:,i);
% 	for j = 1:length(history_lengths)
% 		hl = floor(history_lengths(j)*1e3);
% 		K = ones(hl,1);
% 		% filter the stimulus
% 		shat = abs(filter(K,sum(K),this_stim));

% 		% remove some junk
% 		temp = this_inst_gain(~isnan(this_inst_gain));
% 		shat = shat(~isnan(this_inst_gain));

% 		% use the Spearman rank correlation
% 		rho_fA(j,i) = (spear(shat(1:100:end),temp(1:100:end)));
% 	end

% 	this_inst_gain = inst_gain_LFP(:,i);
% 	for j = 1:length(history_lengths)
% 		hl = floor(history_lengths(j)*1e3);
% 		K = ones(hl,1);
% 		% filter the stimulus
% 		shat = abs(filter(K,sum(K),this_stim));

% 		% remove some junk
% 		temp = this_inst_gain(~isnan(this_inst_gain));
% 		shat = shat(~isnan(this_inst_gain));

% 		% use the Spearman rank correlation
% 		rho_LFP(j,i) = (spear(shat(1:100:end),temp(1:100:end)));
% 	end
% end

% axes_handles(4) = subplot(2,4,4); hold on
% movePlot(gca,'right',.03)
% c = lines(max(orn)+1);
% plot(history_lengths,rho_LFP,'b+');
% plot(history_lengths,rho_fA,'k.');
% xlabel('History Length (s)')
% ylabel('\rho')
% set(gca,'XScale','log')



% % compute mean_stim on a 200ms window
% mean_stim = NaN*inst_gain_LFP;
% K = ones(300,1);
% for i = 1:max(orn)
% 	mean_stim(:,i) = nanmean(PID(:,orn==i),2);
% 	mean_stim(:,i) = abs(filter(K,length(K),mean_stim(:,i)));
% end


% axes_handles(3) = subplot(2,4,3); hold on
% clear data_fA data_LFP
% for i = 1:max(orn)
% 	 [~,data_fA(i)] = plotPieceWiseLinear(mean_stim(20e3:z,i),inst_gain_fA(20e3:z,i),'nbins',50,'make_plot',false);
% 	 [~,data_LFP(i)] = plotPieceWiseLinear(mean_stim(20e3:z,i),inst_gain_LFP(20e3:z,i),'nbins',50,'make_plot',false);
% end
% clear l L
% [ax,p1,p2]= plotyy([data_fA.x],[data_fA.y],[data_LFP.x],[data_LFP.y]);
% % fix the colours, etc
% for i = 1:max(orn)
% 	set(p1(i),'Color','k','Marker','.','LineStyle','none')
% 	set(p2(i),'Color','b','Marker','+','LineStyle','none')
% end
% set(ax(2),'YColor','b','YLim',[1e-2 1e1],'YScale','log','YTick',[ 1e-2 1e-1 1e0 1e1])
% set(ax(1),'YColor','k','YLim',[1e1 1e4],'YScale','log','YTick',[1e1 1e2 1e3 1e4])
% ylabel(ax(1),'Inst. ORN Gain (Hz/V)')
% ylabel(ax(2),'Inst. LFP Gain (mV/V)')
% xlabel('Mean Stimulus in preceding 200ms (V)')

% clear L
% x = nanmean(mean_stim(20e3:10:z,:),2); 
% y = nanmean(inst_gain_fA(20e3:10:z,:),2);
% x(isnan(y)) = []; y(isnan(y)) = [];
% L{1} = ['\rho = ' oval(spear(x,y))];
% x = nanmean(mean_stim(20e3:10:z,:),2); 
% y = nanmean(inst_gain_LFP(20e3:10:z,:),2);
% x(isnan(y)) = []; y(isnan(y)) = [];
% L{2} = ['\rho = ' oval(spear(x,y))];
% legend([p1(1) p2(2)],L)


% % show how we calculate the inst. gain
% inst_ax(1) = subplot(4,4,1); hold on
% inst_ax(2) = subplot(4,4,5); hold on
% inst_ax(3) = subplot(2,4,2); hold on

% y = nanmean(filtered_LFP(2e4:2.1e4,orn==i),2);
% x = nanmean(LFP_pred(2e4:2.1e4,orn==i),2);
% t = 1:length(x);

% plot(inst_ax(1),t,y,'k','LineWidth',1.1)
% plot(inst_ax(2),t,x,'Color',[.6 .1 .1],'LineWidth',1.1)
% plot(inst_ax(1),t(750:800),y(750:800),'k','LineWidth',3)
% plot(inst_ax(2),t(750:800),x(750:800),'Color',[.6 .1 .1],'LineWidth',3)
% plot(inst_ax(3),x(750:5:800),y(750:5:800),'kd','LineWidth',3)

% plot(inst_ax(1),[800 800],[-2 2],'k','LineWidth',0.5)
% plot(inst_ax(2),[800 800],[-2 2],'k','LineWidth',0.5)

% ff = fit(x(750:800),y(750:800),'poly1');
% plot(inst_ax(3),[min(x(750:800)) max(x(750:800))],ff([min(x(750:800)) max(x(750:800))]),'k','LineWidth',.5)

% set(inst_ax(1),'XTick',[],'XColor','w','YColor','k','Position',[.25 .75 .1 .1])
% set(inst_ax(2),'XTick',[],'XColor','w','YColor',[.6 .1 .1],'Position',[.25 .61 .1 .1])
% set(inst_ax(3),'XColor',[.6 .1 .1],'YColor','k','Position',[.3827 .6768 .08 .16])


% prettyFig('fs=18;','plw=1.5;','lw=1.5;','FixLogX=true;')

% % minor prettification
% ylabel(inst_ax(2),['Projected' char(10) ' Stimulus (V)'],'FontSize',12)
% ylabel(inst_ax(1),'\DeltaLFP (mV)','FontSize',12)

% set(inst_ax(3),'XLim',[1.3 1.7],'YLim',[1.3 1.7],'box','off','XTick',[1.3 1.4 1.5 1.6 ])

% movePlot(axes_handles(3),'left',.025)
% movePlot(axes_handles(7),'right',.03)
% movePlot(axes_handles(5),'right',.05)
% movePlot(axes_handles(6),'right',.05)

% if being_published
% 	snapnow
% 	delete(gcf)
% end


% %% Supplementary Figure: LFP gain control in other sensilla:

% figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% subplot(1,2,1), hold on

% % also show the pb3 data
% [pb3_PID, pb3_LFP, ~, pb3_paradigm] = consolidateData('/local-data/DA-paper/palp/pb3/',1);


% % remove baseline from all PIDs
% for i = 1:width(pb3_PID)
% 	pb3_PID(:,i) = pb3_PID(:,i) - mean(pb3_PID(1:5e3,i));
% end

% % throw out trials where we didn't record the LFP, for whatever reason
% not_LFP = find((max(abs(pb3_LFP))) < 0.1);
% pb3_LFP(:,not_LFP) = NaN;

% % throw our bad traces
% bad_trials =  (isnan(sum(pb3_LFP)));
% pb3_LFP(:,bad_trials) = [];
% pb3_PID(:,bad_trials) = [];
% pb3_paradigm(bad_trials) = [];

% % band pass all the LFP
% for i = 1:width(pb3_LFP)
% 	pb3_LFP(:,i) = bandPass(pb3_LFP(:,i),1e4,10);
% end

% a = 10e3; z = 50e3;
% [~,~,pb3_LFP_gain] = extractFilters(pb3_PID,pb3_LFP,'use_cache',true,'a',a,'z',z);


% x = mean(pb3_PID(a:z,:)); x = x(:);
% y = pb3_LFP_gain(:);

% xx = NaN(length(unique(pb3_paradigm)),1);
% yy = NaN(length(unique(pb3_paradigm)),1);
% ww = NaN(length(unique(pb3_paradigm)),1);

% for i = 1:length(yy)
% 	xx(i) = mean(x(pb3_paradigm==i));
% 	yy(i) = mean(y(pb3_paradigm==i));
% 	ww(i) = 1./sem(y(pb3_paradigm==i));
% end
% ww(isinf(ww)) = max(ww(~isinf(ww)));
% ff = fit(xx(:),yy(:),'power1','Weights',ww);
% plot(sort(x),ff(sort(x)),'k');
% l(2) = plot(nanmean(pb3_PID(a:z,:)),pb3_LFP_gain,'kd');

% fo = fitoptions('power1');
% fo.Upper = [NaN -1];
% fo.Lower = [NaN -1];
% fo.Weights = ww;
% ff = fit(xx(:),yy(:),'power1',fo);
% plot(sort(x),ff(sort(x)),'r');
% set(gca,'XScale','log','YScale','log','XLim',[.05 .5],'XTick',[.05 .1 .5])
% xlabel('amyl acetate stimulus (V)')
% ylabel('pb3 LFP gain (mV/V)')


% % now show ab8 data
% use_cache = 1;
% [PID, LFP, ~, paradigm] = consolidateData('/local-data/obp/ab8/wcs',use_cache);

% % clean up data
% rm_this = isnan(sum(LFP));
% PID(:,rm_this) = [];
% LFP(:,rm_this) = [];
% paradigm(rm_this) = [];

% % bandPass LFP
% filtered_LFP = LFP;
% for i = 1:width(LFP)
% 	filtered_LFP(:,i) = bandPass(LFP(:,i),1e4,10);
% end

% a = 10e3; z = 50e3;
% [~,~,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% subplot(1,2,2), hold on
% x = mean(PID(a:z,:)); x = x(:);
% y = LFP_gain(:);

% xx = NaN(length(unique(paradigm)),1);
% yy = NaN(length(unique(paradigm)),1);
% ww = NaN(length(unique(paradigm)),1);

% for i = 1:length(yy)
% 	xx(i) = mean(x(paradigm==i));
% 	yy(i) = mean(y(paradigm==i));
% 	ww(i) = 1./sem(y(paradigm==i));
% end
% ww(isinf(ww)) = max(ww(~isinf(ww)));
% ff = fit(xx(:),yy(:),'power1','Weights',ww);
% plot(sort(x),ff(sort(x)),'k');
% l(2) = plot(nanmean(PID(a:z,:)),LFP_gain,'ko');

% fo = fitoptions('power1');
% fo.Upper = [NaN -1];
% fo.Lower = [NaN -1];
% fo.Weights = ww;
% ff = fit(xx(:),yy(:),'power1',fo);
% plot(sort(x),ff(sort(x)),'r');

% set(gca,'XScale','log','YScale','log','XLim',[.1 10])
% xlabel('ethyl acetate stimulus (V)')
% ylabel('ab8 LFP gain (mV/V)')

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
pFooter;