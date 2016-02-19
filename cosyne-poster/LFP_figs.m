% LFP_figs.m
% makes LFP figures for the cosyne poster
%
% created by Srinivas Gorur-Shandilya at 1:58 , 23 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

%% Mechanism of Gain Control: Local Field Potential 

% make figure placeholders 
fig_handle=figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on

axes_handles(1) = subplot(2,2,1); hold on
axes_handles(2) = subplot(2,4,5); hold on
axes_handles(3) = subplot(2,4,6); hold on

axes_handles(4) = subplot(2,2,2); hold on
axes_handles(5) = subplot(2,2,4); hold on

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

% show the LFP paradigm-wise
c = parula(max(paradigm)+1);
time = 1e-3*(1:length(PID));

% plot the LFP and the distributions of the LFP
for i = 1:max(paradigm)
	temp = filtered_LFP((15e3:45e3),paradigm==i);
	plot(axes_handles(1),time(15e3:45e3),nanmean(temp,2),'Color',c(i,:));
	[hy,hx] = histcounts(nanmean(temp,2),-3:.1:3);
	hx = hx(1:end-1) + mean(diff(hx));
	hy = hy/sum(hy);
	plot(axes_handles(2),hx,hy,'Color',c(i,:));
end


% extract filters and compute gains
[~,~,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% plot gain vs. mean stimulus and colour correctly
ms = mean(PID(a:z,:)); ms = ms(:);
for i = 1:max(paradigm)
	plot(axes_handles(3),ms(paradigm==i),LFP_gain(paradigm==i),'+','Color',c(i,:));
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
plot(axes_handles(3),[.2 1.6],ff([.2 1.6]),'r');

set(axes_handles(3),'XScale','log','YScale','log','XLim',[.09 2.5],'YLim',[.9 25])


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

% plot the distributions of the projected stimulus
ax(1) = axes_handles(4); hold on
ax(2) = axes_handles(5); hold on
temp = LFP_pred(1e3:5e3,:); temp = nonnans(temp(:));
x = min(min(LFP_pred)):0.02:max(max(LFP_pred));
y = histcounts(temp,x);
y = y/sum(y); 
plot(ax(1),x(2:end),y,'r','LineWidth',3)

% integrate it and re-plot as a prediction
plot(ax(2),x(2:end),cumsum(y),'r--')

temp = LFP_pred(6e3:end,:); temp = nonnans(temp(:));
y = histcounts(temp,x);
y = y/sum(y);
plot(ax(1),x(2:end),y,'b','LineWidth',3)

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

set(ax,'XLim',[-1.2 -.3])

for i = 1:length(axes_handles)
	try
		xlabel(axes_handles(i),'A A')
		ylabel(axes_handles(i),'A A')
	end
end

prettyFig('fs=18;','FixLogX=true;','lw=2;','plw=1.4;')

set(axes_handles(3),'TickLength',[.03 .02])


% move the plots on the right to the right
movePlot(axes_handles(4),'right',.05)
movePlot(axes_handles(5),'right',.05)


