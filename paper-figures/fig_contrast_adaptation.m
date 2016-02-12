% fig_contrast_adaptation.m
% makes figure: contrast adaptation in ORNs
% 
% created by Srinivas Gorur-Shandilya at 7:10 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

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


path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, paradigm, orn] = consolidateData(path_name,1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	LFP(a:z,i) = bandPass(LFP(a:z,i),1000,10)*10; % now in mV
end

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

% 

% ######## #### ########  #### ##    ##  ######   
% ##        ##  ##     ##  ##  ###   ## ##    ##  
% ##        ##  ##     ##  ##  ####  ## ##        
% ######    ##  ########   ##  ## ## ## ##   #### 
% ##        ##  ##   ##    ##  ##  #### ##    ##  
% ##        ##  ##    ##   ##  ##   ### ##    ##  
% ##       #### ##     ## #### ##    ##  ######   



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



% ##     ##    ###    ##    ## ########    ########  ##        #######  ######## 
% ###   ###   ## ##   ##   ##  ##          ##     ## ##       ##     ##    ##    
% #### ####  ##   ##  ##  ##   ##          ##     ## ##       ##     ##    ##    
% ## ### ## ##     ## #####    ######      ########  ##       ##     ##    ##    
% ##     ## ######### ##  ##   ##          ##        ##       ##     ##    ##    
% ##     ## ##     ## ##   ##  ##          ##        ##       ##     ##    ##    
% ##     ## ##     ## ##    ## ########    ##        ########  #######     ##    


time = 1e-3*(1:length(reshaped_PID));
s = .5; % shading opacity

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,2,1), hold on
plot(time(1:10:end),reshaped_PID(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel('Time since switch (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 10],'YLim',[0 1.1])
plot([1 5],[1 1],'r','LineWidth',3)
plot([6 10],[1 1],'b','LineWidth',3)


subplot(2,2,3), hold on
plot(time(1:10:end),reshaped_fA(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel('Time since switch (s)')
ylabel('ORN Response (Hz)')
set(gca,'XLim',[0 10],'YLim',[0 85])
plot([1 5],[80 80],'r','LineWidth',3)
plot([6 10],[80 80],'b','LineWidth',3)

% first show the high contrast epochs
ax(1) = subplot(2,2,2); hold on
ax(2) = subplot(2,2,4); hold on
xlabel(ax(2),'Projected Stimulus (V)')
ylabel(ax(2),'Normalised Response')
ylabel(ax(1),'Probability')

% high variance
temp = fA_pred(1e3:5e3,:); 
x = 0:.05:2;
y = NaN(length(x)-1,width(temp)); 
cy = y;
for i = 1:width(temp)
	y(:,i) = histcounts(temp(:,i),x);
	y(:,i) = y(:,i)/sum(y(:,i));
	cy(:,i) = cumsum(y(:,i));
end
l = errorShade(ax(1),x(2:end),mean(y,2),sem(y'),'Color','r');
set(l(1),'LineWidth',2)

% plot integral -- theoretical prediction for best coding
plot(ax(2),x(2:end),mean(cy,2),'r--');

% now plot the actual i/o curve
x = fA_pred(1e3:5e3,:);
y = reshaped_fA(1e3:5e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; cy(:,rm_this)= [];
all_x = -1:.05:2;
all_y = NaN(length(all_x),width(y));
for i = 1:width(x)
	[~,data] = plotPieceWiseLinear(x(:,i),y(:,i),'nbins',40,'make_plot',false);
	data.y = data.y/max(data.y);
	all_y(:,i) = interp1(data.x,data.y,all_x);
end

% plot data
[line_handle2, shade_handle2] = errorShade(ax(2),all_x,nanmean(all_y,2),sem(all_y'),'Color',[1 0 0],'Shading',s);
uistack(shade_handle2,'bottom')
set(ax(1),'XLim',[0 2],'YLim',[0 .16])
set(ax(2),'XLim',[0 2],'YLim',[0 1.1])

% fake some plots for a nice legend
clear l
l(1) = plot(ax(2),NaN,NaN,'k');
l(2) = plot(ax(2),NaN,NaN,'k--');
L = {'ORN Response','Prediction'};
legend(l,L,'Location','southeast')


% now do low variance case
temp = fA_pred(6e3:9e3,:); 
x = 0:.05:2;
y = NaN(length(x)-1,width(temp)); 
cy = y;
for i = 1:width(temp)
	y(:,i) = histcounts(temp(:,i),x);
	y(:,i) = y(:,i)/sum(y(:,i));
	cy(:,i) = cumsum(y(:,i));
end
errorShade(ax(1),x(2:end),mean(y,2),sem(y'),'Color','b','Shading',s);

% plot integral -- theoretical prediction for best coding
plot(ax(2),x(2:end),mean(cy,2),'--','Color',[0 0 1]);

% now plot the actual i/o curve
x = fA_pred(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; cy(:,rm_this)= [];
all_x = -1:.05:2;
all_y = NaN(length(all_x),width(y));
for i = 1:width(x)
	[~,data] = plotPieceWiseLinear(x(:,i),y(:,i),'nbins',40,'make_plot',false);
	data.y = data.y/max(data.y);
	all_y(:,i) = interp1(data.x,data.y,all_x);
end
% plot data
[line_handle2, shade_handle2] = errorShade(ax(2),all_x,nanmean(all_y,2),sem(all_y'),'Color',[0 0 1],'Shading',s);
uistack(shade_handle2,'bottom')
set(ax(1),'XLim',[0 2],'YLim',[0 .16])
set(ax(2),'XLim',[0 2],'YLim',[0 1.1])

% also plot the distributions of the means 
h = axes(); hold(h,'on')
set(h,'Position',[.8 .8 .1 .13])
mean_1 = mean(reshaped_PID(1e3:5e3,:));
mean_2 = mean(reshaped_PID(6e3:end,:));
std_1 = std(reshaped_PID(1e3:5e3,:));
std_2 = std(reshaped_PID(6e3:end,:));
plot(h,mean_1,std_1,'r.')
plot(h,mean_2,std_2,'b.')
set(h,'XLim',[0 0.6],'YLim',[0 0.2])
xlabel('\mu (V)')
ylabel('\sigma (V)')


prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end



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


% calculate instantaneous gain everywhere 
[mean_stim,inst_gain,e] = makeFig6G(reshaped_PID(:),reshaped_fA(:),fA_pred(:),500);

% reshape back
inst_gain = reshape(inst_gain,size(reshaped_PID,1),size(reshaped_PID,2));
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


% now do the gain filter analysis 
history_lengths = logspace(-2,0.5,20);
history_lengths(14) = .5;
gain_K1_rho = NaN*history_lengths;
gain_K2_rho = NaN*history_lengths;

stim = reshaped_PID(:);
y = inst_gain(:);

for i = 1:length(history_lengths)
	% first do the simple box filter
	temp = floor(history_lengths(i)*1e3);
	temp = abs(filter(ones(temp,1),temp,stim));
	subplot(3,2,3), hold on
	gain_K1_rho(i) = spear(temp(1:100:end),y(1:100:end));
	if i == 14
		plot(temp(1:100:end),y(1:100:end),'k.')
		ylabel('Inst. Gain (Hz/V)')
		xlabel(['Stimulus projected ' char(10) 'using integrating filter'])
		legend(['\rho = ' oval(gain_K1_rho(i))])
	end
	


	% now the differentiating filter
	temp = floor(history_lengths(i)*1e3/2);
	temp = abs(filter([ones(temp,1); -ones(temp,1)],2*temp,stim-nanmean(stim)));
	gain_K2_rho(i) = spear(temp(1:100:end),y(1:100:end));
	subplot(3,2,4), hold on
	if i == 14
		plot(temp(1:100:end),y(1:100:end),'k.')
		xlabel(['Stimulus projected ' char(10) 'using differentiating filter'])
		legend(['\rho = ' oval(gain_K2_rho(i))])
	end
	

end


subplot(3,2,5), hold on
xlabel('History Length (s)')
ylabel('\rho')
plot(history_lengths,gain_K1_rho,'Color','k','Marker','d')
set(gca,'YLim',[-.4 .2])
subplot(3,2,6), hold on
plot(history_lengths,gain_K2_rho,'Color',[0 .6 0],'Marker','+')
set(gca,'YLim',[-.4 .2])
xlabel('History Length (s)')
ylabel('\rho')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% 

pFooter;
