% fig_contrast_adaptation.m
% makes figure: contrast adaptation in ORNs
% 
% created by Srinivas Gorur-Shandilya at 7:10 , 03 October 2015. Contact me at http://srinivas.gs/contact/
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

all_offsets = [1 3 6 8];
window_length = 2;

% compute the slopes -- firing rate 
if exist('.cache/VSA_fA_slopes.mat','file') == 2
	load('.cache/VSA_fA_slopes.mat','n')
else
	n = NaN(width(reshaped_PID),length(all_offsets));
	for i = 1:width(reshaped_PID)
		textbar(i,width(reshaped_fA))
		for j = 1:length(all_offsets)
			x = fA_pred(:,i);
			y = reshaped_fA(:,i);

			a = round(1e3*(all_offsets(j) - window_length/2));
			z = round(1e3*(all_offsets(j) + window_length/2));

			x = circshift(x,length(x) - z );
			y = circshift(y,length(y) - z );

			x = x(end-window_length*1e3:end);
			y = y(end-window_length*1e3:end);

			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			y = y(x > 0.33*max(x) & x < .66*max(x));
			x = x(x > 0.33*max(x) & x < .66*max(x));

			try
				ff = fit(x(:),y(:),'poly1');
				n(i,j) = ff.p1;
			catch
			end
		end
	end
	save('.cache/VSA_fA_slopes.mat','n')
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

figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
subplot(2,4,1), hold on
plot(time(1:10:end),reshaped_PID(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel('Time since switch (s)')
ylabel('Stimulus (V)')

subplot(2,4,5), hold on
plot(time(1:10:end),reshaped_fA(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel('Time since switch (s)')
ylabel('ORN Response (Hz)')

% first show the high contrast epochs
ax(1) = subplot(2,4,2); hold on
ax(2) = subplot(2,4,6); hold on
xlabel(ax(2),'Projected Stimulus')
ylabel(ax(2),'Normalised Response')
ylabel(ax(1),'Stimulus Probability')

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
set(ax(1),'XLim',[0 2],'YLim',[-.01 .16])
set(ax(2),'XLim',[0 2],'YLim',[-.01 1.1])

% fake some plots for a nice legend
clear l
l(1) = plot(ax(2),NaN,NaN,'k');
l(2) = plot(ax(2),NaN,NaN,'k--');
L = {'ORN Response','Prediction'};
legend(l,L,'Location','southeast')


% now do low variance case
ax(1) = subplot(2,4,3); hold on
ax(2) = subplot(2,4,7); hold on
xlabel(ax(2),'Projected Stimulus')
ylabel(ax(2),'Normalised Response')
ylabel(ax(1),'Stimulus Probability')

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
set(ax(1),'XLim',[0 2],'YLim',[-.01 .16])
set(ax(2),'XLim',[0 2],'YLim',[-.01 1.1])



% show how gain changes with time
subplot(2,4,4), hold on
all_offsets = 0:0.1:10;
[lhl, shl]= errorShade(all_offsets,nanmean(n),nanstd(n)/sqrt(width(reshaped_fA)),'Color','k');
set(lhl,'LineWidth',2);
xlabel('Time since switch (s)')
ylabel('Inst. Gain (Hz/V)')


% now do the gain filter analysis 
history_lengths = logspace(-2,0.5,20);
gain_K1_rho = NaN*history_lengths;
gain_K2_rho = NaN*history_lengths;

stim = reshaped_PID(:);
global inst_gain
inst_gain = n';
inst_gain = inst_gain(1:100,:);
inst_gain = inst_gain(:);
inst_gain(inst_gain < 0) = NaN;

for i = 1:length(history_lengths)

	% first do the simple box filter
	temp = floor(history_lengths(i)*1e3);
	temp = filter(ones(temp,1),temp,stim);
	temp = (temp(1:100:end));

	gain_K1_rho(i) = spear(temp(1:10:end),inst_gain(1:10:end));

	% now the differentiating filter
	temp = floor(history_lengths(i)*1e3/2);
	temp = filter([ones(temp,1); -ones(temp,1)],2*temp,stim-nanmean(stim));
	temp = abs(temp(1:100:end));

	gain_K2_rho(i) = spear(temp(1:10:end),inst_gain(1:10:end));

end

subplot(2,4,8), hold on
plot(history_lengths,gain_K1_rho,'r')
plot(history_lengths,gain_K2_rho,'b')
legend({'Integrating','Differentiating'},'Location','southeast')
xlabel('History Length (s)')
ylabel('\rho')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end


%  ######## #### ##       ######## ######## ########   ######      ##  
% %  ##        ##  ##          ##    ##       ##     ## ##    ##    #### 
% %  ##        ##  ##          ##    ##       ##     ## ##           ##  
% %  ######    ##  ##          ##    ######   ########   ######          
% %  ##        ##  ##          ##    ##       ##   ##         ##     ##  
% %  ##        ##  ##          ##    ##       ##    ##  ##    ##    #### 
% %  ##       #### ########    ##    ######## ##     ##  ######      ##  

% %    ##       ######## ########  
% %    ##       ##       ##     ## 
% %    ##       ##       ##     ## 
% %    ##       ######   ########  
% %    ##       ##       ##        
% %    ##       ##       ##        
% %    ######## ##       ##        


% % make colour scheme for block analysis
% filter_length = 1000;
% offset = 200;

% % we are going to calculate only one filter/epoch
% sr = 1e3; % sampling rate, Hz
% if exist('.cache/VSA_K.mat','file') == 2
% 	load('.cache/VSA_K.mat','K')
% else
% 	filter_length = 1000;
% 	offset = 200;
% 	K = NaN(2,filter_length-offset,width(reshaped_LFP));
% 	for i = 1:width(reshaped_LFP)
% 		textbar(i,width(reshaped_PID))

% 		% calculate filter for large variance epoch
% 		stim = reshaped_PID(:,i);
% 		resp = reshaped_LFP(:,i);

% 		resp(1:1e3) = NaN;
% 		resp(5e3:end)= NaN;

% 		try
% 			[this_K,ft] = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
% 			K(1,:,i) = this_K(100:end-101);
% 		catch 
% 		end

% 		% calculate filter for low variance epoch
% 		stim = reshaped_PID(:,i);
% 		resp = reshaped_LFP(:,i);

% 		resp(1:6e3) = NaN;

% 		try
% 			[this_K,ft] = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
% 			K(2,:,i) = this_K(100:end-101);
% 		catch 
% 		end
% 	end
% 	mkdir('.cache')
% 	save('.cache/VSA_K.mat','K')
% end

% % make linear predictions on the de-trended data using a mean filter averaged over all cases
% K_mean = mean(squeeze(mean(K,1)),2);
% ft = -99:700;
% LFP_pred = NaN*reshaped_LFP;
% for i = 1:width(reshaped_LFP)
% 	LFP_pred(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K_mean,ft);
% end

% all_offsets = [1 3 6 8];
% window_length = 1;

% % compute the slopes
% if exist('.cache/VSA_LFP_slopes.mat','file') == 2
% 	load('.cache/VSA_LFP_slopes.mat','n')
% else
% 	n = NaN(width(reshaped_PID),length(all_offsets));
% 	for i = 1:width(reshaped_PID)
% 		textbar(i,width(reshaped_LFP))
% 		for j = 1:length(all_offsets)
% 			x = LFP_pred(:,i);
% 			y = reshaped_LFP(:,i);

% 			a = round(1e3*(all_offsets(j) - window_length/2));
% 			z = round(1e3*(all_offsets(j) + window_length/2));

% 			x = circshift(x,length(x) - z );
% 			y = circshift(y,length(y) - z );

% 			x = x(end-window_length*1e3:end);
% 			y = y(end-window_length*1e3:end);

% 			rm_this = isnan(x) | isnan(y);
% 			x(rm_this) = [];
% 			y(rm_this) = [];
% 			try
% 				ff = fit(x(:),y(:),'poly1');
% 				n(i,j) = ff.p1;
% 			catch
% 			end
% 		end
% 	end
% 	save('.cache/VSA_LFP_slopes.mat','n')
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

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
