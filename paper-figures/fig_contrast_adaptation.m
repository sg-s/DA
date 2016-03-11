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

temp = fA_pred(1e3:5e3,:); temp = nonnans(temp(:));
x = min(min(fA_pred)):0.02:max(max(fA_pred));
y = histcounts(temp,x);
y = y/sum(y); 
plot(ax(1),x(2:end),y,'r')

% integrate it and re-plot as a prediction
plot(ax(2),x(2:end),cumsum(y),'r--')

temp = fA_pred(6e3:end,:); temp = nonnans(temp(:));
y = histcounts(temp,x);
y = y/sum(y);
plot(ax(1),x(2:end),y,'b')

% integrate it and re-plot as a prediction
plot(ax(2),x(2:end),cumsum(y),'b--')


% now plot the actual i/o curve: high contrast
x = fA_pred(1e3:5e3,:);
y = reshaped_fA(1e3:5e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | max(y) < 40 | min(y) > 10);
x(:,rm_this) = []; y(:,rm_this) = []; 
[~,data_hi] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);

% low contrast
x = fA_pred(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | max(y) < 40 | min(y) > 10);
x(:,rm_this) = []; y(:,rm_this) = [];
[~,data_lo] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);

% normalise globally
m = min([data_lo.y data_hi.y]);
data_lo.y = data_lo.y - m; data_hi.y = data_hi.y - m;

M = max([data_lo.y data_hi.y]);
data_lo.y = data_lo.y/M; data_hi.y = data_hi.y/M;

% plot data
plot(ax(2),data_hi.x,data_hi.y,'r','LineWidth',2)
plot(ax(2),data_lo.x,data_lo.y,'b','LineWidth',2)

% create  some phantom plots for a nice legend
clear l
l(1) = plot(ax(2),NaN,NaN,'k--');
l(2) = plot(ax(2),NaN,NaN,'k');
legend(l,{'Prediction','Response'},'Location','southeast')

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

X = fA_pred(4000:6000,:);
Y = reshaped_fA(4000:6000,:);
S = reshaped_PID(4000:6000,:);
if ~exist('sc','var')
	[gain,gain_r2,sc] = findEnsembleGain(X,Y,S);
end

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on

% show how gain changes with time
ax = plotyy(time(4000:6000)-5,gain,time(4000:6000)-5,sc);
xlabel(ax(1),'Time since switch (s)')
ylabel(ax(1),'Inst. Gain (Hz/V)')
ylabel(ax(2),'Simulus contrast')

% put a timescale on this change by finding the time to half asymptote 
a = nanmean(gain(1:.5e3));
z = nanmean(gain(1.5e3:end));
tau_fA = find(gain > a+ (z-a)/2,1,'first');

a = nanmean(sc(1:.5e3));
z = nanmean(sc(1.5e3:end));
tau_sc = find(sc < z+ (a-z)/2,1,'first');


prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%% Zero-parameter fits to data
% In this section we attempt to use a divisive gain control model that is sensitive to the contrast of the stimulus in some recent window to see if we can explain this change in gain. First, we calculate the filter only in the low-variance epoch (because the gain here is high), and use it to project the stimulus everywhere. 

% also determine the degree of contrast-sensitive gain control (the beta parameter)

% load the low-variance filter
load('.cache/VSA_K2.mat','K2')
K = (nanmean(squeeze(K2(2,:,:)),2));

% now project all the stimulus using the high-variance filter
X = fA_pred;
for i = 1:width(X)
	X(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K,ft);
end

%%
% The following figure shows the response as a function of stimulus projected using the low-variance filter, during the two epochs. Since we have an estimate of the timescale of contrast gain control, the only thing that remains is to find the degree to which this happens, and multiply this term with the linear projection. We can directly calculate this term from the data, for each point along these two curves as follows.

%%
% Let $X_{lo}$ and $X_{hi}$ be two output nonlinearities that map the projected stimulus onto the response of the neuron, corresponding to the low and high variance stimulus. In general, these two curves overlap each other, and for each point on the Y axis, there are two corresponding points on the X axis, mapped through these two functions. Let these two points be called $x_{hi}$ and $x_{lo}$. 

%%
% Now, we define 
% 
% $$ \hat{x}=\frac{x}{1+\beta\left|K_{g}\otimes s\right|} $$
% 
% where $\hat{x}$ is the gain-corrected projection of the stimulus. For now, we neglect the kinetics of contrast gain control, and assume
%
% $$\beta\left|K_{g}\otimes s\right| := \beta\left|s-\bar{s}\right| $$
%
% which means that
%
% $$ \beta_{lo}=\beta\left|s-\bar{s}\right|_{lo} $$
% 

%%
% Now, we impose the condition that for a single ORN response, the gain-corrected projection is independent of variance condition, and depends only on the linear projection. Solving for $\beta$, we get:
%
% $$ \beta=\frac{x_{lo}-x_{hi}}{x_{hi}\left|s-\bar{s}\right|_{lo}-x_{lo}\left|s-\bar{s}\right|_{hi}} $$

%%
% Everything on the right hand-side is measurable from the data, so we can directly estimate $\beta$ at each point along the input-output curve. The following figure shows this:


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
[~,data_hi] = plotPieceWiseLinear(X(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
[~,data_lo] = plotPieceWiseLinear(X(5e3:end,:),reshaped_fA(5e3:end,:),'nbins',50,'Color','b');
xlabel('K_{high} \otimes s')
ylabel('Response (Hz)')

% match curves
data_hi.x = interp1(data_hi.y,data_hi.x,data_lo.y);
data_hi.y = data_lo.y;

% filter the stimulus using a differentiating filter 92ms long
Kg = [ones(round((tau_fA-tau_sc)/2),1); -ones(round((tau_fA-tau_sc)/2),1)];
Shat = reshaped_PID;
for i = 1:width(reshaped_PID)
	Shat(:,i) = abs(filter(Kg,sum(abs(Kg)),reshaped_PID(:,i)));
end

s_hi = nanmean(nanmean(Shat(1e3:5e3,:)));
s_lo = nanmean(nanmean(Shat(5e3:end,:)));

B = NaN*data_lo.y;
for i = 1:length(data_lo.y)
	B(i) = data_lo.x(i) - data_hi.x(i);
	B(i) = B(i)/(data_hi.x(i)*s_lo -s_hi*data_lo.x(i));
end

subplot(1,2,2), hold on
plot(data_lo.x,B,'k+')
xlabel('K_{high} \otimes s')
ylabel('\beta')

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.999978291937413;

% Fit model to data.
[fitresult, gof] = fit(data_lo.x(:), B(:), ft, opts );
plot(data_lo.x,fitresult(data_lo.x),'r')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%%
% Suprisingly, $\beta$ is not constant but looks like a function of the projected stimulus. Since we can directly measure this in the data, we modify our gain term as follows: 
%
% $$ \hat{x}=\frac{x}{1+\beta(K_{r}\otimes s)\left|K_{g}\otimes s\right|} $$
% 
% We now replot the data, using the gain-corrected model:
% 


% XG = X;

% % pass it through the model
% for i = 1:width(XG)
% 	textbar(i,width(XG))
% 	xg = XG(1e3:5e3,i);
% 	r = reshaped_fA(1e3:5e3,i);
% 	xg(r > max(data_lo.y)) = NaN;
% 	xg(r < min(data_lo.y)) = NaN;
% 	xg = xg./(1+fitresult(xg)*s_hi);
% 	XG(1e3:5e3,i) = xg;

% 	xg = XG(6e3:end,i);
% 	r = reshaped_fA(6e3:end,i);
% 	xg(r > max(data_lo.y)) = NaN;
% 	xg(r < min(data_lo.y)) = NaN;
% 	xg = xg./(1+fitresult(xg)*s_lo);
% 	XG(6e3:end,i) = xg;
% end


% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% plotPieceWiseLinear(XG(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
% plotPieceWiseLinear(XG(6e3:end-100,:),reshaped_fA(6e3:end-100,:),'nbins',50,'Color','b');
% xlabel('Gain-corrected prediction')
% ylabel('Response (Hz)')


% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% plot(data_lo.x,data_lo.y,'b')
% plot(data_hi.x,data_hi.y,'r')

% xlabel('K_{high} \otimes s')
% ylabel('Response (Hz)')

% temp1 = Shat(1e3:end,10:25); temp1 = temp1(:);
% temp2 = X(1e3:end,10:25); temp2 = temp2(:);
% temp3 = reshaped_fA(1e3:end,10:25); temp3 = temp3(:);
% rm_this = isnan(temp1) | isnan(temp2) |isnan(temp3);
% temp1(rm_this) = []; temp2(rm_this) = []; temp3(rm_this) = [];
% clear d
% d.stimulus = [temp1 temp2];
% d.response = temp3; d.response = d.response(:);

% fA_gain_corrected = fA_pred;
% for i = 1:width(fA_pred)
% 	fA_gain_corrected(:,i) = contrastGainBeta([reshaped_PID(:,i) X(:,i)],p);
% end

% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% plotPieceWiseLinear(fA_gain_corrected(1e3:5e3,1:10:end),reshaped_fA(1e3:5e3,1:10:end),'nbins',50,'Color','r');
% plotPieceWiseLinear(fA_gain_corrected(5e3:end,1:10:end),reshaped_fA(5e3:end,1:10:end),'nbins',50,'Color','b');
% xlabel('K_{high} \otimes s')
% ylabel('Response (Hz)')

% fit this one number

%% Version Info
% 

pFooter;
