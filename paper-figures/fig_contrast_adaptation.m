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

X = fA_pred;
Y = reshaped_fA;
S = reshaped_PID;
S(1:1e3,:) = NaN;
X(1:1e3,:) = NaN;
Y(1:1e3,:) = NaN;
if ~exist('sc','var')
	[gain,gain_r2,sc] = findEnsembleGain(X,Y,S,'step_size',10);
end

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on

% show how gain changes with time
[ax,h1,h2] = plotyy(time-5,gain,time-5,sc);
set(h1,'Marker','+','LineStyle','none')
set(h2,'Marker','.','LineStyle','none')
xlabel(ax(1),'Time since switch (s)')
ylabel(ax(1),'Inst. Gain (Hz/V)')
ylabel(ax(2),'Stimulus contrast')

% put a timescale on this change by finding the time to half asymptote 
a = nanmean(gain(1:4e3));
z = nanmean(gain(6e3:end));
tau_fA = 5e3 + find(gain(5e3:end) > a + (z-a)/2,1,'first');

a = nanmean(sc(1:4e3));
z = nanmean(sc(6e3:end));
tau_sc = find(sc < z+ (a-z)/2,1,'first');


prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end


%% Zero-parameter fits to data
% In this section we attempt to use a divisive gain control model that is sensitive to the contrast of the stimulus in some recent window to see if we can explain this change in gain. The model is as follows:
%
% $$ \hat{R}(t)=\frac{\alpha x(t)}{1+\beta\left|K_{g}\otimes s(t)\right|} $$
% 

X = fA_pred;

%%
% In the following figure, we ignore the kinetics for now and calculate the degree of contrast-dependent gain control over the entire epoch to find the $\alpha$ and $\beta$ parameters. Specifically, we solve the following equations: 
% 
% $$ gain_{lo}=\frac{\alpha}{1+\beta\sigma_{lo}} $$
% 
% and 
% 
% $$ gain_{hi}=\frac{\alpha}{1+\beta\sigma_{hi}} $$
%
% for $\alpha$ and $\beta$
% 

%%
% The following figure shows the ORN response vs. the linear projections for the low (blue) and high (variance) contrast epochs. We then plot the response vs. the gain-corrected linear projections, following the formula above. We also plot the gain-corrected linear projections as a function of the linear projections, to show that there is a gain change similar to the first plot. Finally, we plot the $r^2$ of the linear projection and the gain-corrected linear projection vs. the data, to see if adding the gain-correcting term explains more of the variance in the data. 

% measure the gain in each trial in the high and lo variance epoch
temp = NaN(width(X),1);
temp2 = NaN(width(X),1);
for i = 1:length(temp)
	x = X(1e3:5e3,i); y = reshaped_fA(1e3:5e3,i);
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = []; y(rm_this) = [];
	if length(x) > 1e3
		ff = fit(x,y,'poly1');
		temp(i) = ff.p1;
	end
	x = X(6e3:end,i); y = reshaped_fA(6e3:end,i);
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = []; y(rm_this) = [];
	if length(x) > 1e3
		ff = fit(x,y,'poly1');
		temp2(i) = ff.p1;
	end
end

g_lo = nanmean(temp2);
g_hi = nanmean(temp);

s_lo = nanmean(sc(6e3:end));
s_hi = nanmean(sc(1e3:5e3));

% now calculate beta
B = ((g_lo/g_hi) - 1)/(s_hi - s_lo*(g_lo/g_hi));

% now calculate alpha
A = g_hi*(1+B*s_hi);

% fix the gain
XG = X;
for i = 1:width(X)
	XG(1e3:5e3,i) = A*X(1e3:5e3,i)./(1 + B*s_hi);
	XG(6e3:end,i) = A*X(6e3:end,i)./(1 + B*s_lo);
end


figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plotPieceWiseLinear(X(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(X(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('K_{lo} \otimes s')
ylabel('Response (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(XG(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(XG(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('$\hat{R}$','interpreter','latex')
ylabel('Response (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(X(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(X(6e3:end,:),XG(6e3:end,:),'nbins',50,'Color','b');
ylabel('$\hat{R}$','interpreter','latex')
xlabel('K_{lo} \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(X),1);
r2_XG = r2_X;
for i = 1:width(X)
	fp = X([1e3:5e3 6e3:10e3],i); r = reshaped_fA([1e3:5e3 6e3:10e3],i);
	try
		r2_X(i) = rsquare(fp,r);
	catch
	end
	fp = XG([1e3:5e3 6e3:10e3],i);
	try
		r2_XG(i) = rsquare(fp,r);
	catch
	end
end
plot(r2_X,r2_XG,'k+')
plot([0 1],[0 1],'k--')
xlabel('r^2 (Linear Model)')
ylabel('r^2 (Gain-corrected)')
suptitle('Gain corrected by static contrast')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end


%%
% So the divisive term "works" in the sense that it has now the same gain in the high and low contrast case, but it ends up moving the curves relative to each other, so they no longer overlap. Also, the gain correction reduces the overall $r^2$, which is not good. 

%%
% Perhaps the $r^2$ is getting hammered because of this shift. What if we compute the $r^2$ separately in the high and low variance epochs? 

r2_X = NaN(width(X),2);
r2_XG = r2_X;
for i = 1:width(X)
	fp = X(1e3:5e3,i); r = reshaped_fA(1e3:5e3,i);
	try
		r2_X(i,1) = rsquare(fp,r);
	catch
	end
	fp = X(6e3:end,i); r = reshaped_fA(6e3:end,i);
	try
		r2_X(i,2) = rsquare(fp,r);
	catch
	end

	fp = XG(1e3:5e3,i); r = reshaped_fA(1e3:5e3,i);
	try
		r2_XG(i,1) = rsquare(fp,r);
	catch
	end
	fp = XG(6e3:end,i); r = reshaped_fA(6e3:end,i);
	try
		r2_XG(i,2) = rsquare(fp,r);
	catch
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(r2_X(:,1),r2_XG(:,1),'r+')
plot([0 1],[0 1],'k--')
xlabel('r^2 (Linear Model)')
ylabel('r^2 (Gain-corrected)')
suptitle('Gain corrected by static contrast')
title('High contrast')

subplot(1,2,2), hold on
plot(r2_X(:,2),r2_XG(:,2),'b+')
plot([0 1],[0 1],'k--')
xlabel('r^2 (Linear Model)')
ylabel('r^2 (Gain-corrected)')
suptitle('Gain corrected by static contrast')
title('Low contrast')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we take the kinetics into account. We assume that the timescale of the gain filter is around 100ms, as specified by the plot earlier, and that the filter has a simple differentiating shape. 

% filter the stimulus using a differentiating filter 92ms long
Kg = [ones(round((tau_fA-tau_sc)/2),1); -ones(round((tau_fA-tau_sc)/2),1)];
Shat = reshaped_PID;
for i = 1:width(reshaped_PID)
	Shat(:,i) = abs(filter(Kg,sum(abs(Kg)),reshaped_PID(:,i)));
end

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1:2); hold on
plot(time-5,nanmean(Shat,2),'k')
xlabel('Time since switch (s)')
ylabel('$| K_{g} \otimes s |$','interpreter','latex')
set(gca,'XLim',[-4 5])
set(gca,'YLim',[0 0.1])

subplot(1,3,3), hold on
plot(sc,nanmean(Shat,2),'k+')
xlabel('Stimulus contrast')
ylabel('$| K_{g} \otimes s |$','interpreter','latex')
set(gca,'YLim',[0 0.1])
prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%%
% We now correct the linear prediction using a differentiating filter operating on the stimulus. 

XG = X;
for i = 1:width(X)
	XG(:,i) = A*X(:,i)./(1 + B*Shat(:,i));
end


figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plotPieceWiseLinear(X(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(X(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('K_{lo} \otimes s')
ylabel('Response (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(XG(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(XG(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('$\hat{R}$','interpreter','latex')
ylabel('Response (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(X(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(X(6e3:end,:),XG(6e3:end,:),'nbins',50,'Color','b');
ylabel('$\hat{R}$','interpreter','latex')
xlabel('K_{lo} \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(X),1);
r2_XG = r2_X;
for i = 1:width(X)
	fp = X([1e3:5e3 6e3:10e3],i); r = reshaped_fA([1e3:5e3 6e3:10e3],i);
	try
		r2_X(i) = rsquare(fp,r);
	catch
	end
	fp = XG([1e3:5e3 6e3:10e3],i);
	try
		r2_XG(i) = rsquare(fp,r);
	catch
	end
end
plot(r2_X,r2_XG,'k+')
plot([0 1],[0 1],'k--')
xlabel('r^2 (Linear Model)')
ylabel('r^2 (Gain-corrected)')
suptitle('Gain corrected by contrast filter ')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%%
% Perhaps the $r^2$ is getting hammered because of this shift. What if we compute the $r^2$ separately in the high and low variance epochs? 

r2_X = NaN(width(X),2);
r2_XG = r2_X;
for i = 1:width(X)
	fp = X(1e3:5e3,i); r = reshaped_fA(1e3:5e3,i);
	try
		r2_X(i,1) = rsquare(fp,r);
	catch
	end
	fp = X(6e3:end,i); r = reshaped_fA(6e3:end,i);
	try
		r2_X(i,2) = rsquare(fp,r);
	catch
	end

	fp = XG(1e3:5e3,i); r = reshaped_fA(1e3:5e3,i);
	try
		r2_XG(i,1) = rsquare(fp,r);
	catch
	end
	fp = XG(6e3:end,i); r = reshaped_fA(6e3:end,i);
	try
		r2_XG(i,2) = rsquare(fp,r);
	catch
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(r2_X(:,1),r2_XG(:,1),'r+')
plot([0 1],[0 1],'k--')
xlabel('r^2 (Linear Model)')
ylabel('r^2 (Gain-corrected)')
suptitle('Gain corrected by contrast filter')
title('High contrast')

subplot(1,2,2), hold on
plot(r2_X(:,2),r2_XG(:,2),'b+')
plot([0 1],[0 1],'k--')
xlabel('r^2 (Linear Model)')
ylabel('r^2 (Gain-corrected)')
title('Low contrast')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%%
% It still looks bad in each epoch, suggesting that the gain-correcting term isn't correct. Either the timescale of the constants are wrong. 

%%
% But there's a deeper problem here. When we add on the gain-sensitive term, and determine the constants as we did above, the only thing we require is that the gain in the two places is the same. This is a very tolerant restriction, and as we saw, we satisfied it without satisfying a more important one: that the two curves must intersect. Now, we impose a stricter condition: that each point on the input-output curve is identical in the low and high contrast case. This means that we require:
%
% $$ R_{i}=f\left(\frac{\alpha X_{i}^{lo}}{1+\beta\sigma_{lo}}\right)=f\left(\frac{\alpha X_{i}^{hi}}{1+\beta\sigma_{hi}}\right) $$
% 

%%
% Here we ask if it possible for a single $\beta$ to ever give us a condition where the two curves overlap. We can solve for $\beta$ using this stricter condition to get:
% 
% $$ \beta_{i}=\frac{X{}_{i}^{hi}-X{}_{i}^{lo}}{X{}_{i}^{lo}\sigma_{hi}-X{}_{i}^{lo}\sigma_{lo}} $$
% 

% match curves
data_hi.x = interp1(data_hi.y,data_hi.x,data_lo.y);
data_hi.y = data_lo.y;

B = NaN*data_lo.y;
for i = 1:length(data_lo.y)
	B(i) = data_hi.x(i) - data_lo.x(i);
	B(i) = B(i)/(data_lo.x(i)*s_hi - s_lo*data_hi.x(i));
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data_lo.x,B,'k+')
xlabel('K \otimes s')
ylabel('\beta')

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks like no single $\beta$ can make these curves lie on top of each other. 

%%
% To verify our intuition, we now generate synthetic data using the sontrast-sensitive gain-corrected model, and see if we can ever get the two input-output curves to cross.

clear p
p.   s0 = -0.3204;
p.tau_z = 100.1875;
p.  n_y = 2;
p.tau_y = 20.5000;
p.    A = 129.5954;
p.    B = 0;

all_B = linspace(0,6,4);

S = reshaped_PID(:,50:100);
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
	subplot(2,2,i); hold on
	p.B = all_B(i);
	R = NaN*S;
	for j = 1:width(S_lo)
		R(:,j) = DAModel_contrast(S(:,j),p);
	end
	R(R<0) = 0;
	KDA = fitFilter2Data(vectorise(S(1e3:5e3,:)),vectorise(R(1e3:5e3,:)),'reg',1,'filter_length',500);
	% make linear predictions and plot input-response curves
	temp = NaN*S;
	for j = 1:width(S)
		temp(:,j) = filter(KDA,1,S(:,j));
	end
	plotPieceWiseLinear(temp(1e3:5e3,:),R(1e3:5e3,:),'Color','r','nbins',50);
	plotPieceWiseLinear(temp(6e3:end,:),R(6e3:end,:),'Color','b','nbins',50);
	title(['\beta = ' oval(p.B)])
	xlabel('Linear Prediction (a.u.)')
	ylabel('DA-like contrast model (a.u.)')
end

prettyFig('FixLogX=1;','fs=16;')

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
% 

pFooter;
