% fig_contrast_adaptation2.m
% makes figure: contrast adaptation in ORNs
% 
% created by Srinivas Gorur-Shandilya at 7:10 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

dm = dataManager;

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


[PID, LFP, fA, paradigm, orn] = consolidateData(dm.getPath('7955d1ed77512dfe3452b39d71a50e1b'),1);

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

% compute gains per trial
lo_gain = NaN(width(reshaped_PID),1);
hi_gain = NaN(width(reshaped_PID),1);
for i = 1:width(reshaped_PID)
	y = reshaped_fA(1e3:4e3,i);
	x = fA_pred(1e3:4e3,i);
	r = max(y) - min(y);
	l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
	try
		ff = fit(x(:),y(:),'poly1');
		hi_gain(i) = ff.p1;
	catch
	end

	y = reshaped_fA(6e3:9e3,i);
	x = fA_pred(6e3:9e3,i);
	r = max(y) - min(y);
	l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
	try
		ff = fit(x(:),y(:),'poly1');
		lo_gain(i) = ff.p1;
	catch
	end
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

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
clear ax
ax(1) = subplot(3,10,1:7);
ax(2) = subplot(3,10,8:10);
ax(3) = subplot(3,10,11:17);
ax(4) = subplot(3,10,18:20);

ax(5) = subplot(3,3,7); 
ax(6) = subplot(3,3,8);
ax(7) = subplot(3,3,9);
for i = 1:length(ax)
	hold(ax(i),'on');
end  

plot(ax(1),time(1:10:end),reshaped_PID(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel(ax(1),'Time since switch (s)')
ylabel(ax(1),'Stimulus (V)')
set(ax(1),'XLim',[0 10],'YLim',[0 1.1])
plot(ax(1),[1 5],[1 1],'r','LineWidth',3)
plot(ax(1),[6 10],[1 1],'b','LineWidth',3)

% plot the distributions of the projected stimulus
x = fA_pred(1e3:10:5e3,:); x = x(:);
x = x - mean(x);
x = x/std(x);
x = x*mean(std(reshaped_PID(1e3:5e3,:)));
x = x+mean(mean(reshaped_PID(1e3:5e3,:)));
hx = linspace(min(x),max(x),100);
hxx = hx(1:end-1) + mean(diff(hx));
hy = histcounts(x,hx);
hy = hy/sum(hy);
plot(ax(2),hy,hxx,'r')

x = fA_pred(6e3:10:9e3,:); x = x(:);
x = x - mean(x);
x = x/std(x);
x = x*mean(std(reshaped_PID(6e3:9e3,:)));
x = x+mean(mean(reshaped_PID(6e3:9e3,:)));
hxx = hx(1:end-1) + mean(diff(hx));
hy = histcounts(x,hx);
hy = hy/sum(hy);
plot(ax(2),hy,hxx,'b')
xlabel(ax(2),'Probability')


plot(ax(3),time(1:10:end),reshaped_fA(1:10:end,1:10:end),'Color',[.5 .5 .5 .5]);
xlabel(ax(3),'Time since switch (s)')
ylabel(ax(3),'ORN Response (Hz)')
set(ax(3),'XLim',[0 10],'YLim',[0 85])
plot(ax(3),[1 5],[80 80],'r','LineWidth',3)
plot(ax(3),[6 10],[80 80],'b','LineWidth',3)

% and show the response distributions 
x = reshaped_fA(1e3:10:5e3,:); x = x(:);
hx = linspace(min(x),max(x),100);
hxx = hx(1:end-1) + mean(diff(hx));
hy = histcounts(x,hx);
hy = hy/sum(hy);
plot(ax(4),hy,hxx,'r')
xlabel(ax(4),'Probability')

x = reshaped_fA(6e3:10:9e3,:); x = x(:);
hxx = hx(1:end-1) + mean(diff(hx));
hy = histcounts(x,hx);
hy = hy/sum(hy);
plot(ax(4),hy,hxx,'b')
xlabel(ax(4),'Probability')

% integrate stimulus distributions and re-plot as a prediction
temp = fA_pred(1e3:5e3,:); temp = nonnans(temp(:));
x = min(min(fA_pred)):0.02:max(max(fA_pred));
y = histcounts(temp,x);
y = y/sum(y); 
plot(ax(5),x(2:end),cumsum(y),'r--')

temp = fA_pred(6e3:end,:); temp = nonnans(temp(:));
y = histcounts(temp,x);
y = y/sum(y);
plot(ax(5),x(2:end),cumsum(y),'b--')

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
plot(ax(5),data_hi.x,data_hi.y,'r','LineWidth',2)
plot(ax(5),data_lo.x,data_lo.y,'b','LineWidth',2)

% show gain as function of contrast
lo_gain(lo_gain==0) = NaN;
hi_gain(hi_gain==0) = NaN;
x = std(reshaped_PID(1e3:4e3,:));
plot(ax(6),x,hi_gain,'r+')
x = std(reshaped_PID(6e3:9e3,:));
plot(ax(6),x,lo_gain,'b+')
xlabel(ax(6),'\sigma_{Stimulus} (V)')
ylabel(ax(6),'ORN Gain (Hz/V)')

% create  some phantom plots for a nice legend
clear l
l(1) = plot(ax(5),NaN,NaN,'k--');
l(2) = plot(ax(5),NaN,NaN,'k');
lh = legend(l,{'c.d.f.','Measured'},'Location','southeast');
lh.FontSize = 10; lh.Box = 'off';

% now use Laughlin to compute the gains per trial
x = linspace(min(min(reshaped_PID)),max(max(reshaped_PID)),100);
laughlin_hi_gain = NaN*hi_gain;
laughlin_lo_gain = NaN*lo_gain;
for i = 1:width(reshaped_PID)
	s = reshaped_PID(1e3:5e3,i);
	hy = histcounts(s,x);
	hy = cumsum(hy);
	hy = hy/hy(end);
	yy = hy(find((hy>.4),1,'first'):find((hy>.6),1,'first'));
	xx = hx(find((hy>.4),1,'first'):find((hy>.6),1,'first'));
	ff = fit(xx(:),yy(:),'poly1');
	laughlin_hi_gain(i) = ff.p1;

	s = reshaped_PID(6e3:9e3,i);
	hy = histcounts(s,x);
	hy = cumsum(hy);
	hy = hy/hy(end);
	yy = hy(find((hy>.4),1,'first'):find((hy>.6),1,'first'));
	xx = hx(find((hy>.4),1,'first'):find((hy>.6),1,'first'));
	ff = fit(xx(:),yy(:),'poly1');
	laughlin_lo_gain(i) = ff.p1;
end


plot(ax(7),laughlin_lo_gain,lo_gain,'b+')
plot(ax(7),laughlin_hi_gain,hi_gain,'r+')
r2 = rsquare([laughlin_lo_gain; laughlin_hi_gain],[lo_gain; hi_gain]);
ylabel(ax(7),'ORN Gain (Hz/V)')
xlabel(ax(7),'c.d.f slope (a.u.)')
l = plot(ax(7),NaN,NaN,'k+');
lh = legend(l,['r^2 = ' oval(r2)]);
lh.Location = 'southeast';

% cosmetics
ax(1).Position(3) = .53;
ax(3).Position(3) = .53;
set(ax(2),'YTick',[],'YLim',ax(1).YLim)
set(ax(4),'YTick',[],'YLim',ax(3).YLim)

xlabel(ax(5),'Projected Stimulus (V)')
ylabel(ax(5),'Firing Rate (norm)')

prettyFig('fs',14,'lw',1.5)
labelFigure

if being_published
	snapnow
	delete(gcf)
end

%% Supplementary Figure

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
x = std(reshaped_PID(1e3:4e3,:));
y = std(reshaped_fA(1e3:4e3,:))./x; y(y==0) = NaN;
plot(x,y,'r+')
x = std(reshaped_PID(6e3:9e3,:));
y = std(reshaped_fA(6e3:9e3,:))./x; y (y==0) = NaN;
plot(x,y,'b+')
ylabel('\sigma_{Firing Rate}/\sigma_{Stimulus} (Hz/V)')
xlabel(gca,'\sigma_{Stimulus} (V)')

subplot(1,2,2), hold on
x = std(reshaped_PID(1e3:4e3,:));
plot(x,laughlin_hi_gain,'r+')
x = std(reshaped_PID(6e3:9e3,:));
plot(x,laughlin_lo_gain,'b+')
ylabel(gca,'c.d.f Slope (a.u.)')
xlabel(gca,'\sigma_{Stimulus} (V)')

prettyFig('fs',18,'lw',2)
labelFigure

if being_published
	snapnow
	delete(gcf)
end

return

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


prettyFig('FixLogX',true,'fs',16)

if being_published
	snapnow
	delete(gcf)
end

%% Zero-parameter fits to data
% In this section, we attempt to directly measure the parameters of a "zero-parameter" model that includes gain control terms that are sensitive to the mean and the contrast of the stimulus. The model response is given by
% 
% $$ \hat{R}(t)=f\left(g(t)K_{r}\otimes s(t)\right) $$
%
% where 
% 
% $$ g(t)=\frac{1}{1+\beta_{\mu}K_{\mu}\otimes s(t)} $$
%
% is the time-dependent gain of the transduction currents and the nonlinear function is a Hill function. The exponent of the Hill function is also controlled by the stimulus as follows:
%
% $$ n(t)=\frac{n_{0}}{1+\beta_{\sigma}K_{\sigma}\otimes\left[s'(t)\right]_{+}} $$ 
%  

%%
% First, we concentrate on determining the parameters controlling the contrast-sensitive gain change. We fit Hill functions to the measured input-output curves:

[~,data_hi] = plotPieceWiseLinear(fA_pred(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'make_plot',false);
[~,data_lo] = plotPieceWiseLinear(fA_pred(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'make_plot',false);
data_lo.y = data_lo.y/nanmean(max(reshaped_fA));
data_hi.y = data_hi.y/nanmean(max(reshaped_fA));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
ft = fittype('hill2(x,k,n,x_offset)');
ff = fit(data_lo.x(:),data_lo.y(:),ft,'StartPoint',[.5 2 0],'Lower',[0 1 0],'Upper',[10 10 0],'MaxIter',1e4);
plot(data_lo.x,data_lo.y,'k')
l = plot(data_lo.x,ff(data_lo.x),'r');
legend(l,['K = ' oval(ff.k) ', n = ' oval(ff.n)],'Location','southeast')
title('Low contrast')
xlabel('Proj. Stimulus')
ylabel('Response (norm)')
n_lo = ff.n;

subplot(1,2,2); hold on
ft = fittype('hill2(x,k,n,x_offset)');
ff = fit(data_hi.x(:),data_hi.y(:),ft,'StartPoint',[.5 2 0],'Lower',[0 1 0],'Upper',[10 10 0],'MaxIter',1e4);
plot(data_hi.x,data_hi.y,'k')
l = plot(data_hi.x,ff(data_hi.x),'r');
legend(l,['K = ' oval(ff.k) ', n = ' oval(ff.n)],'Location','southeast')
title('High contrast')
xlabel('Proj. Stimulus')
n_hi = ff.n;

prettyFig('FixLogX',true,'fs',16,'EqualiseX',true,'EqualiseY',true)

if being_published
	snapnow
	delete(gcf)
end

%%
% In the following figure, we ignore the kinetics for now and calculate the degree of contrast-dependent gain control over the entire epoch to find the $\n_{0}$ and $\beta$ parameters. $\beta$ is given by
%
% $$ \beta=\frac{n_{hi}-n_{lo}}{n_{lo}\sigma_{lo}-n_{hi}\sigma_{hi}} $$
% 

% compute the derivative everywhere
Sd = reshaped_PID;
for i = 1:width(Sd)
	Sd(:,i) = (filtfilt(ones(10,1),10,[0; diff(Sd(:,i))]));
end
Sd(Sd<0) = 0;


s_hi = nanmean(nanmean(Sd(1e3:5e3,:)));
s_lo = nanmean(nanmean(Sd(6e3:end,:)));

B = (n_hi-n_lo)/(n_lo*s_lo - n_hi*s_hi);
n0 = n_lo*(1+B*s_lo);

%%
% The calculated B and n0 are:

B ,n0

% correct the linear prediction
XG = fA_pred;
for i = 1:width(XG)
	n = n0/(1 + B*s_hi);
	XG(1e3:5e3,i) = hill2(fA_pred(1e3:5e3,i),.79,n,0);
	n = n0/(1 + B*s_lo);
	XG(6e3:end,i) = hill2(fA_pred(6e3:end,i),.79,n,0);
end

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1); hold on
plotPieceWiseLinear(fA_pred(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(fA_pred(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('K \otimes s')
ylabel('Response (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(XG(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(XG(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('$\hat{R}$','interpreter','latex')
ylabel('Response (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(fA_pred(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(fA_pred(6e3:end,:),XG(6e3:end,:),'nbins',50,'Color','b');
ylabel('$\hat{R}$','interpreter','latex')
xlabel('K \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(fA_pred),1);
r2_XG = r2_X;
for i = 1:width(fA_pred)
	fp = fA_pred([1e3:5e3 6e3:10e3],i); r = reshaped_fA([1e3:5e3 6e3:10e3],i);
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

% convert into remaining variance accounted for
[y,x] = histcounts((r2_XG-r2_X)./(1-r2_X),-1:.02:1);
y = y/sum(y);
y = y*length(y);
x = x(1:end-1) + mean(diff(x));
plot(x,y,'k')
plot([0 0],[0 100],'k--')
set(gca,'XLim',[-1 1],'YLim',[0 10])
xlabel(['Additional Variance' char(10) 'accounted for'])
ylabel('pdf')

prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end


%%
% Now we take the kinetics into account. To determine the timescale of contrast gain control, we fit a contrast-sensitive model to the data, keeping the B and n0 parameters fixed at what we measured ealrier. We now correct the linear prediction using a filter operating on the derivative of the stimulus. The parameters of the best-fit contrast sensitive model are:


clear d
i = 8;
ft = -99:700;
fp = convolve(1e-3*(1:length(PID)),PID(:,i),K2_mean,ft);
S = filtfilt(ones(10,1),10,[0; diff(PID(:,i))]);
d.stimulus = [fp(1e4:end-1e4), S(1e4:end-1e4)];
d.response = fA(1e4:end-1e4,i);
d.response(1:1e3) = NaN;
clear p
p. n0 = 8.1800;
p.tau = 215.9647;
p.  K = 0.7623;
p.  A = 65.2611;
p.  B = 679;
p.  n = 2.2969;

disp(p)

XG = X;
for i = 1:width(X)
	temp = [fA_pred(:,i), Sd(:,i)];
	XG(:,i) = contrastLNModel(temp,p);
end


figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plotPieceWiseLinear(X(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(X(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('K \otimes s')
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
xlabel('K \otimes s')

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
% convert into remaining variance accounted for
[y,x] = histcounts((r2_XG-r2_X)./(1-r2_X),-1:.02:1);
y = y/sum(y);
y = y*length(y);
x = x(1:end-1) + mean(diff(x));
plot(x,y,'k')
plot([0 0],[0 100],'k--')
set(gca,'XLim',[-1 1],'YLim',[0 10])
xlabel(['Additional Variance' char(10) 'accounted for'])
ylabel('pdf')


prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
% 

pFooter;
