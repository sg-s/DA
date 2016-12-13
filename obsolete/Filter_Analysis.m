%% Filter Estimation
% How do we best calculate the filters? A summary of different ways of estimating the filters and results from each. 


% clear vars
clearvars -except options

% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;



%% Overview of different regularisation methods. 
% Given some input stimulus _S_ and some output response _f_, the linear filter _K_ that best describes the data is given by 
%
% $$ \hat{C}*K=s'*f $$
% 
% where $\hat{C}$ is the regularised covariance matrix _C_, _s_ is a _N_ x _M_ matrix of the time-shifted stimulus and _f_ is the response vector. _N_ is the filter length, and _M_ is the length of the stimulus and response vectors - _N_. 
%
%%
% _C_ is the unscaled covariance matrix and is given by:
% 
% $$ C=s^{T}*s $$
% 
% can be regularised by different means:
% 
%% 
% 1) Carlotta regularises _C_ using: 
%% 
% $$ \hat{C}=C+\hat{r}I $$
%%
% 2) and Damon suggested:
% 
% 
% $$ \hat{C}=\frac{(C+\hat{r}I)*\mathrm{tr}(C)}{(\mathrm{tr}(C)+\hat{r}N)} $$ 
% 
% where _I_ is the identity matrix and _r_ is a free parameter called the regularisation factor that suppresses the high-frequency components of _K_. $\hat{r}$  is obtained by scaling _r_ by $\mu$,  the mean of the Eigenvalues of the covariance matrix _C_. _s_ is the stimulus vector (e.g. the PID) and _f_ is the response vector (here, the firing rate of the ORN). In practice, _K_ is estimated by a left matrix division:
%
% $$ K=C\setminus(s'*f) $$

%% Synthetic Data with No regularisation 
% Synthetic data is prepared using Gaussian random inputs and an exponential filter, and the output is the convolution of the input with the exponential filter, with 10% Gaussian Random noise. The stimulus is zero mean, as is the response, by construction. 
a = randn(1,10000);
filtertime = 0:3e-3:3e-3*333;
filter_length = 333;
Kexp = exp(-10*filtertime);
b = filter(Kexp,1,a) + 0.1*randn(1,10000); % adding noise to make the reconstruction harder

figure('outerposition',[0 0 800 450],'PaperUnits','points','PaperSize',[800 450]); hold on
subplot(2,1,1), hold on
plot(a)
title('Synthetic Stimulus')
subplot(2,1,2), hold on
plot(b)
title('Output of filter')

prettyFig;
snapnow;
delete(gcf);


%% Synthetic Data 2: effect of regularisation 
% We now create a new synthetic dataset, identical to the old one, except we filter the white noise input with some boxcar filter to remove high frequency components from the input. We find the best filters as before, and look at how regularisation parameter affects the choice of filter. 
a2 = filter(ones(1,10),1,a)/sqrt(10);
b2 = filter(Kexp,1,a2);

figure('outerposition',[0 0 800 450],'PaperUnits','points','PaperSize',[800 450]); hold on
subplot(2,1,1), hold on
plot(a2)
title('Synthetic Stimulus 2')
subplot(2,1,2), hold on
plot(b2)
title('Output of filter')

prettyFig;
snapnow;
delete(gcf);



if ~(exist('K1best_f') == 1)
	[K1best_f, diagnostics1f] = FindBestFilter(a2,b2,[],'regtype=1;');
end

if ~(exist('K2best_f') == 1)
	[K2best_f, diagnostics2f] = FindBestFilter(a2,b2,[],'regtype=2;');
end

%%
% The best filters are shown below for the two reg. methods:
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K1best_f,'r','LineWidth',2)
title('Reg method 1','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])
xlabel('Time')
ylabel('Filter amplitude')

subplot(1,2,2), hold on
plot(filtertime,K2best_f,'b','LineWidth',2)
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])
xlabel('Time')
ylabel('Filter amplitude')


snapnow;
delete(gcf);

%%
% the figure below shows how varying _r_ affects reg method 1:
PlotFilterDiagnostics2(diagnostics1f,marker_size,marker_size2,font_size,'Method 1')

snapnow;
delete(gcf);

%%
% the figure below shows how varying _r_ affects reg method 2:
PlotFilterDiagnostics2(diagnostics2f,marker_size,marker_size2,font_size,'method 2')

snapnow;
delete(gcf);

%    ########  ########    ###    ##          ########     ###    ########    ###          ##   
%    ##     ## ##         ## ##   ##          ##     ##   ## ##      ##      ## ##       ####   
%    ##     ## ##        ##   ##  ##          ##     ##  ##   ##     ##     ##   ##        ##   
%    ########  ######   ##     ## ##          ##     ## ##     ##    ##    ##     ##       ##   
%    ##   ##   ##       ######### ##          ##     ## #########    ##    #########       ##   
%    ##    ##  ##       ##     ## ##          ##     ## ##     ##    ##    ##     ##       ##   
%    ##     ## ######## ##     ## ########    ########  ##     ##    ##    ##     ##     ###### 

%% Real Data 1: 1-octen-3-ol flickering stimulus
% We use actual PID traces as the input, and the firing rate of a ORN as the output. The data looks like this:

% plot the data
load /local-data/DA-paper/data.mat

time = data(7).time;
PID = data(7).PID;
f = data(7).ORN;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(time,PID,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22])
ylabel('PID (a.u.)')

subplot(2,1,2), hold on
plot(time,f,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22])
ylabel('Firing rate (Hz)')

prettyFig;

snapnow;
delete(gcf);

%%
% The covariance of the input matrix _C_ looks like:
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% chop up the stimulus into blocks  
stim = PID- mean(PID);
OnlyThesePoints = 1:length(stim)-filter_length;
s = zeros(length(OnlyThesePoints), filter_length+1);
for i=OnlyThesePoints
    s(i,:) = stim(filter_length+i:-1:i);
end
C = s'*s;
imagesc(C), colorbar, axis tight
set(gca,'LineWidth',2,'FontSize',font_size)

snapnow;
delete(gcf);

%%
% Regularisation adds a large diagonal matrix to this, in effect making it look more like a diagonal matrix. Here, the mean of the eigenvalues of _C_ has been added to the diagonal. This is essentially what Method 1 does.  
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
r = trace(C)/length(C);
imagesc(C+r*eye(length(C)));
colorbar, axis tight
set(gca,'LineWidth',2,'FontSize',font_size)

snapnow;
delete(gcf);

%%
% The effect of regularisation using method 2 is shown on the correlation matrix below. The same _r_ is used, i.e., the mean of the eigenvalues of _C_. 
C2 = (C + r*eye(filter_length+1))*trace(C)/(trace(C) + r*filter_length);
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
imagesc(C2);
colorbar, axis tight
set(gca,'LineWidth',2,'FontSize',font_size)

snapnow;
delete(gcf);

%%
% Once again, we repeat the filter extraction. 
if ~(exist('K1real') == 1)
	[K1real,  diagnostics1r] = FindBestFilter(PID,f,[],'regtype=1;');
end


if ~(exist('K2real') == 1)
	[K2real,  diagnostics2r] = FindBestFilter(PID,f,[],'regtype=2;');
end


%%
% the figure below shows how varying _r_ affects reg method 1:
PlotFilterDiagnostics2(diagnostics1r,marker_size,marker_size2,font_size,'Method 1')

snapnow;
delete(gcf);

%%
% the figure below shows how varying _r_ affects reg method 2:
PlotFilterDiagnostics2(diagnostics2r,marker_size,marker_size2,font_size,'Method 2')

snapnow;
delete(gcf);


%%
% and the best filters look like: 

figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K1real,'r','LineWidth',2)
title('Reg method 1','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])

subplot(1,2,2), hold on
plot(filtertime,K2real,'b','LineWidth',2)
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])

snapnow;
delete(gcf);


%%
% Reg. method 2 leads to a minimum in the value of _S_ at _r_ = 1, which means that the mean of the eigenvalues of _C_ is chosen. This seems natural, and produces nice-looking filters with little high-frequency noise on them. The gain, though, is off, but can be corrected post-hoc by scaling the filter:
K2real_scaled = K2real*diagnostics2r.slope(diagnostics2r.bestfilter);
figure('outerposition',[0 0 350 350],'PaperUnits','points','PaperSize',[800 350]); hold on
plot(filtertime,K2real,'b','LineWidth',2), hold on
plot(filtertime,K2real_scaled,'k','LineWidth',2), hold on
legend Method2 Scaled
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])

snapnow;
delete(gcf);

%%
% This filter now has gain of exactly 1. Here is a figure of the prediction vs the actual data:
figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = convolve(time,PID,K2real,filtertime);

% censor initial prediction
fp(1:filter_length+1)=NaN;
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data ScaledPrediction

snapnow;
delete(gcf);

%%
% It seems to be following the trends of the data, but is off by a constant. That constant is the mean of the response, that we can add back to the prediction: 

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = fp + mean(f);
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend({'Data','Scaled Prediction + mean(response)'});

snapnow;
delete(gcf);


%%
% The r-square of this prediction is 
disp(rsquare(f,fp))

%%
% and the best prediction from all the regularisation factors has a r-square of:
disp(max(diagnostics2r.err))


%% Analysis of Linear Prediction - Does adding another linear filter improve prediction?
% An ideal linear filter should capture all the linear variation in the data. This means that if one were to construct a vector of residuals, a linear filter would be unable to predict the residuals from the data, and would lead to no improvement on the original filter. Is this true? 

%%
% First, we construct a vector of residuals by subtracting the linear prediction from the data.
f = f(:);
fp = fp(:);

res = f - fp;
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(2,1,1), hold on
plot(time,res,'k','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on')
ylabel('Residuals','FontSize',font_size)
subplot(2,1,2), hold on
plot(time,PID,'k','LineWidth',2)
ylabel('PID','FontSize',font_size)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on')
xlabel('Time (s)','FontSize',font_size)

snapnow;
delete(gcf);

%%
% we then construct a filter from the PID to these residuals to try to predict the residuals.
[Kres ,diagnostics_res] = FindBestFilter(PID(filter_length+2:end),res(filter_length+2:end));

figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[800 400]); hold on
plot(filtertime,Kres,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size,'XLim',[min(filtertime) max(filtertime)])
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: PID > residuals','FontSize',20)

snapnow;
delete(gcf);

%%
% and use this filter to predict the residuals from the stimulus. 
resp = convolve(time,PID,Kres,filtertime) + mean2(res);
figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(time,res,'k','LineWidth',2)
plot(time,resp,'b','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Residuals','FontSize',font_size)
legend Residuals 'Predicted Residuals'

snapnow;
delete(gcf);

%%
% Does adding this back to the prediction improve the prediction? The following plot shows the data (black), the linear prediction (red), and the linear prediction corrected by the prediciton of the residual (green).
fpr = fp + resp;
figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
plot(time,fpr,'g','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Residuals','FontSize',font_size)
legend Data 'Linear Prediction' 'Residual Corrected'

snapnow;
delete(gcf);

%%
% The r-square of the corrected prediction is
disp(rsquare(f(filter_length+2:end),fpr(filter_length+2:end)))

%% 
% and the r-square of the simple linear prediction is
disp(rsquare(f(filter_length+2:end),fp(filter_length+2:end)))

%%
% How does this make any sense? If this is true, then simply adding the filters together will yield a better filter, which we should have computed in the very beginning. The filter, the second filter computed from residuals, and their sum is shown in the panel on the left. On the right is a filter with lower regularisation. 

figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K2real,'k','LineWidth',2)
plot(filtertime,Kres,'r','LineWidth',2)
plot(filtertime,Kres+K2real,'g','LineWidth',2)
legend Filter ResidualFilter Sum
set(gca,'FontSize',font_size,'LineWidth',2,'box','on','XLim',[min(filtertime) max(filtertime)])
ylabel('Filter Amplitude (Hz)','FontSize',font_size)

K_lowreg = fitFilter2Data(PID,f,[],'reg',1e-1');
subplot(1,2,2), hold on
plot(0:3e-3:1,K_lowreg,'k','LineWidth',2)
title('Filter with low regularisation','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',2,'box','on','XLim',[min(filtertime) max(filtertime)])
xlabel('Time (s)','FontSize',font_size)

snapnow;
delete(gcf);

%%
% So what we are doing is simply undoing the work we did in regularising the filter. 


%       ########    ###     ######  ########     #######  ########   #######  ########  
%       ##         ## ##   ##    ##    ##       ##     ## ##     ## ##     ## ##     ## 
%       ##        ##   ##  ##          ##       ##     ## ##     ## ##     ## ##     ## 
%       ######   ##     ##  ######     ##       ##     ## ##     ## ##     ## ########  
%       ##       #########       ##    ##       ##     ## ##     ## ##     ## ##   ##   
%       ##       ##     ## ##    ##    ##       ##     ## ##     ## ##     ## ##    ##  
%       ##       ##     ##  ######     ##        #######  ########   #######  ##     ## 

%% Real Data 2: Flickering stimulus for a fast odor 
% We now attempt to back out a filter from a different data set, where the statistics of the stimulus are more tightly correlated, and where the neuron's response goes to zero frequently. This is what the data looks like: 


time = data(5).time;
PID = data(5).PID;
f = data(5).ORN;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(time,PID,'k','LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22])
ylabel('PID (a.u.)')

subplot(2,1,2), hold on
plot(time,f,'k','LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22])
ylabel('Firing rate (Hz)')

prettyFig;

snapnow;
delete(gcf);

%% 
% We now back out the filter using the methods described above and get: 

[K2real,  diagnostics2r, filtertime] = FindBestFilter(PID,f,[],'regtype=2;');

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
filtertime = filtertime*3e-3;
plot(filtertime,K2real,'k')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')
set(gca,'XLim',[min(filtertime) max(filtertime)])

prettyFig;
snapnow;
delete(gcf);

%%
% The following figure shows how the choice of regularisation parameter affects the quality of prediction:

PlotFilterDiagnostics2(diagnostics2r,marker_size,marker_size2,font_size,'Method 2')

snapnow;
delete(gcf);

%%
% And now we can compare the prediction to the actual response: 
figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = convolve(time,PID,K2real,filtertime);

% censor initial prediction
fp(1:filter_length+1)=NaN;
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data ScaledPrediction

prettyFig;

snapnow;
delete(gcf);

%%
% Once again, the linear prediction cannot account for a response with a non-zero mean. Adding the mean of the response back to the data, we get: 

fp = fp + mean(f);

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data ScaledPrediction

prettyFig;

snapnow;
delete(gcf);

%%
% There are some times when the prediction of firing rates goes below 0, which has no physical meaning. 


% ##    ##    ###    ######## ##     ## ########     ###    ##        ######  ######## #### ##     ## 
% ###   ##   ## ##      ##    ##     ## ##     ##   ## ##   ##       ##    ##    ##     ##  ###   ### 
% ####  ##  ##   ##     ##    ##     ## ##     ##  ##   ##  ##       ##          ##     ##  #### #### 
% ## ## ## ##     ##    ##    ##     ## ########  ##     ## ##        ######     ##     ##  ## ### ## 
% ##  #### #########    ##    ##     ## ##   ##   ######### ##             ##    ##     ##  ##     ## 
% ##   ### ##     ##    ##    ##     ## ##    ##  ##     ## ##       ##    ##    ##     ##  ##     ## 
% ##    ## ##     ##    ##     #######  ##     ## ##     ## ########  ######     ##    #### ##     ## 

%% Real Data 3: "Natual Stimuli" Responses 
% Now we try our filter estimation algorithms on a different dataset, where the neuron is driven by a very sparse, broadly fluctuating stimulus, as shown below: 

load('/local-data/DA-paper/mahmut_data.mat')

time = data(5).time;
PID = mean2(data(5).PID);
f = mean2(data(5).ORN);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(time,PID,'k','LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[10 22])
ylabel('PID (a.u.)')

subplot(2,1,2), hold on
plot(time,f,'k','LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[10 22])
ylabel('Firing rate (Hz)')
xlabel('Time (s)')

prettyFig;

snapnow;
delete(gcf);

%%
% We now attempt to fit a filter to this data. The following figure shows the best reconstructed filter for this dataset:

[K2real,  diagnostics2r, filtertime] = FindBestFilter(PID,f,[],'regtype=2;');

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
filtertime = filtertime*3e-3;
plot(filtertime,K2real,'k')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')
set(gca,'XLim',[min(filtertime) max(filtertime)])

prettyFig;

snapnow;
delete(gcf);

%%
% The following figure shows how the choice of regularisation parameter affects the quality of prediction:

PlotFilterDiagnostics2(diagnostics2r,marker_size,marker_size2,font_size,'Method 2')

snapnow;
delete(gcf);

%%
% And now we can compare the prediction to the actual response: 
figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = convolve(time,PID,K2real,filtertime);

% censor initial prediction
fp(1:filter_length+1)=NaN;
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-10 mean(time)+10],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data ScaledPrediction

prettyFig;

snapnow;
delete(gcf);

%%
% Once again, the linear prediction cannot account for a response with a non-zero mean. Adding the mean of the response back to the data, we get: 

fp = fp + mean(f);

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-10 mean(time)+10],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data ScaledPrediction

prettyFig;

snapnow;
delete(gcf);

%%
% There are some times when the prediction of firing rates goes below 0, which has no physical meaning. 



%% Version Info


pFooter;
