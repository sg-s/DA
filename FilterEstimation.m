%% Filter Estimation
% How do we best calculate the filters? A summary of different ways of estimating the filters and results from each. 

%%

% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;

%% Overview of different regularisation methods. 
% The filter _K_ is given by 
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
% Synthetic data is prepared using Gaussian random inputs and an exponential filter, and the output is the convolution of the input with the exponential filter, with 10% noise. 
a = randn(1,10000);
filtertime = 0:3e-3:3e-3*333;
filter_length = 333;
Kexp = exp(-10*filtertime);
b = filter(Kexp,1,a) + 0.1*randn(1,10000);

%%
% The following figure shows filters estimated using the three methods shown above, with zero regularisation. The actual filter is in black, and each method is coloured differently. 

figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
K1 = FitFilter2Data(a,b,'reg=0;','regtype=1;');
plot(filtertime,Kexp,'k','LineWidth',2)
plot(filtertime,K1,'r','LineWidth',2)
title('Reg method 1 (r=0)','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

subplot(1,2,2), hold on
K2 = FitFilter2Data(a,b,'reg=0;','regtype=2;');
plot(filtertime,Kexp,'k','LineWidth',2)
plot(filtertime,K2,'b','LineWidth',2)
title('Reg method 2 (r=0)','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')


%% Synthetic Data: effect of regularisation 
% what is the effect of increasing _r_, the regularisation parameter on each of these methods? Here we systematically vary _r_, and pick the best filter that minimises 
% 
% $$ \arg\min_{r}\left\Vert \frac{e(r)-1}{1},\frac{S(r)-\min(S)}{\min(S)},\frac{\left|m(r)-1\right|}{1}\right\Vert _{\infty} $$
% 
% where _e_ is the coefficient of determination (r-square) of the prediction w.r.t to the actual data, _S_ is the sum of absolute values of the derivative of the filter _K_, and _m_ is the slope of the best fit of the prediction to the data. Minimising (_e_ - 1) minimises the error of the fit. Minimising _S_ minimises the high-frequency components in the filter, and prevents over-fitting. We also want the slope _m_ to be as close to 1 as possible. If _S_ has a minimum, that value of _r_ is automatically chosen. 

clear K1 K2 
regstr = 'regtype=1;';
if ~(exist('K1best') == 1)
	[K1best, diagnostics1] = FindBestFilter(a,b,filter_length,[1e-5 10],regstr);
end

regstr = 'regtype=2;';
if ~(exist('K2best') == 1)
	[K2best, diagnostics2] = FindBestFilter(a,b,filter_length,[1e-5 10],regstr);
end



%%
% In each of the plots below, the variation of error, filter sum, filter height, gain, and condition number of _C_ is show with _r_, the regularisation factor, in units of $\mu$. The error is calculated from the prediction to the actual output. The filter sum is the sum of the absolute values of the filter. Very high values here indicate the filter dominated by high-frequency components. The filter height is the peak of the filter. The gain is the slope of the best fit line of the data to the prediction. The condition number is the ratio of the largest to the smallest eigenvalue of _C_. Values of the condition number close to 1 mean that the matrix is easier to invert.  The figure below shows how varying _r_ affects reg method 1:
PlotFilterDiagnostics2(diagnostics1,marker_size,marker_size2,font_size,'Method 1')

%%
% the figure below shows how varying _r_ affects reg method 2:
PlotFilterDiagnostics2(diagnostics2,marker_size,marker_size2,font_size,'Method 2')

%%
% The best filters are shown below for the two reg. methods:
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K1best,'r','LineWidth',2)
title('Reg method 1','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

subplot(1,2,2), hold on
plot(filtertime,K2best,'b','LineWidth',2)
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')


%% Synthetic Data 2: effect of regularisation 
% We now create a new synthetic dataset, identical to the old one, except we filter the white noise input with some boxcar filter to remove high frequency components from the input. We find the best filters as before, and look at how regularisation parameter affects the choice of filter. 
a2 = filter(ones(1,10),1,a)/sqrt(10);
b2 = filter(Kexp,1,a2);

regstr = 'regtype=1;';
if ~(exist('K1best_f') == 1)
	[K1best_f, diagnostics1f] = FindBestFilter(a2,b2,filter_length,[1e-5 10],regstr);
end

regstr = 'regtype=2;';
if ~(exist('K2best_f') == 1)
	[K2best_f, diagnostics2f] = FindBestFilter(a2,b2,filter_length,[1e-5 10],regstr);
end

%%
% The best filters are shown below for the two reg. methods:
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K1best_f,'r','LineWidth',2)
title('Reg method 1','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

subplot(1,2,2), hold on
plot(filtertime,K2best_f,'b','LineWidth',2)
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

%%
% the figure below shows how varying _r_ affects reg method 1:
PlotFilterDiagnostics2(diagnostics1f,marker_size,marker_size2,font_size,'Method 1')

%%
% the figure below shows how varying _r_ affects reg method 2:
PlotFilterDiagnostics2(diagnostics2f,marker_size,marker_size2,font_size,'method 2')



%% Real Data: effect of regularisation 
% We use PID traces as the input, and the firing rate of a ORN as the output. The data looks like this:

% plot the data
load sample_data.mat
time = 3e-3:3e-3:3e-3*length(PID);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(time,PID,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[8 12])
ylabel('PID (a.u.)')

subplot(2,1,2), hold on
plot(time,f,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[8 12])
ylabel('Firing rate (Hz)')

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

%%
% Regularisation adds a large diagonal matrix to this, in effect making it look more like a diagonal matrix. Here, the mean of the eigenvalues of _C_ has been added to the diagonal. This is essentially what Method 1 does.  
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
r = trace(C)/length(C);
imagesc(C+r*eye(length(C)));
colorbar, axis tight
set(gca,'LineWidth',2,'FontSize',font_size)

%%
% The effect of regularisation using method 2 is shown on the correlation matrix below. The same _r_ is used, i.e., the mean of the eigenvalues of _C_. 
C2 = (C + r*eye(filter_length+1))*trace(C)/(trace(C) + r*filter_length);
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
imagesc(C2);
colorbar, axis tight
set(gca,'LineWidth',2,'FontSize',font_size)



%%
% Once again, we repeat the filter extraction. 
regstr = 'regtype=1;';
if ~(exist('K1real') == 1)
	[K1real,  diagnostics1r] = FindBestFilter(PID,f,filter_length,[1e-5 10],regstr,'Method 1');
end

regstr = 'regtype=2;';
if ~(exist('K2real') == 1)
	[K2real,  diagnostics2r] = FindBestFilter(PID,f,filter_length,[1e-5 10],regstr,'method 2');
end


%%
% the figure below shows how varying _r_ affects reg method 1:
PlotFilterDiagnostics2(diagnostics1r,marker_size,marker_size2,font_size,'Method 1')

%%
% the figure below shows how varying _r_ affects reg method 2:
PlotFilterDiagnostics2(diagnostics2r,marker_size,marker_size2,font_size,'Method 2')


%%
% and the best filters look like: 

figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K1real,'r','LineWidth',2)
title('Reg method 1','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

subplot(1,2,2), hold on
plot(filtertime,K2real,'b','LineWidth',2)
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

%%
% Reg. method 2 leads to a minimum in the value of _S_ at _r_ = 1, which means that the mean of the eigenvalues of _C_ is chosen. This seems natural, and produces nice-looking filters with little high-frequency noise on them. The gain, though, is off, but can be corrected post-hoc by scaling the filter:
K2real_scaled = K2real*diagnostics2r.slope(diagnostics2r.bestfilter);
figure('outerposition',[0 0 350 350],'PaperUnits','points','PaperSize',[800 350]); hold on
plot(filtertime,K2real,'b','LineWidth',2), hold on
plot(filtertime,K2real_scaled,'k','LineWidth',2), hold on
legend Method2 Scaled
title('Reg method 2','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on')

%%
% This filter now has gain of exactly 1. The prediction is pretty good: 
figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = filter(K2real_scaled,1,PID-mean(PID)) + mean(f);
% censor initial prediction
fp(1:filter_length+1)=NaN;
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data ScaledPrediction

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

%%
% we then construct a filter from the PID to these residuals to try to predict the residuals.
[Kres ,diagnostics_res] = FindBestFilter(PID(filter_length+2:end),res(filter_length+2:end));

figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[800 400]); hold on
plot(filtertime,Kres,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size,'XLim',[min(filtertime) max(filtertime)])
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: PID > residuals','FontSize',20)

%%
% and use this filter to predict the residuals from the stimulus. 
resp = filter(Kres,1,PID-mean(PID)) + mean2(res);
figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(time,res,'k','LineWidth',2)
plot(time,resp,'b','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Residuals','FontSize',font_size)
legend Residuals 'Predicted Residuals'

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
plot(filtertime,K,'k','LineWidth',2)
plot(filtertime,Kres,'r','LineWidth',2)
plot(filtertime,Kres+K,'g','LineWidth',2)
legend Filter ResidualFilter Sum
set(gca,'FontSize',font_size,'LineWidth',2,'box','on','XLim',[min(filtertime) max(filtertime)])
ylabel('Filter Amplitude (Hz)','FontSize',font_size)

K_lowreg = FitFilter2Data(PID,f,[],'reg=1e-1;');
subplot(1,2,2), hold on
plot(filtertime,K_lowreg,'k','LineWidth',2)
title('Filter with low regularisation','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',2,'box','on','XLim',[min(filtertime) max(filtertime)])
xlabel('Time (s)','FontSize',font_size)

%%
% So what we are doing is simply undoing the work we did in regularising the filter. 


