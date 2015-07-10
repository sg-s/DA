% Unused code from Analysis_January.m
% no longer relevant.
% 
% created by Srinivas Gorur-Shandilya at 2:47 , 10 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% 
% For a sanity check, let's construct an exponential filter like this:
figure('outerposition',[0 0 350 350],'PaperUnits','points','PaperSize',[800 350]); hold on
Kexp = exp(-10*filtertime);
plot(filtertime,Kexp,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)

%%
% And let's convolve it with the PID to generate a fake neuron response, and add 10% Gaussian noise to it. Now, we can extract the (fake) filter once more using these techniques from this fake dataset. 
ffake = filter(Kexp,1,PID) + 0.1*randn(1,length(time));
if ~(exist('Kfake') == 1)
	[Kfake Kdamonfake diagnostics_fake] = FindBestFilter(PID,ffake,filter_length);
end
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,Kexp,'LineWidth',2)
plot(filtertime,Kfake,'r','LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size,'YLim',[0 1])
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
title('Chichlinsky')

subplot(1,2,2), hold on
plot(filtertime,Kexp,'LineWidth',2)
plot(filtertime(2:end),Kdamonfake,'r','LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size,'YLim',[0 1])
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
title('Damon')

%%
% What are the slopes of best fit from these filters? For Chichilnisky's filter, it is: 
fpfake = filter(Kfake,1,PID);
[fitfake ~] = fit(fpfake(filter_length+2:end)',ffake(filter_length+2:end)','Poly1');
disp(fitfake.p1)

%%
% For Damon's filter, it is:
fpfake = filter(Kdamonfake,1,PID);
[fitfake ~] = fit(fpfake(filter_length+2:end)',ffake(filter_length+2:end)','Poly1');
disp(fitfake.p1)

%%
% Once again, we can see how filter properties vary with the free parameter _r_.

bf = diagnostics_fake.C.bestfilter;
figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(2,4,1), hold on
plot(diagnostics_fake.C.reg,diagnostics_fake.C.err,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.C.reg(bf),diagnostics_fake.C.err(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
ylabel('Chichilnisky','FontSize',font_size)
title('Error','FontSize',font_size)
axis square

subplot(2,4,2), hold on
plot(diagnostics_fake.C.reg,diagnostics_fake.C.filter_sum,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.C.reg(bf),diagnostics_fake.C.filter_sum(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('\Sigma|K|','FontSize',font_size)
axis square

subplot(2,4,3), hold on
plot(diagnostics_fake.C.reg,diagnostics_fake.C.filter_height,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.C.reg(bf),diagnostics_fake.C.filter_height(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Filter Height','FontSize',font_size)
axis square

subplot(2,4,4), hold on
plot(diagnostics_fake.C.reg,diagnostics_fake.C.slope,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.C.reg(bf),diagnostics_fake.C.slope(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Gain','FontSize',font_size)
axis square

bf = diagnostics_fake.D.bestfilter;
subplot(2,4,5), hold on
plot(diagnostics_fake.D.reg,diagnostics_fake.D.err,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.D.reg(bf),diagnostics_fake.D.err(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
ylabel('Damon','FontSize',font_size)
title('Error','FontSize',font_size)
axis square

subplot(2,4,6), hold on
plot(diagnostics_fake.D.reg,diagnostics_fake.D.filter_sum,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.D.reg(bf),diagnostics_fake.D.filter_sum(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('\Sigma|K|','FontSize',font_size)
axis square

subplot(2,4,7), hold on
plot(diagnostics_fake.D.reg,diagnostics_fake.D.filter_height,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.D.reg(bf),diagnostics_fake.D.filter_height(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Filter Height','FontSize',font_size)
axis square

subplot(2,4,8), hold on
plot(diagnostics_fake.D.reg,diagnostics_fake.D.slope,'.','MarkerSize',marker_size2)
plot(diagnostics_fake.D.reg(bf),diagnostics_fake.D.slope(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Gain','FontSize',font_size)
axis square

%% 
% From the variation of _r_ we see that, for Chichlinsky's filter, it is possible to minimise the error, the absolute filter sum, get the filter height correctly, and get the slope = 1. However, this happens in an intermediate position of _r_, and there is no function of _r_ that is at an extremum here. This is why we used the soft optimisation described earlier, where all these were simultaneously minimised. 



%%
% What happens in the ideal case when there is no additive noise? The same procedure is repeated, and now the filter is now exactly reconstructed with no regularisation. 
ffake = filter(Kexp,1,PID);
Kfake2 = FitFilter2Data(PID,ffake,filter_length,0);
figure('outerposition',[0 0 350 350],'PaperUnits','points','PaperSize',[800 350]); hold on
plot(filtertime,Kexp,'LineWidth',2)
plot(filtertime,Kfake2,'r','LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)
axis square
legend ActualFilter Estimation

%%
% and the slope of fit of prediction to data is:
fpf = filter(Kfake2,1,PID);
fpf(1:filter_length) = NaN;
[fall gof2] = fit(fpf(filter_length+2:end)',ffake(filter_length+2:end)','Poly1');
disp(fall.p1)




%%
% The autocorrelation function of the gain shows that it is very tightly constrained, almost as much as the valve, and much less than the PID, which we are trying to use to predict it. 
cf = xcorr(f(filter_length+2:end)-mean(f),'unbiased'); cf=cf/max(cf((length(cf)+1)/2-100:(length(cf)+1)/2+100));
cp = xcorr(PID(filter_length+2:end)-mean(PID),'unbiased'); cp=cp/max(cp((length(cp)+1)/2-100:(length(cp)+1)/2+100));
cv = xcorr(Valve(filter_length+2:end)-mean(Valve),'unbiased'); cv=cv/max(cv((length(cv)+1)/2-100:(length(cv)+1)/2+100));
cg = xcorr(gain(filter_length+2:end)-mean(gain(filter_length+2:end)),'unbiased');  cg=cg/max(cg((length(cg)+1)/2-100:(length(cg)+1)/2+100));
cn2 = xcorr(randn(1,length(f(filter_length+2:end))),'unbiased'); cn2 = cn2/max(cn2);
ctime =  mean(diff(time)):mean(diff(time)):length(cf)*mean(diff(time)); 
ctime = ctime - mean(ctime);
figure('outerposition',[10 10 850 400],'PaperUnits','points','PaperSize',[850 400]); hold on
subplot(1,2,1), hold on
plot(ctime,cf,'b','LineWidth',2)
plot(ctime,cp,'r','LineWidth',2)
plot(ctime,cv,'g','LineWidth',2)
plot(ctime,cg,'k','LineWidth',2)
set(gca,'box','on','XLim',[-0.5 0.5],'FontSize',font_size,'LineWidth',2,'YLim',[0 1.5])
xlabel('Time (s)','FontSize',font_size)
ylabel('Normalised autocorrelation','FontSize',font_size)
legend ORN PID Valve Gain 
subplot(1,2,2), hold on
plot(ctime,cn2,'r','LineWidth',2)
set(gca,'box','on','XLim',[-0.5 0.5],'FontSize',font_size,'LineWidth',2,'YLim',[0 1.5])
xlabel('Time (s)','FontSize',font_size)
legend GaussianNoise
