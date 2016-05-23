%% makeMSGplots
% helper scripts that plots I/O curves, gain, and kinetics for MSG data
% this is meant to be called by webers_generally_observed.m


function makeMSGplots(cdata)

% unpack data
v2struct(cdata)

% show some coarse statistics of the data
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,2,1), hold on
mean_stim = nanmean(PID(a:z,:));
std_stim = nanstd(PID(a:z,:));
plot(mean_stim,std_stim,'k+')
M = max([mean_stim(:); std_stim(:)]);
plot([0 M],[0 M],'k--')
xlabel('\mu_{Stimulus} (V)')
ylabel('\sigma_{Stimulus} (V)')
ff = fit(mean_stim(:),std_stim(:),'poly1');
l = plot(sort(mean_stim),ff(sort(mean_stim)),'r');
legend(l,['m=' oval(ff.p1)])

subplot(2,2,2), hold on
plot(mean_stim,nanmean(fA(a:z,:)),'r+')
plot(mean_stim,nanmean(fB(a:z,:)),'bo')
legend({'A','B'},'Location','southeast')
xlabel('\mu_{Stimulus} (V)')
ylabel('\mu_{Firing rate (Hz)}')

subplot(2,2,3), hold on
plot(mean_stim,nanstd(fA(a:z,:)),'r+')
plot(mean_stim,nanstd(fB(a:z,:)),'bo')
legend({'A','B'},'Location','northeast')
xlabel('\mu_{Stimulus} (V)')
ylabel('\sigma_{Firing rate (Hz)}')

subplot(2,2,4), hold on
plot(nanmean(fA(a:z,:)),nanstd(fA(a:z,:)),'r+')
plot(nanmean(fB(a:z,:)),nanstd(fB(a:z,:)),'bo')
legend({'A','B'},'Location','northeast')
xlabel('\mu_{Firing rate (Hz)}')
ylabel('\sigma_{Firing rate (Hz)}')

prettyFig


% show gain
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on

% make a nice colour scheme
c = parula(length(paradigm));

% also estimate gain using variances of stimulus and response

frac_var = NaN(width(PID),1);
frac_var_LFP = NaN(width(PID),1);
for i = 1:width(PID)
	try
		frac_var(i) = std(fA(a:z,i))/std(PID(a:z,i));
	catch
	end
	try
		frac_var_LFP(i) = std(LFP(a:z,i))/std(PID(a:z,i));
	catch
	end
end
frac_var_LFP = frac_var_LFP(:);
frac_var = frac_var(:);

% plot frac var of LFP
subplot(2,3,2), hold on


[~,idx]=sort(mean_stim,'ascend');


for i = 1:width(PID)
	ii = idx(i);
	plot(mean_stim(ii),frac_var_LFP(ii),'+','Color',c(i,:))
end

options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(nonnans(mean_stim),nonnans(frac_var_LFP),'power1',options);
plot(sort(mean_stim),cf(sort(mean_stim)),'r');
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('\sigma_{LFP}/\sigma_{Stimulus} (Hz/V)')

if exist('fA','var')
	% plot frac. var of firing rate
	subplot(2,3,5), hold on
	for i = 1:width(PID)
		ii = idx(i);
		plot(mean_stim(ii),frac_var(ii),'+','Color',c(i,:))
	end

	options = fitoptions(fittype('power1'));
	options.Lower = [-Inf -1];
	options.Upper = [Inf -1];
	cf = fit(nonnans(mean_stim),nonnans(frac_var),'power1',options);
	plot(sort(mean_stim),cf(sort(mean_stim)),'r');
	set(gca,'XScale','log','YScale','log')
	xlabel('Mean Stimulus (V)')
	ylabel('\sigma_{Firing Rate}/\sigma_{Stimulus} (Hz/V)')

	% show gain changes -- firing gain vs. mean stimulus
	subplot(2,3,4), hold on
	for i = 1:width(PID)
		ii = idx(i);
		plot(mean_stim(ii),fA_gain(ii),'+','Color',c(i,:))
	end

	% fit a power law with exponent -1
	mean_stim = mean_stim(:);
	fA_gain = fA_gain(:);
	options = fitoptions(fittype('power1'));
	options.Lower = [-Inf -1];
	options.Upper = [Inf -1];
	cf = fit(mean_stim(~isnan(fA_gain)),fA_gain(~isnan(fA_gain)),'power1',options);
	l = plot(sort(mean_stim),cf(sort(mean_stim)),'r');
	r2 = rsquare(nonnans(fA_gain),cf(nonnans(mean_stim)));
	set(gca,'XScale','log','YScale','log')
	xlabel('Mean Stimulus (V)')
	ylabel('ORN Gain (Hz/V)')
end

% show gain changes -- LFP gain vs. mean stimulus
if exist('fA','var')
	subplot(2,3,1), hold on
else
	subplot(1,3,1), hold on
end
for i = 1:width(PID)
	ii = idx(i);
	plot(mean_stim(ii),LFP_gain(ii),'+','Color',c(i,:))
end

% fit a power law with exponent -1
mean_stim = mean_stim(:);
LFP_gain = LFP_gain(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(mean_stim(~isnan(LFP_gain)),LFP_gain(~isnan(LFP_gain)),'power1',options);
l = plot(sort(mean_stim),cf(sort(mean_stim)),'r');
r2 = rsquare(nonnans(LFP_gain),cf(nonnans(mean_stim)));
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')

% do kinetics 

dt = 1e-3;
lag_LFP = NaN(1,width(PID));
lag_fA = NaN(1,width(PID));

for i = 1:width(PID)
	s = PID(a:z,i)-mean(PID(a:z,i)); s = s/std(s);
	try
		r = fA(a:z,i)-mean(fA(a:z,i)); r = r/std(r);
	catch
	end
	x = LFP(a:z,i)-mean(LFP(a:z,i)); x = x/std(x); x = -x;

	try
		temp = xcorr(r,s); temp = temp/(z-a);
		[~,lag_fA(i)] = max(temp);
	catch
	end

	temp = xcorr(x,s);  temp = temp/(z-a);
	[~,lag_LFP(i)] = max(temp);

end

try
	lag_fA = lag_fA - (z-a);
	lag_fA(lag_fA<0) = NaN; lag_fA(lag_fA>1e3) = NaN;
catch
end
lag_LFP = lag_LFP - (z-a);
lag_LFP(lag_LFP<0) = NaN; lag_LFP(lag_LFP>1e3) = NaN;



clear ax
if exist('fA','var')
	ax(1) = subplot(2,3,3); hold on
else
	ax(1) = subplot(1,3,3); hold on
end
plot(mean_stim,lag_LFP,'k+')
xlabel('Mean Stimulus (V)')
ylabel('LFP lag (ms)')

if exist('fA','var')
	ax(2) = subplot(2,3,6); hold on
	plot(mean_stim,lag_fA,'k+')
	xlabel('Mean Stimulus (V)')
	ylabel('Firing lag (ms)')
	set(ax(1),'YLim',[0 max([lag_LFP lag_fA])])
	set(ax(2),'YLim',[0 max([lag_LFP lag_fA])])
end

prettyFig
