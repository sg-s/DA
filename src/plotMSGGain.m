%% plotMSGGain.m
% helper scripts that plots I/O curves, gain, and kinetics for MSG data
% this is meant to be called by webers_generally_observed.m


function plotMSGGain(cdata,ax)

% unpack data
v2struct(cdata)


% make a nice colour scheme
c = parula(length(paradigm));
mean_stim = nanmean(PID(a:z,:));
[~,idx]=sort(mean_stim,'ascend');

% show gain changes -- firing gain vs. mean stimulus
for i = 1:width(PID)
	ii = idx(i);
	plot(ax(2),mean_stim(ii),fA_gain(ii),'+','Color',c(i,:))
end

% fit a power law with exponent -1
mean_stim = mean_stim(:);
fA_gain = fA_gain(:);

x = mean_stim(~isnan(fA_gain));
y = fA_gain(~isnan(fA_gain));
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x,y,'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on

l = plot(ax(2),sort(mean_stim),cf(sort(mean_stim)),'r');
r2 = rsquare(nonnans(fA_gain),cf(nonnans(mean_stim)));
set(ax(2),'XScale','log','YScale','log')
xlabel(ax(2),'Mean Stimulus (V)')
ylabel(ax(2),'ORN Gain (Hz/V)')


% show gain changes -- LFP gain vs. mean stimulus
for i = 1:width(PID)
	ii = idx(i);
	plot(ax(1),mean_stim(ii),LFP_gain(ii),'+','Color',c(i,:))
end

% fit a power law with exponent -1
mean_stim = mean_stim(:);
LFP_gain = LFP_gain(:);

x = mean_stim(~isnan(LFP_gain));
y = LFP_gain(~isnan(LFP_gain));
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x,y,'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on

l = plot(ax(1),sort(mean_stim),cf(sort(mean_stim)),'r');
set(ax(1),'XScale','log','YScale','log')
xlabel(ax(1),'Mean Stimulus (V)')
ylabel(ax(1),'LFP Gain (mV/V)')

