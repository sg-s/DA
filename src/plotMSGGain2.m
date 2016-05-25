%% makeMSGGain2
% helper scripts that plots I/O curves, gain, and kinetics for MSG data
% this is meant to be called by webers_generally_observed.m


function makeMSGGain2(cdata,ax)

% unpack data
v2struct(cdata)

% make a nice colour scheme
c = parula(length(paradigm));
mean_stim = nanmean(PID(a:z,:));
[~,idx]=sort(mean_stim,'ascend');


% also estimate gain using variances of stimulus and response

frac_var = NaN(width(PID),1);
frac_var_LFP = NaN(width(PID),1);
for i = 1:width(PID)
	try
		frac_var(i) = nanstd(fA(a:z,i))/nanstd(PID(a:z,i));
	catch
	end
	try
		frac_var_LFP(i) = nanstd(LFP(a:z,i))/nanstd(PID(a:z,i));
	catch
	end
end
frac_var_LFP = frac_var_LFP(:);
frac_var = frac_var(:);

% plot frac var of LFP

[~,idx]=sort(mean_stim,'ascend');
for i = 1:width(PID)
	ii = idx(i);
	plot(ax(1),mean_stim(ii),frac_var_LFP(ii),'+','Color',c(i,:))
end

x = mean_stim(~isnan(frac_var_LFP));
y = frac_var_LFP(~isnan(frac_var_LFP));
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on

plot(ax(1),sort(mean_stim),cf(sort(mean_stim)),'r');
set(ax(1),'XScale','log','YScale','log')
xlabel(ax(1),'Mean Stimulus (V)')
ylabel(ax(1),'\sigma_{LFP}/\sigma_{Stimulus} (Hz/V)')


% plot frac. var of firing rate
for i = 1:width(PID)
	ii = idx(i);
	plot(ax(2),mean_stim(ii),frac_var(ii),'+','Color',c(i,:))
end

x = mean_stim(~isnan(frac_var));
y = frac_var(~isnan(frac_var));
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on

plot(ax(2),sort(mean_stim),cf(sort(mean_stim)),'r');
set(ax(2),'XScale','log','YScale','log')
xlabel(ax(2),'Mean Stimulus (V)')
ylabel(ax(2),'\sigma_{Firing Rate}/\sigma_{Stimulus} (Hz/V)')
