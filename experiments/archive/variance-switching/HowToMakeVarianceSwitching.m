% How to make variance switching stimuli
% 
% created by Srinivas Gorur-Shandilya at 2:56 , 15 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% How to Make Variance Switching Stimuli
% This document shows you how to make Variance Switching Stimuli, where blocks of high-variance Gaussian flickering stimuli are interspersed with blocks of low-variance Gaussian flickering stimulus. We want to do this while keeping the mean constant throughout, and we also don't want any linear trends in the individual blocks. 


calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

load alldata
epoch_length=5;


%% Ansatz Solution
% A good guess on how to make this is to use two MFCs, which vary with different variances, and to use valves to switch between odourised airstreams from them. In the following figure, we show what this ansatz looks like: 

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,7,1:5), hold on
PID = alldata(6).data(2).PID(2,:);
PID = PID(1:10:end);
nepochs = ((length(PID)/epoch_length)*1e-3);
time = 1e-3*(1:length(PID));
plot(time,PID,'k');
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(1,7,6:7), hold on
PID = alldata(6).data(2).PID(2,1:10:end); PID = PID';
PID(1:epoch_length*1e3,:) = [];
PID(end-epoch_length*1e3+1:end,:) = [];
PID = PID(:);
s = repmat([ones(epoch_length*1e3,1); zeros(epoch_length*1e3,1)],nepochs/4,1);
s = ~logical(s);

[x,y] = hist(PID(s),50);
plot(x,y,'r')
[x,y] = hist(PID(~s),50);
plot(x,y,'b')
xlabel('Count')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% As we can see from this figure, we can quickly change between high and low variances, and there appears to be no trend in the data, but the mean is off. For some (unknown) reason, the mean when the variance is high is high, even though we constructed these stimuli so that the means would be the same. 

%%
% Instead of trying to solve this problem, we shall be lazy and work around by adding a free parameter that scales the mean setpoint when the variance is high. Using the powerful automation features of kontroller, we can go through a bunch of these free parameters (called fudge factors). This is in findFudgeFactor.m. 

%%
% The following figure shows the results of this parameter scan. On the left, we compare how the mean changes, and on the right, the variance. On the Y axis, we plot the ratio between the mean (or variance) during the high variance epoch and the low variance epoch. For our target stimulus, we want the value on the left plot to be close to 1. 


stimulus_mean = NaN*fudge_factor_range;
stimulus_mean_err = NaN*fudge_factor_range;
stimulus_std = NaN*fudge_factor_range;
stimulus_std_err = NaN*fudge_factor_range;

nepochs = ((length(PID)/epoch_length)*1e-3);

for i = 1:length(fudge_factor_range)
	PID = alldata(i).data(2).PID(:,1:10:end); PID = PID';
	PID(1:epoch_length*1e3,:) = [];
	PID(end-epoch_length*1e3+1:end,:) = [];
	PID = PID(:);
	s = repmat([ones(epoch_length*1e3,1); zeros(epoch_length*1e3,1)],nepochs/2,1);
	s = ~logical(s);

	hi_means = mean(reshape(PID(s),epoch_length*1e3,nepochs/2));
	lo_means = mean(reshape(PID(~s),epoch_length*1e3,nepochs/2));
	stimulus_mean(i) = mean(hi_means./lo_means);
	stimulus_mean_err(i) = std(hi_means./lo_means);

	hi_var = std(reshape(PID(s),epoch_length*1e3,nepochs/2));
	lo_var = std(reshape(PID(~s),epoch_length*1e3,nepochs/2));
	stimulus_std(i) = mean(hi_var./lo_var);
	stimulus_std_err(i) = std(hi_var./lo_var);

end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
errorbar(fudge_factor_range,stimulus_mean,stimulus_mean_err,'k')
ylabel('Ratio of means')
xlabel('Fudge Factor')

subplot(1,2,2), hold on
errorbar(fudge_factor_range,stimulus_std,stimulus_std_err,'k')
ylabel('Ratio of variances')
xlabel('Fudge Factor')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% Well that looks pretty good. What does the raw trace look like with the best fudge factor?

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,7,1:5), hold on
PID = alldata(1).data(2).PID(2,:);
PID = PID(1:10:end);
time = 1e-3*(1:length(PID));
plot(time,PID,'k');
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(1,7,6:7), hold on
PID = alldata(1).data(2).PID(2,1:10:end); PID = PID';
PID(1:epoch_length*1e3,:) = [];
PID(end-epoch_length*1e3+1:end,:) = [];
PID = PID(:);
s = repmat([ones(epoch_length*1e3,1); zeros(epoch_length*1e3,1)],nepochs/2,1);
s = ~logical(s);

[x,y] = hist(PID(s),50);
plot(x,y,'r')
[x,y] = hist(PID(~s),50);
plot(x,y,'b')
xlabel('Count')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end



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

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
