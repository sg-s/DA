% Variance Flickering
% analysis of response to stimulus where variance fluctuates over time
% 
% created by Srinivas Gorur-Shandilya at 2:41 , 07 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Variance Flickering
% In this document, we analyse signals where the mean is constant, but the variance changes over time. In practice, what we do is oscillate the MFC very quickly (~20ms), and configure the parameters (P->1000, I->0, D->6000) for maximum speed, at the cost of reproducibility. We then choose the variance from one of two values, randomly, every 500ms. 

% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%%
% The following figure shows an example trace showing the flow rates (top) through the odour, and the corresponding stimulus measurements (bottom). Note the periods of alternating variance. 


load('/local-data/DA-paper/variance-flickering/2015_08_07_variance_flickering_test2_var_.25.mat')

% recover the variance classes
T = 60; % length of trial
dt = 1e-4; % sampling rate
tau2 = 5e-1; % 100ms, switching time of variance. 

msteps = floor(T/tau2);

% set rand stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 

variances = double(rand(msteps,1) > .5);
variances = repmat(variances,1,(T/dt)/msteps)';
variances = variances(:);

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,1,1), hold on
time = 1e-4*(1:length(data.PID));
plot(time(1:100:end),100*data.MFC500(2,1:100:end),'k')
ylabel('Flow Rate (mL/min)')
set(gca,'XLim',[30 40])
subplot(2,1,2), hold on
plot(time(1:100:end),data.PID(2,1:100:end),'k')
ylabel('Stimulus (V)')
xlabel('Time (s)')
set(gca,'XLim',[30 40],'YLim',[0 1])

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% In the following figure, we plot the distributions of the control signals (left), the flow rates (middle) and the PID (right) at times when the variance is low and when the variance is high. Error bars indicate standard error of the mean.  


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,3,1), hold on
control_signals = ControlParadigm.Outputs(1,variances<0.5);
dil = (control_signals*100)./((control_signals*100)+2000);
dil = dil*100;
[y,x] = hist(dil,1:0.1:8);
y = y/sum(y);
plot(x,y,'b');
control_signals = ControlParadigm.Outputs(1,variances>0.5);
dil = (control_signals*100)./((control_signals*100)+2000);
dil = dil*100;
[y,x] = hist(dil,1:0.1:8);
y = y/sum(y);
plot(x,y,'r')
legend({'Low Variance','High Variance'})
xlabel('Setpoint Dilution (%)')
ylabel('Probability')


subplot(1,3,2), hold on
x = 0:0.05:2;
y = zeros(length(x),width(data.MFC500)-1);
for i = 2:width(data.MFC500)
	data.MFC500(i,1:10e4) = NaN;
	stim = data.MFC500(i,variances<.5);
	stim(isnan(stim)) = [];
	y(:,i-1) = hist(stim,x);
	y(:,i-1) = y(:,i-1)/sum(y(:,i-1));
end
shadedErrorBar(100*x,mean2(y),sem(y),{'Color','b'})
for i = 2:width(data.MFC500)	
	stim = data.MFC500(i,variances>.5);
	stim(isnan(stim)) = [];
	y(:,i-1) = hist(stim,x);
	y(:,i-1) = y(:,i-1)/sum(y(:,i-1));
end
shadedErrorBar(100*x,mean2(y),sem(y),{'Color','r'})

xlabel('Flow Rate (mL/min)')

subplot(1,3,3), hold on
x = 0:0.05:1.5;
y = zeros(length(x),width(data.PID)-1);
for i = 2:width(data.PID)
	data.PID(i,1:10e4) = NaN;
	stim = data.PID(i,variances<.5);
	stim(isnan(stim)) = [];
	y(:,i-1) = hist(stim,x);
	y(:,i-1) = y(:,i-1)/sum(y(:,i-1));
end
shadedErrorBar(x,mean2(y),sem(y),{'Color','b'})
for i = 2:width(data.PID)	
	stim = data.PID(i,variances>.5);
	stim(isnan(stim)) = [];
	y(:,i-1) = hist(stim,x);
	y(:,i-1) = y(:,i-1)/sum(y(:,i-1));
end
shadedErrorBar(x,mean2(y),sem(y),{'Color','r'})

xlabel('Stimulus (V)')

PrettyFig()

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
disp(DataHash(strcat(mfilename,'.m'),Opt))

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

