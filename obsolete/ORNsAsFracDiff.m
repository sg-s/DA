% ORNsAsFracDiff.m
% ORNs as fractional differentiators
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
% check if there is a data cache and load it
if exist(strcat(mfilename,'.mat'))
	load(strcat(mfilename,'.mat'))

else
	load('/local-data/DA-paper/flickering-stim/data.mat')
end




%% ORNs as fractional differentiators
% Are ORNs fractional differentiators? 

%% Responses to binary flickering stimuli
% The following figure shows a typical response of a neuron (ab3A) to flickering binary odor stimulus (top panel). Note the transient artefact at the end of every pulse that may be the effect of the electrical noise generated by valve switching. 

td = 3;
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID)
ylabel('Stimulus (V)')
set(gca,'XLim',[15 20])
subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
set(gca,'XLim',[15 20])

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% We can fit a fractional differentiator to this data set, and compare it to a LN model: 

if ~exist('fd_pred')
	[alpha,d,fd_pred] = FitFractionalDModel(data(td).PID,data(td).ORN);
end


% build a simple linear model
[data(td).K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
data(td).LinearFit = data(td).LinearFit + mean(data(td).ORN);

xdata = data(td).LinearFit;
ydata = data(td).ORN;

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
% save this for later
LN_pred = hill(x,data(td).LinearFit);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,fd_pred)
plot(data(td).time,LN_pred)
set(gca,'XLim',[15 20],'YLim',[-5 1.1*max(data(td).ORN)])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','Frac. Diff.','LN Model'})
PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% We see that the fractional differentiator is responding strongly to the transient artefact of the valve switching, which may affect the choice of order of fractional differentiation. To work around this, we sand the artefacts on the PID down:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).PID,'k')
data = SandPIDArtifacts(data,td);
plot(data(td).time,data(td).PID2,'r')
legend({'Original Data','Artifacts removed'})
set(gca,'XLim',[15 20])
xlabel('Time (s)')
ylabel('Stimulus (V)')
PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% We now recompute the best order of a fractional differentiator and fit the data. To keep the comparison fair, we also recompute the LN model using the artifact-free PID trace. 

% build a simple linear model
[data(td).K,~,filtertime] = FindBestFilter(data(td).PID2(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = convolve(data(td).time,data(td).PID2,data(td).K,data(td).filtertime);
data(td).LinearFit = data(td).LinearFit + mean(data(td).ORN);

xdata = data(td).LinearFit;
ydata = data(td).ORN;

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
% save this for later
LN_pred = hill(x,data(td).LinearFit);

if ~exist('order_frac_diff')
	[order_frac_diff,~,fd_pred2] = FitFractionalDModel(data(td).PID2,data(td).ORN);
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,fd_pred2)
plot(data(td).time,LN_pred)
legend({'Data','Frac. Diff.','LN Model'})
set(gca,'XLim',[15 20],'YLim',[-5 1.1*max(data(td).ORN)])
xlabel('Time (s)')
ylabel('Firing Rate (Hz) (V)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% 
% Note that the output of the fractional differentiation operation is only rescaled by the mean and thresholded to throw out negative values. There is no additional output non-linearity as in the LN model. 

%%
% The r-square of the LN model is:

disp(rsquare(LN_pred,data(td).ORN))

%%
% cf. the r-square of the fractional diffentiator:

disp(rsquare(fd_pred2,data(td).ORN))

%%
% And the order of the fractional differentiation used here is: 

disp(order_frac_diff)



%% Responses to pulses of different heights
% The responses of ORNs to odor pulses of different heights has been well characterised, without (left) and on top of a background (right). In the following figure, the panels on the top show the stimulus presented to the neurons, and the bottom panels show the the response of the neurons. The odor presented is ethyl acetate, and the neuron is ab3A. 

if ~exist('dr_data')
	dr_data=load('/local-data/DA-paper/dose-response/carlotta/dr_data.mat');
	dr_data = dr_data.data;


	% throw out weird looking data
	dr_data(5).PID(:,8:10) = [];
	dr_data(5).ORN(:,8:10) = [];
	dr_data(6).PID(:,9:10) = [];
	dr_data(6).ORN(:,9:10) = [];
end




figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,2,1), hold on
title('Pulses w/o background')
plot(dr_data(7).time,dr_data(7).PID)
ylabel('Stimulus (V)')
set(gca,'XLim',[0 3])
subplot(2,2,2), hold on
title('Pulses with background')
plot(dr_data(5).time,dr_data(5).PID)
set(gca,'XLim',[0 3])
subplot(2,2,3), hold on
plot(dr_data(7).time,dr_data(7).ORN)
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[0 3])
subplot(2,2,4), hold on
plot(dr_data(5).time,dr_data(5).ORN)
set(gca,'XLim',[0 3])
xlabel('Time (s)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% In the following section, we will attempt to fit a fractional differentiator to each pulse, and extract a order for each pulse. We start by fitting dr_data to the pulses in the absence of a background. 



a = 350;
z = 600;
if ~isfield(dr_data,'fd_pred')
	for td = 5:7
		dr_data(td).PID(1,:) = 0;
		dr_data(td).fd_pred = dr_data(td).ORN;
		dr_data(td).alpha = NaN(width(dr_data(td).ORN),1);
		for i = 1:width(dr_data(td).ORN)
			textbar(i,width(dr_data(td).ORN))
			[dr_data(td).alpha(i),d,dr_data(td).fd_pred(:,i)] =  FitFD2DoseResponse(dr_data(td).PID(:,i), dr_data(td).ORN(:,i),a,z);
			
		end
		clear i
		for i = 1:width(dr_data(td).ORN)
			dr_data(td).fd_pred(isnan(dr_data(td).fd_pred(:,i)),i) = 0;
			dr_data(td).fd_pred(:,i) = filtfilt(ones(1,10)/10,1,dr_data(td).fd_pred(:,i));
		end
	end
end

td=7;
dr_data(td).alpha(1) = 0;
figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
subplot(1,3,1), hold on
title('Data')
plot(dr_data(td).time,dr_data(td).ORN)
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
set(gca,'XLim',[0.5 3])
subplot(1,3,2), hold on
title('Frac. Diff. Prediction')
plot(dr_data(td).time,dr_data(td).fd_pred)
set(gca,'XLim',[0.5 3])
subplot(1,3,3), hold on
plot(max(dr_data(td).PID),dr_data(td).alpha)
set(gca,'XScale','log')
ylabel('Order of frac. diff.')
xlabel('Stimulus Peak (V)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% We now repeat the fit for responses on top of a background, and then compare how the order of fractional differentiation varies with background concentration (right panel)

warning off
td=5;
dr_data(td).alpha(1) = 0;
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,3,1), hold on
title('Pulses on background 2')
plot(dr_data(td).time,dr_data(td).ORN)
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
set(gca,'XLim',[0.5 3])
subplot(1,3,2), hold on
title('Frac. Diff. Prediction')
plot(dr_data(td).time,dr_data(td).fd_pred)
set(gca,'XLim',[0.5 3])
subplot(1,3,3), hold on

for td = 7:-1:5
	dr_data(td).alpha(dr_data(td).alpha >.9) = 0;
	plot(max(dr_data(td).PID),dr_data(td).alpha)
end
set(gca,'XScale','log')
ylabel('Order of frac. diff.')
xlabel('Stimulus Peak (V)')
legend({'No Back','back 1','back 2'},'Location','southeast')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end
warning on



%% Responses to long steps of odor 
% How well does fractional differentiation explain the responses to long steps of odor? In the following section we fit a fractional differentiation operator to long steps of odor. 

% load the data
load('/local-data/DA-paper/long-steps/mahmut/mahmut_ls.mat')
figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
z = 8000;
subplot(2,1,1), hold on
for i = 1:length(ls_data)
	plot(ls_data(i).time(1:z),mean(ls_data(i).PID(1:z,:)'))
end
subplot(2,1,2), hold on
for i = 1:length(ls_data)
	plot(ls_data(i).time(1:z),mean(ls_data(i).ORN(1:z,:)'))
end

% cache data
being_published = 0;
save(strcat(mfilename,'.mat'))

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





