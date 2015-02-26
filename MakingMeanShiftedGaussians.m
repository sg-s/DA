% MakingMeanShiftedGaussians.m
% in this document, we investigate methods to make mean shifted gaussians. 
% 
% created by Srinivas Gorur-Shandilya at 9:39 , 26 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end


%% Section 1: Presenting uniformly distributed noise 
% In this section, we drive the MFC with some uniformly distributed noise with a 100ms correlation time, being careful to choose the noise in the dilution space. 

load('/local-data/DA-paper/mean-shifted-gaussians/2015_02_26_odour_test_1.mat')
% subsample data
time = 1e-4*(1:length(data(2).PID));
t = time(1:10:end);
PID = zeros(width(data(2).PID),length(data(2).PID)/10);
MFC = zeros(width(data(2).PID),length(data(2).PID)/10);
MFC_Control = (ControlParadigm(2).Outputs(2,1:10:end));
for i = 1:width(data(2).PID)
	PID(i,:) = interp1(time,data(2).PID(i,:),t);
	MFC(i,:) = interp1(time,data(2).MFC500(i,:),t);
end
clear time
time = t;

figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on


subplot(2,9,1:5), hold on
plot(time,MFC')
set(gca,'XLim',[0 60])
xlabel('Time (s)')
ylabel('MFC Flow V)')

subplot(2,9,6:7), hold on
[r2,s] = rsquare(MFC);
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,8:9), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


subplot(2,9,10:14), hold on
plot(time,PID')
set(gca,'XLim',[0 60])
xlabel('Time (s)')
ylabel('PID (V)')

subplot(2,9,15:16), hold on
[r2,s] = rsquare(PID);
imagescnan(r2)
caxis([0 1])
axis image
colorbar
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,17:18), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% First off, there is a long initial transient in the PID for the first trial, but subsequent trials are fine. Armed with this knowledge, we ignore the first trial. 

%% LN Model for MFC 
% We now construct a LN model for the MFC for each trial. The following figure shows the filters backed out, together with the nonlinear residuals. 

% extract all MFC filters
if ~exist('K_MFC','var')
	K_MFC = zeros(width(PID),600);
	for j = 1:width(PID)
		temp = FindBestFilter(MFC_Control,MFC(j,:),[],'regmax=1;','regmin=1;','filter_length=799;','offset=100;');
		K_MFC(j,:) = temp(101:end-100);
	end
end

% make linear predictions
MFC_pred = NaN*MFC;
s =[];
for i = 2:5
	MFC_pred(i,:) = filter(K_MFC(i,:),1,MFC_Control);

	% fix trivial scaling
	x = MFC_pred(i,:); x = x(:);
	y = MFC(i,:); y = y(:);
	f = fit(x,y,'poly1');
	MFC_pred(i,:) = f(MFC_pred(i,:));
	s = [s f.p1];
end
s = mean(s);
K_MFC = K_MFC*s;

figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
for i = 2:5
	subplot(2,4,i-1), hold on
	plot(K_MFC(i-1,:));
	title(strcat('Trial#',mat2str(i)))
	xlabel('Lag (ms)')
	ylabel('Filter')

	subplot(2,4,4+i-1), hold on
	plot(MFC_pred(i,:),MFC(i,:),'k.')
	xlabel('Linear Prediction')
	ylabel('Data')
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% LN Model for PID 
% We now construct a LN model for the PID for each trial. The following figure shows the filters backed out, together with the nonlinear residuals. The model depends of course on the odour and the PID (the instrument).

% extract all PID filters
if ~exist('K_PID','var')
	K_PID = zeros(width(PID),600);
	for j = 1:width(PID)
		temp = FindBestFilter(MFC(j,:),PID(j,:),[],'regmax=1;','regmin=1;','filter_length=799;','offset=100;');
		K_PID(j,:) = temp(101:end-100);
	end
end

% make linear predictions
PID_pred = NaN*MFC;
s =[];
for i = 2:5
	PID_pred(i,:) = filter(K_PID(i,:),1,MFC(i,:));

	% fix trivial scaling
	x = PID_pred(i,:); x = x(:);
	y = PID(i,:); y = y(:);
	f = fit(x,y,'poly1');
	PID_pred(i,:) = f(PID_pred(i,:));
	s = [s f.p1];
end
s = mean(s);
K_PID = K_PID*s;


figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
for i = 2:5
	subplot(2,4,i-1), hold on
	plot(K_PID(i-1,:));
	title(strcat('Trial#',mat2str(i)))
	xlabel('Lag (ms)')
	ylabel('Filter')

	subplot(2,4,4+i-1), hold on
	plot(PID_pred(i,:),PID(i,:),'k.')
	xlabel('Linear Prediction')
	ylabel('Data')
end

PrettyFig;

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