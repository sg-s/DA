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
	plot(MFC_pred(i,1000:10:end-1000),MFC(i,1000:10:end-1000),'k.')
	xlabel('Linear Prediction')
	ylabel('Data')
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% What causes those bizzare excursions? Well, what's happening there is that when very low values are written to the MFC, it shuts down to zero, and when it does, it stays there. (a quick inspection by eye shows that it stays shut for ~500ms). So this is a problem. We work around this by never writing a value lower than the maximum flow/turndown ratio to the MFC. 

load('/local-data/DA-paper/mean-shifted-gaussians/2015_02_26_odour_test_3.mat')
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
% remove PID baseline
PID = PID - mean(data(1).PID(6,:));

% extract all MFC filters
if ~exist('K_MFC3','var')
	K_MFC3 = zeros(width(PID),600);
	for j = 1:width(PID)
		temp = FindBestFilter(MFC_Control,MFC(j,:),[],'regmax=1;','regmin=1;','filter_length=799;','offset=100;');
		K_MFC3(j,:) = temp(101:end-100);
	end
end

% make linear predictions
MFC_pred = NaN*MFC;
s =[];
for i = 2:5
	MFC_pred(i,:) = filter(K_MFC3(i,:),1,MFC_Control);

	% fix trivial scaling
	x = MFC_pred(i,:); x = x(:);
	y = MFC(i,:); y = y(:);
	f = fit(x,y,'poly1');
	MFC_pred(i,:) = f(MFC_pred(i,:));
	s = [s f.p1];
end
s = mean(s);
K_MFC3 = K_MFC3*s;

figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
for i = 2:5
	subplot(2,4,i-1), hold on
	plot(K_MFC3(i-1,:));
	title(strcat('Trial#',mat2str(i)))
	xlabel('Lag (ms)')
	ylabel('Filter')

	subplot(2,4,4+i-1), hold on
	plot(MFC_pred(i,1000:10:end-1000),MFC(i,1000:10:end-1000),'k.')
	xlabel('Linear Prediction')
	ylabel('Data')
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% That's pretty nice. We have significantly reduced the number of weird excursions. Now, we repeat the experiment with odour, and proceed to build a LN model for the PID signal. 


%% LN Model for PID 
% We now construct a LN model for the PID for each trial. The following figure shows the filters backed out, together with the nonlinear residuals. The model depends of course on the odour and the PID (the instrument).

% convert MFC flows into a dilution
dil = MFC*100; % 1V = 100mL/min
dil = dil./(dil + 2000);

% extract all PID filters
if ~exist('K_PID','var')
	K_PID = zeros(width(PID),600);
	for j = 1:width(PID)
		temp = FindBestFilter(dil(j,:),PID(j,:),[],'regmax=1;','regmin=1;','filter_length=799;','offset=100;');
		K_PID(j,:) = temp(101:end-100);
	end
end

% make linear predictions
PID_pred = NaN*dil;
s =[];
for i = 2:5
	PID_pred(i,:) = filter(K_PID(i,:),1,dil(i,:));

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
	plot(PID_pred(i,1000:10:end-1000),PID(i,1000:10:end-1000),'k.')
	xlabel('Linear Prediction')
	ylabel('Data')
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% 
% We can now combine everything and make a combined linear model for the odour delivery system. The following figure shows the plot of this prediction to each of the four trials we consider. 



K_MFC = mean2(K_MFC3(2:end,:));
K_PID = mean2(K_PID(2:end,:));
MFC_Scale = 100;
Total_Flow = 2000;

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
PID_pred = DeliverySystemModel(MFC_Control);
clear a
for i = 2:5
	a(i) = autoplot(4,i-1); hold on
	
	plot(PID_pred(1000:10:end-1000),PID(i,1000:10:end-1000),'k.')
	title(strcat('Trial#',mat2str(i)))
	xlabel('Linear Model')
	ylabel('Data')
end

if ~exist('p_LN')
	d.stimulus = PID_pred(1000:10:end-1000);
	d.response = mean2(PID(2:5,1000:10:end-1000));
	p_LN = fit(d.stimulus(:),d.response(:),'poly3');
	for i = 2:5
		plot(a(i),0:1e-3:max(PID_pred),p_LN(0:1e-3:max(PID_pred)),'r')
	end
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we plot the output of the combined model. 

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
PID_pred = DeliverySystemModel_LN(MFC_Control);
for i = 2:5
	a(i) = autoplot(4,i-1); hold on
	plot(time,PID(i,:),'k')
	l = plot(time,PID_pred,'r');
	legend(l,strcat('r^2=',oval(rsquare(PID_pred,PID(i,:)))))
	set(gca,'XLim',[10 20])
	xlabel('Time (s)')
	ylabel('PID (V)')
	title(strcat('Trial#',mat2str(i)))
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% So our model of the delivery system is pretty neat (for this particular odour, dose, etc.)

%% Distribution Analysis
% In this section, we determine if we can fine-tune the shape of the distribution. First, what does it look like? The following figure shows the MFC and PID distributions. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
for i = 2:5
	[~,~,x,y] = splinehist(MFC(i,:));
	plot(x,y)
end
ylabel('pdf')
xlabel('MFC (V)')

subplot(1,2,2), hold on
for i = 2:5
	[~,~,x,y] = splinehist(PID(i,:));
	plot(x,y)
end
xlabel('PID (V)')
ylabel('pdf')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% 
% Now, we build a model for the Odour Delivery System, and then feed it with some control signal to get the simulated PID output. We also build another wrapper model that accepts a parametrised distribution and predicts the distribution of PID values.

%%
% How well does this model do? In the following figure, we feed the model with the actual MFC Control signal, and compare the distribution of PID values with what the model predicts. 


[cy,cx] = hist(MFC_Control,50);
cy(1:5) = 0;
cy = interp1(cx,cy,0:1e-3:5);
cx = 0:1e-3:5;
cy(isnan(cy)) = 0;
cy=cy/max(cy);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(cx,cy,'k')
xlabel('MFC Control (V)')
ylabel('pdf')

subplot(1,2,2), hold on
[py,px]=hist(mean2(PID(2:5,:)),50);
py = interp1(px,py,0:1e-3:5);
px = 0:1e-3:5;
py(isnan(py)) = 0;
py=py/max(py);
plot(px,py,'k')
py_hat = BestDistribution([],MFC_Control);
plot(px,py_hat,'r')
xlabel('PID (V)')
ylabel('pdf')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% That's pretty solid. Now, we numerically solve the problem of choosing the best distribution of MFC Control Signals to get a target PID distribution. Here, let's say we want a target PID distribution that looks like a Gaussian. 


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