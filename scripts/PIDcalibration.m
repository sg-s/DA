% PIDcalibration.m
% 
% created by Srinivas Gorur-Shandilya at 9:52 , 30 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

%% PID Calibration
% In this document, we calibrate the PID and odour delivery system with known quantities of odorant, so that we can map out the relationships between flow rate, PID values, and gas-phase concentration. 

%% Strategy
% The way we do this is to introduce a known amount of ethyl acetate (here, 100uL) into the scintillation vial containing odour, and driving air through it at a fixed rate till we deplete all odour. Assuming that the overwhelming majority of the odourant in the vial is flushed out and measured, our problem is solved. 

%% 
% Ethyl acetate has a density of .902g/mL, so 100uL of the odourant weighs .0902 g. The molar mass of ethyl acetate is 88.11 g/mol, so 100uL of ethyl acetate contains approximately 1.023 milli moles of the compound. 

%%
% Similarly, we can calibrate the PID's sensitivity to 2-butanone. 2-butanone has a density of 0.8050 g/mL, so 100 uL of 2-butanone (that I used in this calibration) weighs 0.0805 g, which corresponds to 1.1164 milli moles of 2-butanone. 

%%
% In the following figure, we plot the raw traces of the PID vs. time, showing how changing the flow rate changes the dynamics of odourant depletion. As we expect, higher flow rates cause quicker depletion. In the second panel, we plot the PID value as a function of the air pushed through the bottle, which we obtain by multiplying the elapsed time by the flow rate. 

load('/local-data/PID-calibration/ethyl-acetate/calibration_data.mat')

% reformat x axis to a more sensical format
min_pid = Inf;
for i = 1:length(alldata)
	t = alldata(i).t;
	t = t - t(1);
	t = datevec(t);
	alldata(i).t = t(:,5)*60 + t(:,6);
	min_pid = min([min_pid; alldata(i).PID]);
end

% remove minimum
for i = 1:length(alldata)
	alldata(i).PID = alldata(i).PID - min_pid;
end

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(2,1,1), hold on
c = parula(length(alldata)+1);
L = {};

for i = 1:length(alldata)
	t = alldata(i).t;
	plot(t-t(1),alldata(i).PID,'Color',c(i,:))
	L{i} = [mat2str(alldata(i).flowrate) 'mL/min'];

end
legend(L)
xlabel('Time (s)')
ylabel('PID (V)')

subplot(2,1,2), hold on
cum_air_flow = NaN(length(alldata),1);
for i = 1:length(alldata)
	t = alldata(i).t;
	t = t-t(1);
	t = t*alldata(i).flowrate/60;
	cum_air_flow(i) = max(t);
	plot(t,alldata(i).PID,'Color',c(i,:))
end
xlabel('Cumulative Air Flow (mL)')
ylabel('PID (V)')

ax = axes;
ax.Position = [0.6558    0.2776    0.2110    0.1657];
plot([alldata.flowrate],cum_air_flow,'k+')
set(gca,'XLim',[0 250],'YLim',[0 3100])
xlabel('Flow rate (mL/min)')
ylabel(['Air required' char(10) 'to deplete (mL)'])

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% Why do we need different cumulative air flows to completely deplete the odour? That's a bit weird -- I would expect that you need the same amount of air. 

%%
% We now plot the peak, mean, and integrated PID value as a function of the flow rate. We see that the peak and mean PID increase with flow rate, as expected. However, we see something weird with the integrated PID signal -- it decreases with flow rate. This definitely breaks my assumption that all the odor is measured. Why does it decrease? One possibility is that as flow rates increase, the PID captures a smaller and smaller fraction of the total odor flux (because the PID has a fixed suction rate of 750mL/min). This model predicts the data well (red curve), but the best-fit parameter in this model for the total flow rate is ten times lower than in reality. So I fall back to ignoring why this happens, and simply fitting the data with a smoothing function, and use this function to guess what integrated PID value we expect to see, given a certain (constant) flow rate. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
clear l 
for i = 1:length(alldata)
	l(1) = plot(alldata(i).flowrate/60,mean(alldata(i).PID),'k+');
	l(2) = plot(alldata(i).flowrate/60,max(alldata(i).PID),'r+');
end
xlabel('Flow Rate (mL/s)')
ylabel('PID (V)')
legend(l,{'Mean','Peak'},'Location','northwest')

subplot(1,2,2), hold on
x = []; y = [];
for i = 1:length(alldata)
	x = [alldata(i).flowrate/60 x]; % units now in mL/s
	y = [mean(diff(alldata(i).t))*sum((alldata(i).PID)) y]; % this is the time integrated PID signal
end
ff = fit(x(:),y(:),'smoothingspline');
xx = linspace(0.5,max(x),100);
clear l
l(1) = plot(xx,ff(xx),'b');
plot(x,y,'k+')
ylabel('\int PID (V s)')
xlabel('Flow Rate (mL/s)')
set(gca,'XLim',[0 max(x)],'YLim',[0 1.1*max(y)])

% compute the fraction captured 
ft = fittype('A*(12.5/(B+x))');
ff2 = fit(x(:),y(:),ft,'StartPoint',[190 3],'Lower',[0 0]);
l(2) = plot(xx,ff2(xx),'r');

legend(l,{'Smoothing spline','Model'},'Location','southwest')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now use this expected value of integrated PID signal to back-calculate the odour flux as a function of time for these long pulses. The bottom panel also shows the cumulative odour delivered as a function of time, for the different flow rates. The dashed line is the expected # of moles delivered, based on the volume of odourant put into the vial in the first place. 

peak_flux = NaN(length(alldata),1);
asym_flux = NaN(length(alldata),1);
peak_PID = NaN(length(alldata),1);
asym_PID = NaN(length(alldata),1);

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1), hold on
for i = 1:length(alldata)
	odour_flux = mean(diff(alldata(i).t))*alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate/60);
	peak_flux(i) = max(odour_flux);
	at = find((alldata(i).t - alldata(i).t(1)) > 400,1,'first');
	asym_flux(i) = odour_flux(at);
	asym_PID(i) = alldata(i).PID(at);
	plot(alldata(i).t,odour_flux,'Color',c(i,:))
end
legend(L)
ylabel('Odour flux (mol/s)')

subplot(2,1,2), hold on
odour_flux = [];
for i = 1:length(alldata)
	plot(alldata(i).t,mean(diff(alldata(i).t))*cumsum(alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate/60)),'Color',c(i,:))
	x = alldata(i).t;
	y = mean(diff(alldata(i).t))*cumsum(alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate/60));
	x = x(:); y = y(:);
	temp = fit(x,y,'poly1');
	odour_flux(i) = temp.p1;
end
ylabel('cumulative odour (mol)')
xlabel('Time (s)')
plot([0 3500],[1.023e-3 1.023e-3],'k--')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% 
% The fact that all the cumulative sums end up close to the dashed line means that our estimation techniques are somewhat OK. But what we really want is a map from flow rates to odour flux, so that we can use this widely. This also has the advantage of being independent of PID sensitivity, and we can use PID values to correct small details of the signal, instead of the overall scale, which we ant to be accurate and instrument-independent. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
x = [alldata.flowrate]; x = x(:);
fit_flow_rate_2_odour_flux = fit([0; -10; -100; x],1e6*[0; 0; 0; odour_flux(:)],'smoothingspline');
plot(1:200,fit_flow_rate_2_odour_flux(1:200),'r')
plot(x,1e6*odour_flux,'k+')
prettyFig()
ylabel('Odour Flux (\mu mol /s)')
xlabel('Mean Flow Rate (mL/min)')

if being_published
	snapnow
	delete(gcf)
end


%% Linearity of Response
% In this section, I check for the linearity of the PID response by plotting calculated odor flux as a function of the measured PID for the two odors I calibrated. You can see that the PID is very sensitive to 2-butanone, but less sensitive to ethyl-acetate. Both calibration curves look very linear. Fits are forced to go through the origin. 

peak_PID = cellfun(@max,{alldata.PID});
lf = fit(peak_PID(:),peak_flux,'poly1','Upper',[Inf 0],'Lower',[-Inf 0]);
x = linspace(0,max(peak_PID),10);
clear l


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(peak_PID,peak_flux,'r+');
plot(asym_PID,asym_flux,'ro');
plot(x,lf(x),'r');
r2_2ac = rsquare(lf(peak_flux),peak_PID);

ylabel('odorant flux (mol/s)')
xlabel('PID response (V)')

% now also do the 2-butanone calibration

load('/local-data/PID-calibration/2-butanone/alldata.mat','alldata')

% reformat x axis to a more sensical format
min_pid = Inf;
for i = 1:length(alldata)
	t = alldata(i).t;
	t = t - t(1);
	t = datevec(t);
	alldata(i).t = t(:,4)*60*60 + t(:,5)*60 + t(:,6);
	min_pid = min([min_pid; alldata(i).PID]);
end

% remove minimum
for i = 1:length(alldata)
	alldata(i).PID = alldata(i).PID - min_pid;
end

x = []; y = [];
for i = 1:length(alldata)
	x = [alldata(i).flowrate/60 x]; % units now in mL/s
	y = [mean(diff(alldata(i).t))*sum((alldata(i).PID)) y]; % this is the time integrated PID signal
end
ff = fit(x(:),y(:),'smoothingspline','SmoothingParam',1-1e-6);


peak_flux = NaN(length(alldata),1);
asym_flux = NaN(length(alldata),1);
asym_PID = NaN(length(alldata),1);

for i = 1:length(alldata)
	odour_flux = mean(diff(alldata(i).t))*alldata(i).PID*(1.1164e-3)/ff(alldata(i).flowrate/60);
	peak_flux(i) = max(odour_flux);
	at = find((alldata(i).t - alldata(i).t(1)) > 400,1,'first');
	asym_flux(i) = odour_flux(at);
	asym_PID(i) = alldata(i).PID(at);
end
peak_PID = cellfun(@max,{alldata.PID});
lf = fit(peak_PID(:),peak_flux,'poly1','Upper',[Inf 0],'Lower',[-Inf 0]);
x = linspace(0,max(peak_PID),10);

plot(peak_PID,peak_flux,'b+');
plot(asym_PID,asym_flux,'bo');
plot(x,lf(x),'b');
r2_2but = rsquare(lf(peak_flux),peak_PID);

clear l L
l(1) = plot(NaN,NaN,'r');
l(2) = plot(NaN,NaN,'b');

L{1} = ['ethyl acetate, r^2 = ' oval(r2_2ac)];
L{2} = ['2-butanone, r^2 = ' oval(r2_2but)];
legend(l,L,'Location','northeast');

set(gca,'YLim',[0 1.4e-5],'XLim',[0  8])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Version Info

pFooter;
