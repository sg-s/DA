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
% In the following figure, we plot the raw traces of the PID vs. time, showing how changing the flow rate changes the dynamics of odourant depletion. As we expect, higher flow rates cause quicker depletion. In the second panel, we plot the PID value as a function of the air pushed through the bottle, which we obtain by multiplying the elapsed time by the flow rate. 

load('/local-data/PID-calibration/2ac/calibration_data.mat')

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
for i = 1:length(alldata)
	t = alldata(i).t;
	t = t-t(1);
	t = t*alldata(i).flowrate/60;
	plot(t,alldata(i).PID,'Color',c(i,:))
end
xlabel('Cumulative Air Flow (mL)')
ylabel('PID (V)')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%%
% We now plot the peak, mean, and integrated PID value as a function of the flow rate. We see that the peak and mean PID increase with flow rate, as expected. However, we see something weird with the integrated PID signal -- it decreases with flow rate. Why? I don't know. I fit a smooth function to the data, and use this function to guess what integrated PID value we expect to see, given a certain (constant) flow rate. 

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
plot(linspace(0.5,max(x),100),ff(linspace(0.5,max(x),100)),'r')
plot(x,y,'k+')
ylabel('\int PID (V s)')
xlabel('Flow Rate (mL/s)')
set(gca,'XLim',[0 max(x)],'YLim',[0 1.1*max(y)])

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now use this expected value of integrated PID signal to back-calculate the odour flux as a function of time for these long pulses. The bottom panel also shows the cumulative odour delivered as a function of time, for the different flow rates. The dashed line is the expected # of moles delivered, based on the volume of odourant put into the vial in the first place. 

peak_flux = NaN(length(alldata),1);

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1), hold on
for i = 1:length(alldata)
	odour_flux = mean(diff(alldata(i).t))*alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate/60);
	peak_flux(i) = max(odour_flux);
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
% In this section, we check for the linearity of the PID response by plotting the PIF value as a function of the odor flux. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
y = cellfun(@max,{alldata.PID});
plot(peak_flux,y,'k+');
lf = fit(peak_flux,y(:),'poly1');
l = plot(sort(peak_flux),lf(sort(peak_flux)),'r');
L = ['r^2 = ' oval(rsquare(y,lf(peak_flux)))];
legend(l,L,'Location','southeast')
xlabel('ethyl acetate flux (mol/s)')
ylabel('PID response (V)')
set(gca,'XLim',[0 1.3e-5],'YLim',[0  2.2])
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info

pFooter;
