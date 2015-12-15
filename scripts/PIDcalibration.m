% PIDcalibration.m
% 
% created by Srinivas Gorur-Shandilya at 9:52 , 30 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic

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


figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(2,1,1), hold on
c = parula(length(alldata)+1);
L = {};

for i = 1:length(alldata)
	t = alldata(i).t;
	plot(t-t(1),alldata(i).PID - min_pid,'Color',c(i,:))
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
	plot(t,alldata(i).PID-min_pid,'Color',c(i,:))
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
	l(1) = plot(alldata(i).flowrate,mean(alldata(i).PID),'k+');
	l(2) = plot(alldata(i).flowrate,max(alldata(i).PID),'r+');
end
xlabel('Flow Rate (mL/min)')
ylabel('PID (V)')
legend(l,{'Mean','Peak'},'Location','northwest')

subplot(1,2,2), hold on
x = []; y = [];
for i = 1:length(alldata)
	x = [alldata(i).flowrate x];
	y = [mean(diff(alldata(i).t))*sum((alldata(i).PID)) y];
end
ff = fit(x(:),y(:),'smoothingspline');
plot(1:max(x),ff(1:max(x)),'r')
plot(x,y,'k+')
ylabel('\int PID (V s)')
xlabel('Flow Rate (mL/min)')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now use this expected value of integrated PID signal to back-calculate the odour flux as a function of time for these long pulses. The bottom panel also shows the cumulative odour delivered as a function of time, for the different flow rates. The dashed line is the expected # of moles delivered, based on the volume of odourant put into the vial in the first place. 

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1), hold on
for i = 1:length(alldata)
	plot(alldata(i).t,mean(diff(alldata(i).t))*alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate),'Color',c(i,:))
end
legend(L)
ylabel('Odour flux (mol/s)')

subplot(2,1,2), hold on
odour_flux = [];
for i = 1:length(alldata)
	plot(alldata(i).t,mean(diff(alldata(i).t))*cumsum(alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate)),'Color',c(i,:))
	x = alldata(i).t;
	y = mean(diff(alldata(i).t))*cumsum(alldata(i).PID*(1.023e-3)/ff(alldata(i).flowrate));
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

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
