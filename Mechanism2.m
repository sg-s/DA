% Mechanism2.m
% in this document, we attempt to understand the mechanism behind gain adaptation
% 
% created by Srinivas Gorur-Shandilya at 1:15 , 06 April 2015. Contact me at http://srinivas.gs/contact/
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
tic


%% What is the mechanism behind gain adaptation?
% Two (non-exclusive) possibilities are that the gain is controlled at the receptor level, and that the gain is controlled at the firing machinery. To disambiguate the two, we record the neuron's responses to a flickering odour stimulus. We then repeat the experiment, but also drive the neuron optically through light activating ReaChR channels. 

%% LEDs can deliver lots of power through the objective
% We deliver light through the objective, as shown below:
%
% <</Users/sigbhu/code/da/images/led.jpg>>
%

%%
% The following figure shows the light levels measured at the fly position vs the voltage driving the LED circuit. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
load('LED_calib.mat')
plot(LED.x,LED.y,'k+')
xlabel('Driving Voltage (V)')
ylabel('Power @ 627nm (mW)')
clear p
p.A= 2.7288;
p.k= 1.9863;
p.n= 2.8827;
l=plot(LED.x,hill(LED.x,p),'r');
legend(l,'Hill Fit','location','southeast')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% This is comparable to the light delivered through the microscope, and is much brighter than the LED not through the objective:

LightLevels = [287 315 315; 297 323 308; 294 313 316; 304 309 313; 246 241 241 ; 127 162 165; 771 727 725; 1280 1320 1250; 2210 2060 2070];

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(mean(LightLevels'),std(LightLevels'),'kx')
set(gca,'XTickLabel',{'LED:5V','LED:4V','LED:3V','LED:2V','LED:1V','M: 2V','M: 3V','M: 3.5V','M: 4V'})
set(gca,'XTick',[1:9],'XTickLabelRotation',45)
ylabel('Power @ 591nm (\muW)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% LED at maximum power can barely elicit spikes
% The following figure shows the responses of the neuron to the LED at maximum power delivered through the objective:

load('/local-data/DA-paper/reachr/2015_04_08_RR_F2_ab3_2_EA.mat')
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
time = 1e-4*(1:length(data(6).PID));
fA=spiketimes2f(spikes(6).A([2 3 5 6],:),time,1e-3);
tA = 1e-3*(1:length(fA));
plot(tA,mean2(fA),'r')
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Direct LED illumination is very bright
% By mounting the LED on a micromanipulator, we can introduce it under the objective, a few millimetres above the fly. This leads to a lot of light being dumped on the fly. We measure ~35mW of power at the fly location, a more than 10x improvement. 
%
% <</Users/sigbhu/code/da/images/led2.jpg>>
%





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
