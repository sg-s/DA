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

%% Direct LED illumination is very effective in eliciting spikes
% By sticking the LED a few mm above the fly, we can activate the neuron strongly: 

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,2,1), hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_1_EA.mat')
time = 1e-4*(1:length(data(6).PID));
[fA,tA]=spiketimes2f(spikes(5).A,time,1e-3);
plot(tA,fA)
[fA,tA]=spiketimes2f(spikes(6).A,time,1e-3);
plot(tA,fA)
legend({'3V','4V'})
xlabel('Time (s)')
title('Neuron 1')
ylabel('Firing Rate (Hz)')

subplot(2,2,2), hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_2_EA_2.mat')
[fA1,tA]=spiketimes2f(spikes(3).A,time,1e-3);
[fA2,tA]=spiketimes2f(spikes(4).A,time,1e-3);
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_2_EA.mat')
temp=spiketimes2f(spikes(3).A,time,1e-3);
fA1 = [fA1 temp];
temp=spiketimes2f(spikes(4).A,time,1e-3);
fA2 = [fA2 temp];
plot(tA,mean2(fA1))
plot(tA,mean2(fA2))
[fA,tA]=spiketimes2f(spikes(5).A,time,1e-3);
plot(tA,fA)
[fA,tA]=spiketimes2f(spikes(6).A,time,1e-3);
plot(tA,mean2(fA))
legend({'1V','2V','3V','4V'})
xlabel('Time (s)')
title('Neuron 2')
ylabel('Firing Rate (Hz)')

subplot(2,2,3), hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_3_EA.mat')
[fA,tA]=spiketimes2f(spikes(7).A,time,1e-3);
plot(tA,fA)
legend({'2.2V'})
xlabel('Time (s)')
title('Neuron 3')
ylabel('Firing Rate (Hz)')

subplot(2,2,4), hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_4_EA.mat')
[fA,tA]=spiketimes2f(spikes(6).A,time,1e-3);
plot(tA,mean2(fA))
legend({'4V'})
xlabel('Time (s)')
title('Neuron 4')
ylabel('Firing Rate (Hz)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% However, moving the LED in and out is very tricky and hard to do.

%% Long Pulses of light cause ORN to oscillate 
% Now we present the light in the same way (around 30mM of light just above the fly) and present long pulses of light to the fly. As the following figure shows, instead of responding with a  steady firing rate, the neuron oscillates:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_2_EA.mat')
time = 1e-4*(1:length(data(9).PID));
[fA,tA]=spiketimes2f(spikes(9).A,time,1e-3);
plot(tA,fA)
legend({'2V'})
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% If we fit a sine wave to this oscillation, we get the time period of the oscillation to be: 


x = tA(tA>10&tA<50);
y = fA(tA>10&tA<50);
y = y-mean(y);
y = y/(2*std(y));
cf=fit(x(:),y(:),'sin1');
disp(cf.b1/(2*pi))

%%
% which is suspiciously close to 1s.

%% Driving the neuron with light and odour decreases spike amplitude but leaves firing rate unchanged. 
% In the following figure, this neuron responds strongly to light activation:

figure('outerposition',[0 0 700 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F2_ab3_2_EA_2.mat')
time = 1e-4*(1:length(data(3).PID));
[fA,tA]=spiketimes2f(spikes(3).A,time,1e-3);
plot(tA,fA)
[fA,tA]=spiketimes2f(spikes(4).A,time,1e-3);
plot(tA,fA)
legend({'1V','2V'})
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% The following plot shows the responses to odour flicker with and without supplemental light activation:


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
time = 1e-4*(1:length(data(12).PID)); 
ax(1)=subplot(1,3,1:2); hold on
ax(2)=subplot(1,3,3); hold on
[fA0,tA]=spiketimes2f(spikes(12).A,time,1e-3);
l(1) = plot(ax(1),tA,mean2(fA0));
[y,x] = hist(mean2(fA0),50);
plot(ax(2),y,x)
[fA1,tA]=spiketimes2f(spikes(13).A(2:end,:),time,1e-3);
l(2) = plot(ax(1),tA,mean2(fA1));
[y,x] = hist(mean2(fA1),50);
plot(ax(2),y,x)
xlabel(ax(1),'Time (s)')
ylabel(ax(1),'Firing Rate (Hz)')
legend(l,{'Odour + 0V Light','Odour + 1V Light'})
set(ax(1),'XLim',[20 60])
PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% In the following plot, we show how the gain changes by plotting the two responses against each other:

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
x = mean2(fA0);y= mean2(fA1);
ss = 20;
plot([0 100],[0 100],'k--')
plot(x(1:ss:end),y(1:ss:end),'k.')
xlabel('Response to odour flicker (Hz)')
ylabel('Response to odour flicker+light(Hz)')
set(gca,'XLim',[0 max(y)*1.1],'YLim',[0 max(y)*1.1])

cf = fit(x(:),y(:),'poly1');
l=plot(sort(x),cf(sort(x)),'r');
legend(l,strcat('Rel. gain = ',oval(cf.p1)))

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% ORN pinches even at low firing rates
% In the following figure, we show that the ORN pinches (spike amplitude becomes very small) even at very modest firing rates (here, the firing rate never exceeds 50Hz)
%
% <</Users/sigbhu/code/da/images/low-f-pinch.png>>
%





%% Feeding flies with 2mM ATR greatly increases their sensitivity to light
% In all previous experiments, flies were fed with 400uL of 400uM ATR. In these experiments, they were fed with 400uL of 2mM ATR. The neurons were so sensitive to light that we could activate them with room lights! In the following figure, we moved the LEDs far away from the fly, and very weakly turned them on. Even then, we can elicit strong responses:
%
% <</Users/sigbhu/code/da/images/led_2mM.jpg>>
%


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
load('/local-data/DA-paper/reachr/2015_04_10_RR_F3_ab3_1_EA_2mM.mat')
time = 1e-4*(1:length(data(3).PID));
[fA,tA]=spiketimes2f(spikes(3).A,time,1e-3);
plot(tA,fA)
[fA,tA]=spiketimes2f(spikes(8).A,time,1e-3);
plot(tA,fA)
legend({'1V','1.5V'})
xlabel('Time (s)')
title('Neuron from fly fed with 2mM ATR')
ylabel('Firing Rate (Hz)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% They respond just as strong as in the previous case, which means upping the ATR concentration led to a 5000x improvement in sensitivity! 

%% Flies fed with 2mM ATR are more sensitive to odour
% While attempting to perform the experiments, we observed that odour stimuli that elicited a nice response in earlier experiemnts quickly saturated the neuron, and were unusable. We traced this to the fact that neurons in flies fed with 2mM ATR are much, much more sensitive to odour: 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,2,4), hold on

load('/local-data/DA-paper/reachr/2015_04_10_RR_F3_ab3_5_EA_2mM.mat')
time = 1e-4*(1:length(data(11).PID));
[fA,tA]=spiketimes2f(spikes(11).A,time,1e-3);
plot(tA,mean2(fA))
xlabel('Time (s)')
set(gca,'YLim',[0 150])
subplot(2,2,2), hold on
title('fly fed with 2mM ATR')
plot(time,mean2(data(11).PID))

set(gca,'YLim',[0 2])


subplot(2,2,1), hold on
title('fly fed with 400uM ATR')
load('/local-data/DA-paper/reachr/2015_04_03_ReaChR_F3_ab3_3_EA.mat')
time = 1e-4*(1:length(data(2).PID));
plot(time,mean2(data(2).PID))
set(gca,'XLim',[0 10])
set(gca,'YLim',[0 2])
ylabel('PID (V)')

[fA,tA]=spiketimes2f(spikes(2).A,time,1e-3);
subplot(2,2,3), hold on
plot(tA,mean2(fA))
set(gca,'XLim',[0 10])
set(gca,'YLim',[0 150])
ylabel('ab3A response (Hz)')

PrettyFig();

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
