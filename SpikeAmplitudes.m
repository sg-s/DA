% SpikeAmplitudes.m
% 
% created by Srinivas Gorur-Shandilya at 12:58 , 30 January 2015. Contact me at http://srinivas.gs/contact/
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


%% Understanding the factors controlling Spike Amplitudes
% In this document, we analyse extracellular recordings of Drosophila ORNs. These ORNs (typically in the ab3 sensillum) have two ORNs, called the A and B neuron. Usually, the A neuron has a bigger amplitude. However, the observed extraceullar spike amplitude has been known to reduce when the ORN is strongly activated, a phenomenon called pinching (which is AFAIK completely unexplained).

%%
% In the following figure, a ab3 sensillum is being recorded from, while presenting a fluctuating ethyl acetate odor. This fly has not been (intentionally) exposed to odor before; this is the first trial of the experiment. 

uiopen('/local-data/DA-paper/large-variance-flicker/2015-01-28-CS-ab3-2-EA.mat-s=0.3-Trial:1.fig',1)

PrettyFig('plw=1;');
if being_published
	snapnow
	delete(gcf)
end

%%
% % The red dots indicate the minima of the A spikes, and the blue dots the minima of the B spikes. What do subsequent trials look like? 

uiopen('/local-data/DA-paper/large-variance-flicker/2015-01-28-CS-ab3-2-EA.mat-s=0.3-Trial:3.fig',1)

PrettyFig('plw=1;');
if being_published
	snapnow
	delete(gcf)
end

uiopen('/local-data/DA-paper/large-variance-flicker/2015-01-28-CS-ab3-2-EA.mat-s=0.3-Trial:5.fig',1)

PrettyFig('plw=1;');
if being_published
	snapnow
	delete(gcf)
end

%%
% From the first trial, it looks as though the minima of the A spikes is drifting up with time, while the minima of the B spikes are more or less invariant. Let's plot this:

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
A = full(spikes(4).A(1,:));
B = full(spikes(4).B(1,:));

V = filter_trace(data(4).voltage(1,:));
time = 1e-4*(1:length(V));

Vmin_A = V(find(A));
Vmin_B = V(find(B));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
clear l
l(1)=plot(time(A==1),Vmin_A,'r');
l(2)=plot(time(B==1),Vmin_B,'b');
ylabel('Spike Minimum (\muV)')
xlabel('Time (s)')

legend({strcat('CV=',oval(std(abs(Vmin_A))/mean(abs(Vmin_A)),2)),strcat('CV=',oval(std(abs(Vmin_B))/mean(abs(Vmin_B)),2))})

PrettyFig('plw=1;');
if being_published
	snapnow
	delete(gcf)
end

%%
% The conventional wisdom on "pinching" states that spike amplitudes shrink when the ORN is active. Let's plot the spike amplitude vs. the firing rate. 

[fA,tA] = spiketimes2f(A,time);
Vmin_A_t = interp1(time(A==1),Vmin_A,tA);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(tA,fA,'r')
ylabel('Firing Rate (Hz)')

subplot(2,1,2), hold on
xlabel('Time (s)')
ylabel('Spike Minimum (\muV)')
plot(tA,Vmin_A_t,'r')
PrettyFig('plw=1;');

if being_published
	snapnow
	delete(gcf)
end

%%
% The two traces are somewhat correlated with a coefficient of determination of 

disp(rsquare(fA,Vmin_A_t))

%% 
% How does the minimum of the spike correlate with the spike amplitude (preceding maxima-minima)? 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
spike_amplitudes_A = spikes(4).amplitudes_A(1,A==1);
spike_amplitudes_B = spikes(4).amplitudes_B(1,B==1);
plot(Vmin_A,spike_amplitudes_A,'k.')
xlabel('Voltage Minima (\muV)')
ylabel('Spike Amplitude (\muV)')
title('ab3A spikes')
subplot(1,2,2), hold on
xlabel('Voltage Minima (\muV)')
title('ab3B spikes')
plot(Vmin_B,spike_amplitudes_B,'k.')
PrettyFig('EqualiseY=1;','EqualiseX=1;');


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
