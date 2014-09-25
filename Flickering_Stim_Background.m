% Flickering Stimulus on a Background
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

if ~exist('data')
	load('/local-data/DA-paper/flickering-stim-back/2014_09_23_srinivas_2ac_flickering_background.mat')
end

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end



%%
% Carlotta's data is almost entirely in the liner regime, and the effects of fast gain control, if any, are very small. Mahmut's "natural" stimuli are too poorly controlled, and the neuron response can't be fit by any model we have, making analysis challenging.

%%
% This document shows a middle path: we use a random flickering stimulus as in Carlotta's experiments, but now the flickering stimuli have varying heights, in addition to varying durations. Furthermore, they occur on top of a background, whose height can be varied too. 

%%
% The following figure shows one such paradigm, where pulses are presented on top of a background. The odor used is Ethyl Acetate. The different colors are different trials of the same experiment, and the black line is a fit to the background, showing good trial-to-trial reproducibility. 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
time = 1e-4:1e-4:1e-4*length(data(1).PID);
c = jet(width(data(1).PID));
for i = 1:width(data(1).PID)
	plot(time,data(1).PID(i,:),'Color',c(i,:))
end
set(gca,'XLim',[4 14])

% figure out the background level and show it
valve = ControlParadigm(1).Outputs(5,:);
[ons,offs]=ComputeOnsOffs(valve);
% for each segment, find the minimum
minPID = NaN*ons;
t_min = NaN*ons;
PID = mean2(data(1).PID);
for i = 1:length(ons)
	[minPID(i),loc]=min(PID(ons(i):offs(i)));
	t_min(i) = loc+ons(i);
end

c =mean(minPID)+.5*std(minPID);
t_min(minPID>c)=[];
minPID(minPID>c)=[];
% fit a line
ff = fit(t_min*1e-4,minPID,'poly1');
plot(time,ff(time),'k')
xlabel('Time (s)')
ylabel('PID (V)')

title(strcat('m=',oval(ff.p1,2),'V/s'))

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%% 
% The following figure shows the same stimulus on top of backgrounds of different amplitudes. Once again, we see that the background is mostly constant, and the pulse height is not really affected by the background. The lines are best fits to the estimated background concentrations. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = jet(length(data));
for i = 1:length(data)
	PID = mean(data(i).PID);
	plot(time,PID,'Color',c(i,:))

	% figure out the background level and show it
	valve = ControlParadigm(i).Outputs(5,:);
	[ons,offs]=ComputeOnsOffs(valve);
	% for each segment, find the minimum
	minPID = NaN*ons;
	t_min = NaN*ons;
	for j = 1:length(ons)
		[minPID(j),loc]=min(PID(ons(j):offs(j)));
		t_min(j) = loc+ons(j);
	end

	cutoff =mean(minPID)+.5*std(minPID);
	t_min(minPID>cutoff)=[];
	minPID(minPID>cutoff)=[];
	% fit a line
	ff = fit(t_min*1e-4,minPID,'poly1');
	plot(time,ff(time),'Color',c(i,:))

end

xlabel('Time (s)')
ylabel('PID (V)')
title('Ethyl acetate pulses on top of backgrounds')
set(gca,'XLim',[12 18],'box','on')
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end


