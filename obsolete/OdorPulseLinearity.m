% OdorPulseLinearity.m
% How the hell does odor delivery work? What does the PID
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

filename = '/local-data/DA-paper/odor-linearity/2014_11_04_ethyl_butyrate_pulse_addition_2sec_pulse_1.mat';

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

if ~exist('data')
	load(filename)
end

%%
% The file used is:
disp(filename)

%% How linearly do the Mass Flow Controllers behave?
% When we set our MFCs to a certain flow rate and dilution, do they actually go to that flow rate and dilution? The following figure shows how the actual flow rates (and dilution ratios) of the MFCs varies with setpoint. 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on

desired_dil_f = [];
actual_dil_f = [];
desired_dil_b = [];
actual_dil_b = [];
for i = 1:16
	if ~isempty(data(i).PID)
		for j = 1:width(data(i).PID)
			% do the foreground
			this_desired_dil = mean(ControlParadigm(i).Outputs(3,end-1000:end-100))/mean(ControlParadigm(i).Outputs(1,end-1000:end-100));
			this_actual_dil = mean(data(i).odor(j,end-1000:end-100))/mean(data(i).odor_dil(j,end-1000:end-100));
			desired_dil_f = [desired_dil_f this_desired_dil];
			actual_dil_f = [actual_dil_f this_actual_dil];

			% do the background
			this_desired_dil = mean(ControlParadigm(i).Outputs(4,end-1000:end-100))/mean(ControlParadigm(i).Outputs(2,end-1000:end-100));
			this_actual_dil = mean(data(i).bck(j,end-1000:end-100))/mean(data(i).bck_dil(j,end-1000:end-100));
			desired_dil_b = [desired_dil_b this_desired_dil];
			actual_dil_b = [actual_dil_b  this_actual_dil];
		end
	end
end

% correc for MFCs having differnet max flows
desired_dil_b = desired_dil_b*(2/5);
actual_dil_b = actual_dil_b*(2/5);
desired_dil_f = desired_dil_f*(1/2);
actual_dil_f= actual_dil_f*(1/2);

scatter(desired_dil_b,actual_dil_b,32,'b','filled')
xlabel('Desired Concentration')
ylabel('Actual Concentration')
title('Foreground Pulses')


subplot(1,2,2), hold on
scatter(desired_dil_f,actual_dil_f,32,'r','filled')
xlabel('Desired Concentration')
ylabel('Actual Concentration')
title('Background Pulses')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end


%%
% It looks like they behave very linearly. 

%% Linearity of PID response to individual pulses. 
% How linear is the PID? What does it scale with? The total flow of odor, or the dilution in the total flow? In the following figure, we plot all cases where only one pulse is presented, either as a "foreground" or as a "background". These labels are purely cosmetic, and don't determine anything about the delivery. 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on

foreground_flow = [];
background_flow = [];

foreground_dil = [];
background_dil = [];

foreground_pid = [];
background_pid = [];

total_flow = 2200; %ml/min, because only one valve turns on

for i = 1:16
	if ~isempty(data(i).PID)
		if any(ControlParadigm(i).Outputs(1:4,1) == 0)
			% only one flow through odor
			for j = 1:width(data(i).PID)
				% figure out which flow it is
				if ControlParadigm(i).Outputs(3,1)/ControlParadigm(i).Outputs(1,1) > 0

					foreground_flow = [foreground_flow 100*ControlParadigm(i).Outputs(3,1) ];
					foreground_dil = [foreground_dil foreground_flow(end)/total_flow];
					foreground_pid = [foreground_pid max(data(i).PID(j,:)) - mean(data(i).PID(j,1:1000))];

				else
					background_flow = [background_flow 40*ControlParadigm(i).Outputs(4,1)];
					background_dil = [background_dil background_flow(end)/total_flow];
					background_pid = [background_pid max(data(i).PID(j,:)) - mean(data(i).PID(j,1:1000))];
				end

			end
		end
		
	end
end

scatter(foreground_flow,foreground_pid,32,'r','filled')
scatter(background_flow,background_pid,32,'b','filled')
xlabel('Flow through odor (mL/min)')
ylabel('Peak PID (V)')


subplot(1,2,2), hold on
scatter(foreground_dil,foreground_pid,32,'r','filled')
scatter(background_dil,background_pid,32,'b','filled')
xlabel('Effective Concentration')
ylabel('Peak PID (V)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end


%% How do odor pulses add up?
% If we present two odor pulses from each of the valves together, at the same time, what does the PID measure? Do the odor pulses add up? In this figure, we consider combinations of two pulses, both of equal amplitude, that are added together, and compare them to cases when they were presented individually (on the right).

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

expected_PID = [];
actual_PID = [];
effective_dil = [];


total_flow = 2400; %ml/min, because only both valves turn on

for i = 1:16
	if ~isempty(data(i).PID)
		if ~any(ControlParadigm(i).Outputs(1:4,1) == 0)
			% both pulses on together
			for j = 1:width(data(i).PID)
				this_foreground_flow = 100*ControlParadigm(i).Outputs(3,1);
				this_background_flow = 40*ControlParadigm(i).Outputs(4,1);
				effective_dil = [effective_dil (this_background_flow+this_foreground_flow)/total_flow];

				expected_PID = [expected_PID (mean(foreground_pid(foreground_flow==this_foreground_flow)) + mean(background_pid(background_flow==this_background_flow)))];

				actual_PID = [actual_PID max(data(i).PID(j,:)) - mean(data(i).PID(j,1:1000))];


			end
		end
		
	end
end

subplot(1,2,1), hold on
scatter(expected_PID,actual_PID)
xlabel('Sum of Individual pulses (V)')
ylabel('Actual PID (V)')

% draw a line of unity
plot([0 max(expected_PID)],[0 max(expected_PID)],'k--')
plot([0 max(expected_PID)],[0 max(expected_PID)/2],'k--')

subplot(1,2,2), hold on
scatter(effective_dil,actual_PID,32,'r','filled')
xlabel('Effective Concentration')
ylabel('Actual PID (V)')

% also plot single pulses
single_pulses_dil = [foreground_dil background_dil];
single_pulses_pid = [foreground_pid background_pid];
scatter(single_pulses_dil,single_pulses_pid,32,'b','filled')

% fit a line to the double pulses
[ff,gof1] = fit(effective_dil(:),actual_PID(:),'poly1');
lh(1)=plot(unique(effective_dil),ff(unique(effective_dil)),'r');

[ff,gof2] = fit(single_pulses_dil(:),single_pulses_pid(:),'poly1');
lh(2)=plot(unique(single_pulses_dil),ff(unique(single_pulses_dil)),'b');

legend(lh,{oval(gof1.rsquare,2), oval(gof2.rsquare,2)} ,'Location','NorthWest')

title('Single Pulses (Blue), Double Pulses (red)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end


%% Are airspeeds at the outlet consistent?
% In the following figure, we plot airspeeds measured at the outlet when one of the valves opens (blue) compared to when both open (red)

load('/local-data/DA-paper/odor-linearity/2014_11_04_pulse_addition_airspeeds.mat')

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on


for i = 1:16
	if ~isempty(data(i).PID)
		if ~any(ControlParadigm(i).Outputs(1:4,1) == 0)
			% both pulses on together
			for j = 1:width(data(i).PID)
				plot(data(i).Airspeed(j,1:10:end),'r')

			end
		else
			% only one pulse on
			for j = 1:width(data(i).PID)
				plot(data(i).Airspeed(j,1:10:end),'b')

			end
		end
		
	end
end


xlabel('Time (ms)')
ylabel('-Airspeed (a.u.)')

PrettyFig;

if being_published
	snapnow;
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




