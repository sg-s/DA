% LargeVarianceFlickering.m
% 
% created by Srinivas Gorur-Shandilya at 3:30 , 19 January 2015. Contact me at http://srinivas.gs/contact/
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

%% Response of ORNs to flickering odor stimuli with large variances
% From previous experiments, we see that the response of ORNs to largely varying stimuli that mimics the natural odor plumes is particularly interesting: in that no model we have can precisely account for the data, and in that, from other results, we expect to see a large variation of gain of the ORNs to these stimuli. 

%%
% In this document, we generate odor stimuli flickers over a large range, like the "natural" stimuli, but never goes to zero, so that the neuron should never silence (allowing us to accurately follow its response). 

%%
% The odor stimulus looks like this:


load('/local-data/DA-paper/natural-flickering/2015_01_19_large_variance_flickering.mat')
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
haz_data = find_data(data);
time = 1e-4*(1:length(data(haz_data(1)).PID));
for i = 1:length(haz_data)-1
	plot(time,mean2(data(haz_data(i)).PID))
end

set(gca,'XLim',[40 60])
xlabel('Time (s)')
ylabel('PID (V)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end
