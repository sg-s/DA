% 
% 
% created by Srinivas Gorur-Shandilya at 1:53 , 02 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Fast Gain Control w/o receptors?
% In this document, we try to see if we can see evidence of fast gain control in our data where we stimulate ORNs with light, but no odour. This document uses the new ORNData class, and all plots are generated using overloaded plot functions defined within this class. 

pFooter;

% load the data
load('/Users/sigbhu/code/da/data/Chrimson_data.mat') 

% specify that we only want to use data from 5 to 55s
temp = zeros(length(light_data.stimulus),1);
temp(5e3:55e3) = 1;
light_data.use_this_segment = logical(temp);

% back out filters trial-wise
light_data = backOutFilters(light_data);

%%
% First, we show all the LN models for all the data we have where we use just light, showing the diversity of gains that the ORN responds with:

plot(light_data,'LN.firing_rate');
prettyFig('fs=13;','FixLogX=true;','lw=1.3;','plw=1.2;')

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, we re-plot all the data as before, but grouping by paradigm, to see if the variability we observe is linked to the different experimental paradigms we subjected these ORns to. 

plot(light_data,'LN.firing_rate','grouping',light_data.paradigm);
prettyFig('fs=13;','FixLogX=true;','lw=1.3;','plw=1.2;')

if being_published
	snapnow
	delete(gcf)
end

%%
% We see that the different curves have different colours, suggesting that the variability in the data reflects different experimental perturbations, which is good. 

%%
% We now find the paradigms with large variance, as this is most likely to have evidence for fast gain control. 

do_this_paradigm = 2;
light_data.use_these_trials = (light_data.paradigm == do_this_paradigm);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax(1) = subplot(1,2,1); hold on
ax(2) = subplot(1,2,2); hold on
plot(light_data,ax,'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3);


%% Version Info
%
pFooter;