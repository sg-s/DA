% EstimatingGainTimescaleByFastFlickers.m
% 
% created by Srinivas Gorur-Shandilya at 3:36 , 27 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Estimating gain timescales using fast flickers
% In this document, we carefully re-analyse some of Carlotta's old data, where she present flickering odor stimuli of the same odour to the same ORN, changing only the correlation time of the stimulus. The nice thing about this dataset is that it was all taken on the same day, removing many experimental confounds and sources of variation. 

pHeader;

% load the data
load('Carlotta_Data.mat');


do_these = [3 4 2];
corr_label = {'\tau = 30ms','\tau = 50ms','\tau = 100ms'};
figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on

for i = 1:length(do_these)
	clear h
	h(1) = subplot(2,3,i); hold on
	h(2) = subplot(2,3,i+3); hold on

	plot(orn_data(do_these(i)),h,'binned_gain_analysis')
	title(h(1),corr_label{i})
end

prettyFig()

if being_published		
	snapnow		
	delete(gcf)	
end

%%
% We now fit a DA model to the 100ms data, and explicitly model a gain control mechanism that depends on a timescale > 100ms. We now repeat the analysis on this dataset, to see if we get a similar result. 

clear p
p.   s0 = 3.8594e-04;
p.  n_z = 2;
p.tau_z = 100;
p.  n_y = 2;
p.tau_y = 10;
p.    C = 0.0438;
p.    A = 4.1276e+04;
p.    B = 648.4775;

sim_data = ORNData;

% generate the data
for i = 1:length(do_these)
	this_stim = orn_data(do_these(i)).stimulus;
	this_resp = 0*this_stim;
	for j = 1:width(this_stim)
		temp = fit((1:length(this_stim))',this_stim(:,j),'poly2');
		this_stim(:,j) = this_stim(:,j) - temp((1:length(this_stim))) + mean(this_stim(:,j));
		this_resp(:,j) = DAModelv2(this_stim(:,j),p);
	end
	sim_data(i).stimulus = this_stim;
	sim_data(i).firing_rate = this_resp;

	sim_data(i) = backOutFilters(sim_data(i));

end

`

prettyFig()



%% Version Info
%
pFooter;