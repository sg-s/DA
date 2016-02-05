% FastGainControlLFP.m
% 
% created by Srinivas Gorur-Shandilya at 1:53 , 02 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Fast Gain Control in Receptors?
% In this document, we try to see if we can see evidence of fast gain control at the receptor level. This dataset has 5 ORNs, with 5-10 trials each where we record both the LFP and the firing rate. 

pHeader;

% od = ORNData;
% for i = 1:5
% 	od(i).stimulus = PID(20e3:end,orn==i);
% 	od(i).firing_rate = fA(20e3:end,orn==i);
% 	od(i).LFP = filtered_LFP(20e3:end,orn==i);
% end

% load the data
if ~exist('od','var')
	load('/Users/sigbhu/code/da/data/LFP_Inst_gain.mat','od') 

	for i = 1:length(od)
		od(i).use_this_segment = 0*od(i).stimulus(:,1);
		od(i).use_this_segment(1e3:end) = true;
		% find the inst. gain
		disp(i)
		od(i) = computeInstGain(od(i));
	end

end



%%
% In the following figures, we look at four different things in the data. First, we plot the firing rate vs. the projected stimulus, and look at the $r^2$ to ensure that our linear models are doing a good job (they are, equally well for the LFP and the firing rate). Second, we plot the distribution of instantaneous gain to see how much the gain varies over this experiment. Having established that the gain does vary a bit (or not), we then plot the inst. gain vs. the mean stimulus in some preceding history length, for three different history lengths (10ms, 500ms, 10s). Finally, we plot the Spearman correlation vs. the history length to see if there is any window in which the mean stimulus predicts the inst. gain well. 


for i = 1:5
	figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
	for j = 1:8
		ax(j) = subplot(2,4,j); hold on
	end

	% show the prediction vs. the data for the LFP
	plot(od(i),ax(1),'ioCurve.LFP','data_bin_type','dots','plot_type','mean','showr2',true);


	% show the dist. of the inst. gain for the LFP
	plot(od(i),ax(2),'pdf.inst_gain_LFP','nbins',100);

	% 
	plot(od(i),ax(3:4),'instGainAnalysis.LFP.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',[10 500 1e4],'data_bin_type','dots','nbins',10);

	set(ax(3),'XScale','log','YScale','log')

	% now do the same for the firing rate
		% show the prediction vs. the data for the firing rate
	plot(od(i),ax(5),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);


	% show the dist. of the inst. gain for the firing_rate
	plot(od(i),ax(6),'pdf.inst_gain_firing','nbins',100);

	% 
	plot(od(i),ax(7:8),'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',[10 500 1e4],'data_bin_type','dots','nbins',10);

	set(ax(7),'XScale','log','YScale','log')

	suptitle(['ORN #' oval(i)])

	prettyFig('fs=13;','FixLogX=true;','lw=1.3;','plw=1.2;');

	if being_published
		snapnow
		delete(gcf)
	end
end


%% Version Info
%
pFooter;