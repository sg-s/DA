% 
% 
% created by Srinivas Gorur-Shandilya at 10:57 , 05 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Fast Gain Control Analysis of Natural Stimulus
% In this document we perform a fast gain control analysis on our naturalistic stimulus data, using the new ORNData class and associated methods, and using the same analysis we use in our analysis of fast gain control. This analysis is unbiased re. the structure of the signal, and doesn't care about whiffs, etc. All it does is fit lines in 50ms windows, and report on the slope and try to explain the variation in slopes using the stimulus averaged over some recent past.

pHeader;

%%
% In the following figure, we plot the inst. gain vs. the mean stimulus for three different hsitory lengths. We observe that the correlation between the two seems to be maximally negative away from either extreme. Finally, plotting the Spearman correlation vs the histroy length, we see there is a sharp peak at ~200ms, and static gain control (the $\rho$ at $\tau_{H} =0$) can only explain less than 50% of what we see. 

load('/Users/sigbhu/code/da/data/nat_stim_parametric_fits.mat','od')

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
hl = [10 500 1e4];
for i = 1:3
	plot_handles = plot(od,[ax(i) ax(4)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i),'data_bin_type','dots','nbins',3);
	title(ax(i),['\tau_{H} = ' oval(hl(i)) 'ms'])
	if i < 3
		delete(plot_handles(2).f2)
	end
end

set(ax(1:3),'XScale','log','YScale','log')
set(ax(4),'YLim',[-1 0])

prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

pFooter;

