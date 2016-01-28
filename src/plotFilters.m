% plotFilters
% subfunction called by the overloaded plot method in class ORNData
%
% created by Srinivas Gorur-Shandilya at 10:24 , 28 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles] = plotFilters(plot_here,filtertime,K,plot_type)
	switch plot_type
	case 'trial-wise'
		plot_handles = plot(plot_here,filtertime,K);
	case 'sem'
		plot_handles = errorShade(plot_here,filtertime,nanmean(K,2),sem(K'));
	end	
	xlabel(plot_here,'Filtertime (s)')
	ylabel(plot_here,'Filter')
end