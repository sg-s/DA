% plotNonlinearity.m
% subfunction called by overloaded plot function for class ORNData
%  
% created by Srinivas Gorur-Shandilya at 10:24 , 28 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles] = plotNonlinearity(plot_here,pred,resp,plot_type,nbins)
	switch plot_type
	case 'trial-wise'
		for i = 1:width(pred)
			if any(~isnan(pred(:,i)))
				[~,data] = plotPieceWiseLinear(pred(:,i),resp(:,i),'nbins',nbins,'make_plot',false);
				plot_handles(i) = plot(plot_here,data.x,data.y);
			end
		end
	case 'sem'
		for i = 1:width(pred)
			[~,data(i)] = plotPieceWiseLinear(pred(:,i),resp(:,i),'nbins',nbins,'make_plot',false);
		end
		plot_handles = errorShade(plot_here,nanmean([data.x],2),nanmean([data.y],2),nanmean([data.ye],2));
	end	
	
end