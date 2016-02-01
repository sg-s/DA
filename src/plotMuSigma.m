% plotMuSigma.m
% 
% created by Srinivas Gorur-Shandilya at 5:43 , 01 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles,r2] = plotMuSigma(plot_here,X,history_length,make_plot)

history_length = ceil(history_length);

% pad with NaNs to divide by history_length
X = [X; NaN(history_length - rem(length(X),history_length),1)];
X = reshape(X,history_length,length(X)/history_length);
if make_plot
	plot_handles = plot(plot_here,nanmean(X),nanstd(X),'Marker','.','Color',[.5 .5 .5],'LineStyle','none');

	% equalise axes and plot a guide
	temp = [get(plot_here,'XLim') get(plot_here,'YLim')];
	temp = temp(:); 
	set(plot_here,'XLim',[min(temp) max(temp)],'YLim',[min(temp) max(temp)])
	plot(plot_here,[1e-6 max(temp)],[1e-6 max(temp)],'k--')
else
	plot_handles = [];
	mx = nanmean(X);
	sx = nanstd(X);

	if length(mx) > 1e3
		r2 = rsquare(mx(1:100:end),sx(1:100:end));
	else
		r2 = rsquare(mx,sx);

	end
		
end