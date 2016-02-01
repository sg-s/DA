% plotTimeSeries.m
% 
% created by Srinivas Gorur-Shandilya at 6:05 , 01 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plot_handles = plotTimeSeries(plot_here,time,data,grouping,plot_options)

if isempty(grouping)
	% no grouping
	grouping = ones(width(data),1);
end

groups = unique(grouping);
% make a colour scheme
if length(groups) == 1
	c = [0 0 0];
else
	c = parula(length(groups)+1);
end

for i = 1:length(groups)
	this_group = groups(i);

	this_data = data(:,grouping == this_group);
	if strcmp(plot_options.plot_type,'trial-wise')
		plot_handles(i).h = plot(plot_here,time,this_data,'Color',[c(i,:) .3]);
		plot_handles(i).h(end+1) = plot(plot_here,time,nanmean(this_data,2),'Color',c(i,:),'LineWidth',2);
	elseif strcmp(plot_options.plot_type,'sem')
		axis(plot_here);
		plot_handles = errorShade(time,nanmean(data,2),sem(data'),'Color',c(i,:));
	elseif strcmp(plot_options.plot_type,'mean')
		plot_handles = plot(plot_here,time,nanmean(data,2),'Color',c(i,:));
	end


end