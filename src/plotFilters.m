% plotFilters.m
% 
% created by Srinivas Gorur-Shandilya at 2:18 , 02 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plot_handles = plotFilters(plot_here,filtertime,K,grouping,plot_options)
if isempty(grouping)
	% no grouping
	grouping = ones(width(K),1);
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

	% get the data for this group
	this_K = K(:,grouping == this_group);

	if strcmp(plot_options.plot_type,'trial-wise')
		for j = 1:width(this_K)
			 plot_handles(i).line(j) = plot(plot_here,filtertime,this_K(:,j),'Color',[c(i,:) .3]);
		end
		% also super-impose the mean
		plot_handles(i).line(j+1) = plot(plot_here,filtertime,nanmean(this_K,2),'Color',c(i,:),'LineWidth',3);
	elseif strcmp(plot_options.plot_type,'sem')
		axes(plot_here);
		plot_handles(i).h = errorShade(filtertime,nanmean(this_K,2),sem(this_K'),'Color',c(i,:));
	elseif strcmp(plot_options.plot_type,'mean')
		plot_handles(i) = plot(filtertime,nanmean(this_K,2),'Color',c(i,:));
	end


end


