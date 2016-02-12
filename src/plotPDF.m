% 
% 
% created by Srinivas Gorur-Shandilya at 12:42 , 01 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles] = plotPDF(plot_here,X,grouping,plot_options)

if isempty(grouping)
	% no grouping
	grouping = ones(width(X),1);
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
	this_X = X(:,grouping == this_group);

	if strcmp(plot_options.plot_type,'trial-wise')
		for j = 1:width(this_X)
			[y,x] = histcounts(this_X(:,j),plot_options.nbins);
			x = mean(diff(x)) + x(1:end-1);
			y = y/sum(y);
			plot_handles(i).line(j) = plot(plot_here,x,y,'Color',[c(i,:) .33]);
		end
		% also do the mean
		[y,x] = histcounts(this_X,plot_options.nbins);
		x = mean(diff(x)) + x(1:end-1);
		y = y/sum(y);
		plot_handles(i).line(j+1) = plot(plot_here,x,y,'Color',c(i,:),'LineWidth',3);
	elseif strcmp(plot_options.plot_type,'sem')
		% first get the bins right
		[y,x] = histcounts(this_X,plot_options.nbins);
		y = NaN(length(y),width(this_X));
		for j = 1:width(this_X)
			[y(:,j),~] = histcounts(this_X(:,j),x);
			y(:,j) = y(:,j)/sum(y(:,j));
			
		end
		x = mean(diff(x)) + x(1:end-1);
		axis(plot_here);
		plot_handles(i).line = errorShade(plot_here,x,nanmean(y,2),sem(y'),'Color',c(i,:));
	elseif strcmp(plot_options.plot_type,'mean')
		% first get the bins right
		[y,x] = histcounts(this_X,plot_options.nbins);
		y = NaN(length(y),width(this_X));
		for j = 1:width(this_X)
			[y(:,j),~] = histcounts(this_X(:,j),x);
			y(:,j) = y(:,j)/sum(y(:,j));
		end
		x = mean(diff(x)) + x(1:end-1);
		plot_handles(i).line = plot(plot_here,x,nanmean(y,2),'Color',c(i,:));
	end


end
