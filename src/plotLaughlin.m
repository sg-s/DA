% plotLaughlin.m
% 
% created by Srinivas Gorur-Shandilya at 9:59 , 02 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plot_handles = plotLaughlin(plot_here,X,Y,grouping,plot_options)

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

	this_X = X(:,grouping == this_group);
	this_Y = Y(:,grouping == this_group);

	if strcmp(plot_options.plot_type,'trial-wise')
		% plot the CDF
		for j = 1:width(this_X)
			[y,x] = histcounts(this_X(:,j),plot_options.nbins);
			x = mean(diff(x)) + x(1:end-1);
			y = cumsum(y/sum(y));
			plot_handles(i).line(j) = plot(plot_here,x,y,'Color',[c(i,:) .33],'LineStyle','--');
		end
		% also do the mean
		[y,x] = histcounts(this_X,plot_options.nbins);
		x = mean(diff(x)) + x(1:end-1);
		y = cumsum(y/sum(y));
		plot_handles(i).line(j+1) = plot(plot_here,x,y,'Color',c(i,:),'LineWidth',3,'LineStyle','--');

		% now plot the normalised neuron I/o curve
		plot_options.normalise = true;
		plot_options.ioCurve_type = 'pwlinear';
		plot_handles(i).io_curve = plotNonlinearity(plot_here,X,Y,grouping,plot_options);

	elseif strcmp(plot_options.plot_type,'sem')
		% first get the bins right
		[y,x] = histcounts(this_X,plot_options.nbins);
		y = NaN(length(y),width(this_X));
		for j = 1:width(this_X)
			[y(:,j),~] = histcounts(this_X(:,j),x);
			y(:,j) = cumsum(y(:,j)/sum(y(:,j)));
			
		end
		x = mean(diff(x)) + x(1:end-1);
		axis(plot_here);
		plot_handles(i).h = errorShade(plot_here,x,nanmean(y,2),sem(y'),'Color',c(i,:));
		set(plot_handles(i).h,'LineStyle','--')

		% now plot the normalised neuron I/o curve
		plot_options.normalise = true;
		plot_options.ioCurve_type = 'pwlinear';
		plot_options.plot_type = 'sem';
		plot_handles(i).io_curve = plotNonlinearity(plot_here,X,Y,grouping,plot_options);
	end
end