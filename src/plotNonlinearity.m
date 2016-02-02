% plotNonlinearity.m
% subfunction called by overloaded plot function for class ORNData
%  
% created by Srinivas Gorur-Shandilya at 10:24 , 28 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

	
function plot_handles = plotNonlinearity(plot_here,pred,resp,grouping,plot_options)
if isempty(grouping)
	% no grouping
	grouping = ones(width(pred),1);
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
	this_pred = pred(:,grouping == this_group);
	this_resp = resp(:,grouping == this_group);

	if strcmp(plot_options.ioCurve_type,'pwlinear')
		if strcmp(plot_options.plot_type,'trial-wise')
			for j = 1:width(this_pred)
				 plot_handles(i).line(j) = plotPieceWiseLinearN(plot_here,this_pred(:,j),this_resp(:,j),[c(i,:) .3],plot_options.nbins);
			end
			% also super-impose the mean
			plot_handles(i).line(j+1) = plotPieceWiseLinearN(plot_here,nanmean(this_pred,2),nanmean(this_resp,2),c(i,:),plot_options.nbins);
			set(plot_handles(i).line(j+1),'LineWidth',3)
		elseif strcmp(plot_options.plot_type,'sem')
			axes(plot_here);
			plot_handles(i) = plotPieceWiseLinear(this_pred,this_resp,'nbins',plot_options.nbins,'make_plot',true,'Color',c(i,:),'normalise',plot_options.normalise);
		elseif strcmp(plot_options.plot_type,'mean')
			plot_handles(i) = plotPieceWiseLinearN(plot_here,nanmean(this_pred,2),nanmean(this_resp,2),c(i,:),plot_options.nbins);
		end
	else strcmp(plot_options.ioCurve_type,'dots')
		if strcmp(plot_options.plot_type,'trial-wise')
			for j = 1:width(this_pred)
				 plot_handles(i).line(j) = plot(plot_here,this_pred(1:plot_options.nbins:end,j),this_resp(1:plot_options.nbins:end,j),'Color',c(i,:),'Marker','.','LineStyle','none');
			end
		elseif strcmp(plot_options.plot_type,'sem')
			error('Cannot plot SEM with dot plot type')
		elseif strcmp(plot_options.plot_type,'mean')
			plot_handles(i) = plot(plot_here,nanmean(this_pred(1:plot_options.nbins:end,:),2),nanmean(this_resp(1:plot_options.nbins:end,:),2),'Color',c(i,:),'Marker','.','LineStyle','none');
		end
	end

end

	function plot_handle = plotPieceWiseLinearN(plot_here,x,y,colour,nbins)
		[~,data] = plotPieceWiseLinear(x,y,'nbins',nbins,'make_plot',false);
		plot_handle = plot(plot_here,data.x,data.y,'Color',colour);

	end




end

