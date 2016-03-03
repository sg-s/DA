% fitPieceWiseNL.m
% method for ORNData class
% 
% created by Srinivas Gorur-Shandilya at 1:37 , 23 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [o] = fitPieceWiseNL(o)

assert(strcmp(class(o),'ORNData'),'Argument should be a ORNData class object')

% recursively call this to work on arrays of objects
original_dimensions = size(o);
o = o(:);
if length(o) > 1
	for i = 1:length(o)
		o(i) = fitPieceWiseNL(o(i));
	end
	o = reshape(o,original_dimensions);
	return
end


% respect constrains on where to calculate this 
if isempty(o.use_these_trials)
	utt = 1:o.n_trials;
else
	utt = o.use_these_trials;		
end
if isempty(o.use_this_segment)
	uts = 1:length(o.stimulus);
else
	uts = o.use_this_segment;
end	
% check if data is there
if ~isempty(o.firing_rate) && ~isempty(o.firing_projected) && ~isempty(o.stimulus)
	r = o.firing_rate(uts,:);
	p = o.firing_projected(uts,:);

	for i = 1:width(r)
		if ismember(i,find(utt))
			temp1 = p(:,i); temp2 = r(:,i);
			rm_this = isnan(temp1) | isnan(temp2);
			temp1(rm_this) = []; temp2(rm_this) = [];
			if length(temp1) > 10
				[~,data] = plotPieceWiseLinear(temp1,temp2,'nbins',50,'make_plot',false);
				ff = fit(data.x,data.y,'cubicinterp');
				o.firing_projected(:,i) = ff(o.firing_projected(:,i));
				o.firing_projected((o.firing_projected(:,i) < 0),i) = 0;
			end
		end
	end

end
if ~isempty(o.LFP) && ~isempty(o.LFP_projected) && ~isempty(o.stimulus)
	r = o.LFP(uts,utt);
	p = o.LFP_projected(uts,utt);

	for i = 1:width(r)
		textbar(i,width(r))
		if ismember(i,find(utt))
			temp1 = p(:); temp2 = r(:);
			rm_this = isnan(temp1) | isnan(temp2);
			temp1(rm_this) = []; temp2(rm_this) = [];
			keyboard
		end
	end
end	