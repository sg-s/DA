% MakeFig6H.m
% makes data for a plot as in fig 6H in this paper:
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003289
% 
% created by Srinivas Gorur-Shandilya at 11:20 , 23 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function c = MakeFig6H_trialwise(stimulus,response,history_length)
if ~nargin
	help MakeFig6H_trialwise
	return
end
if isvector(stimulus) || isvector(response) 
	error('First two arguments should be matrices')
else
	if size(stimulus,2) > size(stimulus,1)
		stimulus = stimulus';
	end
	if size(response,2) > size(response,1)
		response = response';
	end

end
if length(history_length) > 1
	error('History Length has to be a scalar')
end

c = NaN*stimulus;

parfor j = 1:size(stimulus,2)
	disp(j)
	this_c = NaN(length(stimulus),1);
	for i = (history_length+1):length(stimulus)
		x = stimulus(i-history_length:i,j);
		y = response(i-history_length:i,j);
		x = x - mean(x);
		y = y - mean(y);
		y = y/max(y);
		x = x/max(x);
		if any(isnan(x)) || any(isnan(y)) 
			% || max(y) == y(end) || max(y) == y(1) || max(x) == x(end) || max(x) == x(1)
		else
			xc = xcorr(y,x);
			[~,this_c(i)] = max(xc);
			this_c(i) = this_c(i) - history_length;
		end

	end
	c(:,j) = this_c;
end
