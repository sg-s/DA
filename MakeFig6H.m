% MakeFig6H.m
% makes data for a plot as in fig 6H in this paper:
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003289
% 
% created by Srinivas Gorur-Shandilya at 11:20 , 23 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function c = MakeFig6H(stimulus,response,history_length)
if ~nargin
	help MakeFig6H
	return
end
if ~isvector(stimulus) || ~isvector(response) 
	error('First two arguments should be vectors')
else
	stimulus = stimulus(:);
	response =  response(:);
	if length(stimulus) ~= length(response)
		error('Stimulus and response have to have the same length')
	end
end
if length(history_length) > 1
	error('History Length has to be a scalar')
end

c = NaN*stimulus;

parfor i = (history_length+1):length(stimulus)
	x = stimulus(i-history_length:i);
	y = response(i-history_length:i);
	x = x - mean(x);
	y = y - mean(y);
	y = y/max(y);
	x = x/max(x);
	if any(isnan(x)) || any(isnan(y)) || max(y) == y(end) || max(y) == y(1) || max(x) == x(end) || max(x) == x(1)
	else
		xc = xcov(x,y);
		[~,c(i)] = max(xc(history_length/2:history_length*(1.5)));
		c(i) = c(i) - history_length/2;
		if abs(c(i)) < 5
			
		end
	end

end
