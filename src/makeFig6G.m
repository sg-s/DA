% makeFig6G.m
% computes instantenous gain vs. preceding stimulus intensity, as in fig 6G in this paper:
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003289
% 
% created by Srinivas Gorur-Shandilya at 10:12 , 23 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [x,y] = makeFig6G(stimulus,response,prediction,history_length)
if ~nargin
	help makeFig6G
	return
end
if ~isvector(stimulus) || ~isvector(response) || ~isvector(prediction)
	error('First two argumnets should be vectors')
else
	stimulus = stimulus(:);
	response =  response(:);
	prediction = prediction(:);
	if length(stimulus) ~= length(response)
		error('Stimulus and response have to have the same length')
	end
	if length(prediction) ~= length(response)
		error('prediction and response have to have the same length')
	end
end
if length(history_length) > 1
	error('History Length has to be a scalar')
end

% check cache for results
temp.stimulus = stimulus;
temp.response = response;
temp.prediction = prediction;
temp.history_length = history_length;
hash = dataHash(temp);
cached_data = cache(hash);
if ~isempty(cached_data)
	x = cached_data.x;
	y = cached_data.y;
	return
end

x = NaN*stimulus;
y = NaN*stimulus;

parfor i = (history_length+1):length(stimulus)
	x(i) = mean(stimulus(i-history_length:i));
	if any(isnan(prediction(i-history_length:i))) || any(isnan(response(i-history_length:i)))
	else
		cf = fit(prediction(i-history_length:i),response(i-history_length:i),'Poly1');
		y(i) = cf.p1;
	end

end

clear temp
temp.x = x;
temp.y = y;
cache(hash,temp);

