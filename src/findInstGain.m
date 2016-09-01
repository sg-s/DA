% findInstGain.m
% computes instantenous gain vs. preceding stimulus intensity, as in fig 6G in this paper:
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003289
% 
% created by Srinivas Gorur-Shandilya at 10:12 , 23 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [gain,gain_err] = findInstGain(response,prediction,window_size,skip)
if ~nargin
	help findInstGain
	return
end
if ~isvector(response) || ~isvector(prediction)
	error('First two arguments should be vectors')
else
	response =  response(:);
	prediction = prediction(:);
	if length(prediction) ~= length(response)
		error('prediction and response have to have the same length')
	end
end
if length(window_size) > 1
	error('History Length has to be a scalar')
end

% check cache for results
temp.response = response;
temp.prediction = prediction;
temp.window_size = window_size;
hash = dataHash(temp);
cached_data = cache(hash);
if ~isempty(cached_data)
	gain = cached_data.gain;
	gain_err = cached_data.gain_err;
	return
end

disp('cache miss:')
disp(hash)

gain = NaN*response;
gain_err = NaN*response;

for i = (window_size+1):skip:length(response)
	if any(isnan(prediction(i-window_size:i))) || any(isnan(response(i-window_size:i)))
	else
		textbar(i,length(response))
		[cf,gof] = fit(prediction(i-window_size:i),response(i-window_size:i),'Poly1');
		gain(i) = cf.p1;
		gain_err(i) = gof.rsquare;
	end

end

clear temp
temp.gain = gain;
temp.gain_err = gain_err;
cache(hash,[])
cache(hash,temp);

