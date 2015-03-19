% EstimateGain2
% estimates the absolute gain for a dataset
% usage:
% 
% g = EstimateGain2(stim,resp)
%
% where 
% stim and resp are matrices of equal lengths
% g is a vector as long as resp is wide. 
%
% EstimateGain works by fitting a nonparametric filter to the stimulus and the response, normalising the filter, and then fitting a line to the residuals. the slope of the line is defined as the gain, and the units are units of response/units of stimulus. 
% 
% created by Srinivas Gorur-Shandilya at 10:25 , 19 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function g = EstimateGain2(stim,resp)

if ~nargin
	help EstimateGain2
	return
else
	if length(stim) ~= length(resp)
		error('stim and length should be of equal length')
	end
	if width(stim) ~= width(resp)
		error('stim and length should be of equal width')
	end
end	

% preallocate output
g = NaN(width(stim),1);

% orient inputs correctly
if size(stim,2) > size(stim,1)
	stim = stim';
end
if size(stim,2) > size(stim,1)
	stim = stim';
end

% check cache
x.stim = stim;
x.resp = resp;
cached_data = cache(DataHash(x));
if ~isempty(cached_data)
	g = cached_data;
	return
end

for i = 1:width(stim)
	% throw out NaNs
	rm_this = isnan(stim(:,i)) | isnan(resp(:,i));
	this_stim = stim(:,i);
	this_resp = resp(:,i);
	this_stim(rm_this) = [];
	this_resp(rm_this) = [];

	% back out the linear filter
	[K, ~, filtertime_full] = FindBestFilter(this_stim,this_resp,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full;
	filtertime = (-200:900);
	K = interp1(filtertime_full,K,filtertime);
	K = K/max(K);

	% make a linear prediction
	fp = convolve(1:length(this_resp),this_stim,K,filtertime);

	% trim NaNs again
	rm_this = isnan(fp);
	fp(rm_this) = [];
	this_resp(rm_this) = [];

	temp=fit(fp(:),this_resp(:),'poly1');
	g(i) = temp.p1;

end


% cache data
cache(DataHash(x),g);

