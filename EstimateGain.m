% EstimateGain.m
% estimates the gain of a system by comparing response to stimulus
% by taking small snippets, fitting lines to them, and reporting the slopes
% it needs two parameters: the window over which to compute the slope (tw)
% and the sliding it uses (s)
% use it as:
% [gain,r2] = EstimateGain(stimulus,response,tw,s);
% 
% created by Srinivas Gorur-Shandilya at 10:05 , 13 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [gain,r2,t] = EstimateGain(stimulus,response,tw,s)
switch nargin
case 0
	help EstimateGain
	return
case 1
	error('Not enough input arguments')
case 2
	tw = 100;
	s = 10;
case 3
	s = tw/10;
end

if tw < 10 || s < 1
	error('Either the window or the sliding are too small.')
else
	s = floor(s);
	tw = 2*floor(tw/2);
end

if length(stimulus) ~= length(response)
	error('Stimulus and response have to be the same length')
end

% cache support 
temp.stimulus = stimulus;
temp.response = response;
temp.tw = tw;
temp.s = s;
temp.method = 'PCA';
hash = DataHash(temp);

cached_data = cache(hash);
if ~isempty(cached_data)
	temp = cached_data;
	gain = temp.gain;
	r2 = temp.r2;
	t = temp.t;
	return
end

% do the grunt work
t = tw/2:s:(length(stimulus)-tw/2);
gain = NaN(length(t),1);
r2 = gain;


parfor i = 1:length(gain)
	x = stimulus((t(i)-tw/2+1):(t(i)+tw/2-1));
	y = response((t(i)-tw/2+1):(t(i)+tw/2-1));
	if any(isnan(x)) | any(isnan(y))
	else

		% use PCA to get slopes of clouds of points
		[coeff,~,latent] = pca([x y]);
		gain(i) = coeff(2,1)/coeff(1,1);
		low_gof(i) = latent(1)/sum(latent);


		% [f,g] = fit(x,y,'poly1');
		% gain(i) = f.p1;
		% r2(i) = g.rsquare;
	end
end



% cache results 
clear temp
temp.gain = gain;
temp.r2 = r2;
temp.t = t;
cache(hash,temp);
