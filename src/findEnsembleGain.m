% findEnsembleGain.m
% calculates gain in small bins over the ensemble of data
% usage:
% G = findEnsembleGain(X,Y,varargin)
% where X and Y are M x N matrixes (M time points, N trials), where N >> 1
% and G is a vector as long as M
% 
% created by Srinivas Gorur-Shandilya at 9:29 , 10 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [gain,gain_r2,stimulus_contrast] = findEnsembleGain(X,Y,S,varargin)

% options and defaults
options.window_size = 50;
options.R_window = [.33 .66]; % window over which to calculate the gain 
options.step_size = 10;

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options = setfield(options,temp,varargin{ii+1});
    	end
    end
end
else
	error('Inputs need to be name value pairs')
end

hash1 = dataHash(options);
hash2 = dataHash([X;Y;S]);
hash = dataHash([hash1 hash2]);
data = cache(hash);
if ~isempty(data)
	gain = data.gain;
	gain_r2 = data.gain_r2;
	stimulus_contrast = data.stimulus_contrast;
	return
end


gain = NaN(size(X,1),1); gain_r2 = gain;
stimulus_contrast = gain;
a = round(options.window_size/2) + 1;
z = length(gain) - round(options.window_size/2) - 1;
ws = round(options.window_size/2)-1;

options.step_size = round(options.step_size);
assert(options.step_size > 0,'Step size must be positive')
assert(options.step_size < z/10,'Step size too large')

for i = a:options.step_size:z
	textbar(i,z)
	% grab the snippet of data around this point across the ensemble
	x = X(i-ws:i+ws,:); x = x(:);
	y = Y(i-ws:i+ws,:); y = y(:);
	s = S(i-ws:i+ws,:); s = s(:);
	[~,data] = plotPieceWiseLinear(x,y,'nbins',100,'make_plot',false);
	x = data.x; y = data.y; clear data
	% find the correct portion to calculate the gain on
	lb = find((y-min(y))/(max(y)-min(y)) > options.R_window(1),1,'first');
	ub = find((y-min(y))/(max(y)-min(y)) < options.R_window(2),1,'last');
	if any(~isnan(y))
		ff = fit(x(lb:ub),y(lb:ub),'poly1');
		gain(i) = ff.p1;
		gain_r2(i) = rsquare(x(lb:ub),y(lb:ub));

		% also calculate the stimulus contrast
		stimulus_contrast(i) = nanstd(s)/nanmean(s);
	end
end

data.stimulus_contrast = stimulus_contrast;
data.gain = gain;
data.gain_r2 = gain_r2;
cache(hash,data);




