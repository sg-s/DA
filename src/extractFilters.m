% extractFilters
% extracts filters and nonlinearities and gains from consolidated data
% usage:
% [K, prediction, gain, gain_err] = extractFilters(X,Y,...)
%
% options (and defualts are:)
%
% band_pass_y 	: 	false
% band_pass_y	: 	false
% filter_length :  	700
% filter_offset	: 	100
% use_cache		: 	true
% a 			: 	1
% z 			: 	(length of data)
%
% created by Srinivas Gorur-Shandilya at 10:30 , 14 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [K, prediction, gain, gain_err] = extractFilters(X,Y,varargin)

% defaults
band_pass_y 	= false;
band_pass_x 	= false;
filter_length 	= 700;
filter_offset 	= 100;
filter_buffer 	= 300;
dt 				= 1e-3;
use_cache = true;
a = 1;
z = length(X);

if ~nargin
	help extractFilters
	return
else
    if iseven(nargin)
    	for ii = 1:2:length(varargin)-1
        	temp = varargin{ii};
        	if ischar(temp)
        		% disp(strcat(temp,'=varargin{ii+1};'))
            	eval(strcat(temp,'=varargin{ii+1};'));
        	end
    	end
	else
    	error('Inputs need to be name value pairs')
	end
end

% orient data correctly
if size(X,1) < size(X,2) 
	X = X';
end
if size(Y,1) < size(Y,2)
	Y = Y';
end

% defensive programming
assert(length(X) == length(Y),'X and Y should be matrices of the same size');
assert(width(X) == width(Y),'X and Y should be matrices of the same size');
assert(iseven(filter_buffer),'filter_buffer should be an even number')

% back out the filters
K = [];
if use_cache
	K = cache(dataHash([X Y]));
end
if isempty(K)
	K = NaN(filter_length,width(X));
	for i = 1:width(X)
		textbar(i,width(X))
		resp = Y(a:z,i);
		rm_this = isnan(resp);
		resp(rm_this) = [];

		if length(resp) 
			try
				if band_pass_y
					resp = bandPass(resp,1e3,10);
				end
				
				stim = X(a:z,i);
				stim(rm_this) = [];
				if band_pass_x
					stim = bandPass(stim,1e3,10);
				end

				temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',filter_length+filter_buffer,'offset',filter_offset+filter_buffer/2);

				K(:,i) = temp(1+(filter_buffer/2):end-(filter_buffer/2));
			catch err
				
			end
		end
	end
	cache(dataHash([X Y]),[]);
	cache(dataHash([X Y]),K);
end

% make the linear prediction and compute the gain
time = dt*(1:length(X));
filtertime = dt*(1:filter_length) - dt*filter_offset;
prediction = NaN*Y;
gain = NaN(width(Y),1);
gain_err = gain;

for i = 1:width(X)
	prediction(:,i) = convolve(time,X(:,i),K(:,i),filtertime);
	% fit lines to estimate gains
	x = prediction(a:z,i);
	if band_pass_y
		y = bandPass(Y(a:z,i),1e3,10);
	else
		y = Y(a:z,i);
	end
	try
		[ff,gof] = fit(x(:),y(:),'poly1');
		gain(i) = ff.p1;
		% get the weights for the each
		temp = confint(ff);
		gain_err(i) = diff(temp(:,1))/2;
	catch
	end
end

