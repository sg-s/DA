% [K] = FitFilter2Data(stim, response,filter_length,reg,n,OnlyThesePoints)
% fits a linear filter to data.
% based on Chichlinsky's STA method
% created by Srinivas Gorur-Shandilya at 16:04 , 15 January 2014. 
% Contact me at http://srinivas.gs/contact/
% 
% Usage: 
% 
% K = FitFilter2Data(stim, response);
% this is the minimal usage. you must specify stim and response, and they must be vectors of the same length.
% 
% Optional arguments 
% 
% filter_length:    how many points do you want K to have? scalar, default is 333
% reg:              regularisation factor. scalar, default is 50. makes filter ring less. 
% n:                0 or 1. Do you want to subtract mean before computing filter? default is 0. 
% OnlyThesePoints:  vector that contains matrix coordinates of points in time that you want to 
%                   restrict filter computation on. max(OnlyThesePoints) must be < 
%                   length(stim)-filter_length
% regtype:			1: addition of r I to the covariance matrix. 2: addition of rI/(1+r). 3: refer
% 					to code.
%
% Defaults:
% filter_length = 333;
% reg = 50;
% n = 1;
% OnlyThesePoints = 1:length(stim)-filter_length; (everything)
% regtype = 3;
% 
% Example: 
% 
% K=FitFilter2Data(stim, response,'filter_length=500;','n=1;','OnlyThesePoints=1:2000;')
% calculates a 500-point filter from the data after removing mean and regularising, but only at the first 2000 points of the data.
function [K] = FitFilter2Data(stim, response,varargin)

% defaults
filter_length = 333;
reg = 50;
n = 1;
OnlyThesePoints = 1:length(stim)-filter_length;
regtype = 3;

% evaluate optional inputs
for i = 1:nargin-2
	eval(varargin{i})
end


% check that stimulus and response are OK
if ~isvector(stim)
	error('Stimulus is not a vector')
end
if ~isvector(response)
	error('Response is not a vector')
end

% ensure column
stim = stim(:);
response = response(:);

% subtract mean in both
if n
	response = response - mean(response);
	stim = stim - mean(stim);
end


% throw away first bit of response, because we don't have the stimulus before it
response = response(filter_length+1:end);

% chop up the stimulus into blocks  
s = zeros(length(OnlyThesePoints), filter_length+1);
for i=OnlyThesePoints
    s(i,:) = stim(filter_length+i:-1:i);
end

switch regtype 
	case 1
		C = s'*s + reg*eye(filter_length+1); % Carlotta's reg.
	case 2
		C = s'*s + reg*eye(filter_length+1)/(1+reg);  % my modification
	case 3
		C = s'*s; % Damon's reg
		C = (C + reg*eye(filter_length+1))*trace(C)/(trace(C) + reg*filter_length);

end
   
K = C\(s'*response);
        



