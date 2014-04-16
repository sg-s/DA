% FindDAModelParameters
% finds the parameters for the DA model 
% Usage:
% function [] = FindDAModelParameters() run with default values
% or 
% function [] = FindDAModelParameters(filename,x0)
function [x] = FindDAModelParameters(varargin)
if nargin < 1
	filename = '/data/random-stim/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand.mat';
else
	filename = varargin{1};
end
if nargin < 2
	x0 = [100   0  0    0.1234    0.4111    1.0656   0    2.5405]; 
else
	x0 = varargin{2};
end

% load data and prepare
foptions = optimset('Display','iter');

[PID, time, f] = PrepData3(filename);

% make all vectors consistent
PID = PID(:)*100;
time = time(:);
f = f(:)/100;


% run fminsearch on the cost function 
x = fminsearch(@(x) DA_cost_function(x,PID,f),x0,foptions);



