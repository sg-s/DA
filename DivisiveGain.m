% DivisiveGain.m
% this model divides a previous estimate of the response by a divisive gain term (like in the DA model)
% and generates a new prediction
% 
% created by Srinivas Gorur-Shandilya at 1:59 , 02 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [f,gain,Ka] = DivisiveGain(s,p)
switch nargin
case 0
	help DGModelv1
	return
case 1
	error('No parameters specified')
end

% specify bounds
lb.tau = 0; lb.n = 0; 

ub.n = 10;

if isvector(s)
	error('Divisive Gain requires that the first argument be a 2xN matrix, where the first row is the stimulus and the second row is the prediction to be gain-corrected')
end
if size(s,1) > size(s,2)
	s = s';
end


fp = s(2,:);
s = s(1,:);

% make the filters
t = 1:300;
Ka = filter_gamma(p.tau,p.n,1,t);

% run through response filter
gain = filter(Ka,1,s-mean2(s));

% scale
gain = gain*p.beta;

f = fp./(1+gain);


