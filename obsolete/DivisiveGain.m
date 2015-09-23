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

% list parameters for readability
p.tau1;
p.tau2;
p.n;
p.A;
p.beta;


% specify bounds
lb.tau1 = 1;  
lb.tau2 = 1;  
ub.n = 10;
lb.n = 2;
lb.beta = 0;

if isvector(s)
	error('Divisive Gain requires that the first argument be a 2xN matrix, where the first row is the stimulus and the second row is the prediction to be gain-corrected')
end
if size(s,1) > size(s,2)
	s = s';
end


fp = s(2,:);
s = s(1,:);

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 5*max([p.n*p.tau1 p.n*p.tau2]);
if filter_length < length(s)/10
else
	filter_length = length(s)/10; % ridiculously long filters
end
if filter_length < 300
	filter_length = 300;
end
t = 0:filter_length; 

Ka = filter_gamma2(t,p);

% run through adaptation filter
gain = filter(Ka,1,s-mean(s))+1;

% scale
gain = gain*p.beta;

f = fp./(1+gain);



