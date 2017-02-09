% pLFP.m
% 
% created by Srinivas Gorur-Shandilya at 5:40 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% Parametric Linear LFP filter
% the filter is modelled by a single gamma filter (two lobed)
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,K,shat] = pLFP(s,p)

% set lower bounds
lb.n = 2; 
lb.tau = 10;
lb.x_offset = -100;

% set upper bounds
ub.n = 4;
ub.tau = 400;
ub.x_offset = 10;


% show parameters for readability
p.n;
p.tau;
p.A;
p.offset;
p.x_offset;

% make the filters
filter_length = 10*p.tau*p.n;
if filter_length < length(s)/10
else
	filter_length = length(s)/10; % ridiculously long filters
end
t = 0:filter_length; 
K = filter_gamma(t,p);

% offset the stimulus
s = circshift(s,round(p.x_offset));

% filter the input
shat = -filter(K,sum(K),s);

% add the offset
f = shat + p.offset;
