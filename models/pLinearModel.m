% Parametric Linear Model
% the filter is modelled by a double gamma filter (two lobed)
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,K,shat] = pLinearModel(s,p)

% set lower bounds
lb.n = 1; 
lb.tau1 = 1;
lb.tau2 = 2;
lb.A = 0;

% set upper bounds
ub.n = 4;
ub.tau1 = 400;
ub.tau2 = 600;
ub.A = 1;



% show parameters for readability
p.n;
p.A;
p.tau1;
p.tau2;
p.scale;
p.offset;

% make the filters
filter_length = 10*max([p.n*p.tau2  p.n*p.tau1]);
if filter_length < length(s)/10
else
	filter_length = length(s)/10; % ridiculously long filters
end
t = 0:filter_length; 
K = filter_gamma2(t,p);


% filter the input
shat = filter(K,1,s);

% add the offset
shat = shat + p.offset;

% scale
f = shat*p.scale;