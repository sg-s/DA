% Parametric LN Model
% the filter is modelled by a double gamma filter (two lobed)
%
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,K,shat] = pLNModel(s,p)

% explicitly state the parameters

% for the filter
p.tau1;
p.tau2;
p.A;
p.n;

% offset, 
p.offset;

% hill nonlinearity
p.Hill_A;
p.Hill_Kd;
p.Hill_n;

% bounds
lb.n = 2;
ub.n = 4;
lb.A = 1;
lb.tau1 = 1;
lb.tau2 = 2;
lb.Hill_n = 1;
ub.Hill_n = 10;
ub.Hill_A = 1e3;
lb.Hill_A = 50;
ub.tau2 = 500;
ub.tau1 = 200;
lb.Hill_Kd = 1;

% make the filters
filter_length = 4*max([p.n*p.tau2  p.n*p.tau1]);
if filter_length < length(s)/10
else
	filter_length = length(s)/10; % ridiculously long filters
end
t = 0:filter_length; 
K = filter_gamma2(t,p);

% filter the input
shat = filter(K,1,s-mean(s));

% add the offset
shat = shat + p.offset;

% pass through the non-linearity
x = [p.Hill_A p.Hill_Kd p.Hill_n];
f = hill(x,shat);
