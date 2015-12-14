% findBestGainFilter.m
% this is meant to be optimised by fitModel2Data, and is a function that accepts a stimulus vector, and parameterically filters it and returns a number that describes how good it is at explaining a vector of instantenous gain (which is stored as a global)
%
% created by Srinivas Gorur-Shandilya at 11:04 , 13 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function r2 = findBestGainFilter(S,p)

r2 = 0;

global inst_gain	

% parameters for clarity
p.A;
p.tau1;
p.tau2;
p.n;

% bounds
ub.tau1 = 300;
ub.tau2 = 400;
ub.A = 1;
ub.n = 10;

lb.tau1 = 1;
lb.tau2 = 1;
lb.A = 0;
lb.n = .1;

% make the filter
K = filter_gamma2(1:2e3,p);

if any(isnan(K))
	return
end

% filter the stimulus
shat = filter(K,1,S);

% remove some junk
rm_this = isnan(inst_gain) | isnan(shat);
temp = inst_gain(~rm_this);
shat(rm_this) = [];

try
	% use the Spearman rank correlation
	r2 = abs(spear(temp(1:10:end),shat(1:10:end)));

catch
end
