% Parametric filter
% the filter is modelled by a double gamma filter (two lobed)
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function Khat = pFilter(K,p)

% set bounds
lb.tau1 = 1;
lb.tau2 = 1;
lb.n = 1;
lb.A = 0;
ub.tau1 = 1e2;
ub.tau2 = 1e3;


% show parameters for readability
p.n;
p.A;
p.tau1;
p.tau2;

% make the filters
if isempty(K)
	Khat = filter_gamma2(1:300,p);
else	
	Khat = filter_gamma2(1:length(K),p);
end

