function [Ky,Kz] = DA_Filters(p,t)
%% function [Ky,Kz] = DA_Filters(p)
% plots the Ky and Kz filters from DA model parameters
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

if ~nargin
	help DA_Filters
	return
end

Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);



function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




