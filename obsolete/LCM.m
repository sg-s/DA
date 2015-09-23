function [R,y,z] = LCM(S,p)
%% function [R,y,z] = LCM(S,p)
% function takes argument of stimulus and a parameter set and implements the Linearised Cascade Model (LCM)
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
if ~nargin
	help LCM
	return
end

S = (S - p.s0); 

t = [0:1000]; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

R = (p.A*y./(1+p.B*z));

% % pass through rectifier
% R(R<0) = 0;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




