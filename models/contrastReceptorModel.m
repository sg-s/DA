% contrastReceptorModel.m
% simple model of receptor binding and unbinding +
% diffusible factor that depends as a function on the mean stimulus over some tau
% and also on the variance of the stimulus over some tau'
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,b,d] = contrastReceptorModel(S,p)

% receptor parameters
p.r_b; % binding rate
p.theta_b; % ratio constants

% diffusible factor
p.tau_mean;
p.m_mean;

% output nonlinearity 
p.A;
p.R0;



% bounds
lb.r_b = 0;
lb.theta_b = 0;
ub.theta_b = 45;


lb.A = 0;


lb.tau_mean = 1;
lb.m_mean = 0;


b = 0*S;

tau_mean = ceil(p.tau_mean);

d = 0*S;
for i = tau_mean+1:length(S)
	d(i) = std(S(i-tau_mean:i))*mean(S(i-tau_mean:i));
end


for i = 2:length(S)
	fx = (1-b(i-1))*p.r_b*S(i) - d(i-1)*b(i-1)*p.r_b*p.theta_b;
	b(i) = fx + b(i-1);

	if b(i) < 0
		b(i) = 0;
	end

	if b(i) > 1
		b(i) = 1;
	end


end


% R = p.hill_A*(b.^2)./(b.^2 + p.hill_K^2);
R = p.A*b + p.R0;



