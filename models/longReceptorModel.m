% longReceptorModel.m
% simple model of receptor binding and unbinding +
% diffusible factor that increases the rate of unbinding 
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [b,d] = longReceptorModel(S,p)

% parameters
p.r_b; % binding rate
p.r_d; % diffusible factor rate

% ratio constants
p.theta_b;
p.theta_d;

p.k;

% define initial conditions
ic = [.1; .1];

time = 1e-3*(1:length(S));
Tspan = [min(time) max(time)];

options = odeset('MaxStep',.1);
[T, Y] = ode23t(@(t,y) lRM(t,y,time,S,p),Tspan,ic,options); % Solve ODE

% re-interpolate the solution to fit the stimulus
Y = interp1(T,Y,time);
b = Y(:,1); d = Y(:,2);

	function dy = lRM(t,y,time,odor,p)
		% calculate the odor at the time point
		O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

		b = y(1);
		d = y(2);

		db = (1-b)*p.r_b*O - p.r_b*p.theta_b*b*d;
		dd = p.r_d*((1-b)/b)*exp(b/p.k) - p.r_d*p.theta_d*d;

		dy = [db; dd];
	end

end




