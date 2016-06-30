% simpleReceptorModelX.m
% simple model of receptor binding and unbinding +
% diffusible factor that decreases the rate of binding 
% diffusible factor depends on bound fraction only, not on stimulus
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%% SchulzeLouisModel.m
% model specified in this paper:
% https://www.ncbi.nlm.nih.gov/pubmed/26077825
% http://dx.doi.org/10.7554/eLife.06694
% https://elifesciences.org/content/4/e06694



function [R,D] = simpleReceptorModelX(S,p)

	% list parameters for legibility
	p.r_b;
	p.theta_b;
	p.r_d;
	p.theta_d;

	ic = [.1 .1];

	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	options = odeset('MaxStep',.1);
	[T, Y] = ode23t(@(t,y) simpleReceptorModelX_ode(t,y,time,S,p),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	R = interp1(T,Y(:,1),time);
	D = interp1(T,Y(:,2),time);

		function dY = simpleReceptorModelX_ode(t,Y,time,odor,p)
			% calculate the odor at the time point
			O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

			b = Y(1);
			d = Y(2);

			dY = 0*Y;
			dY(1) = ((1-b)*p.r_b*O)/(d) - b*p.r_b*p.theta_b;
			dY(2) = b*p.r_d - d*p.r_d*p.theta_d;
		end
end