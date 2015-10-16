% NagelWilsonModel.m
% specification of the Nagel-Wilson model for receptor binding, channel opening, and adptation via diffusable factor
% liberally modified from the original paper as there is no consistent, complete definition of a model in their paper
% the model has 3 parts:
% 
% 1. a two-state receptor-ligand binding system (4 ODEs)
% 2. a diffusable factor that affects channel opening probabilites (2 ODEs)
% 3. a simple point neuron model that converts channel opening probabilities to a voltage 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function dy = NagelWilsonModel(t,y,time,odor,p)

dy = NaN*y;

% calculate the odor at the time point
O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

% receptor-ligand binding

M = [-(p.ka*p.sa+O*p.kb*p.sb)	p.sa 						p.sb  						0;
	p.ka*p.sa 					-(p.sa+O*p.theta*p.kb*p.sb)	0 							p.sb;
	O*p.kb*p.sb 				0     						-(p.sb*p.theta*p.ka*p.sa) 	p.sa;
	0 							O*p.theta*p.kb*p.sb			p.theta*p.ka*p.sa 			-(p.sa+p.sb)];

dy(1:4) = M*y(1:4);


% diffusible factor and channel opening
dy(5) = (y(2)+y(4))*p.ko*(1/y(6))*(1-y(5)) - p.kc*y(5); % according to Nagel and Wilson
dy(6) = p.a*y(5) - p.b*y(6);

% channel opening -> voltage
dy(7) = (p.V0 + p.R*y(5)*(p.Ec-y(7)) - y(7))/p.tau;
