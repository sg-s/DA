% SolveNWModel.m
% numerically integrates the Nagel-Wilson Model
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = SolveNWModel(time,odor,p,plothere)

if nargin < 1
	% specify time
	time = 0:.001:4;
end
if nargin < 2
	% specify stimulus
	odor = 0*time;
	% odor(1000:2000) = 0;
end

if nargin < 3
	% parameters
	p.ka = .001; 
	p.sa = 10; % 1/s
	p.kb = 1;
	p.sb = 100; % 1/s
	p.theta = 100; % scaling factor
	p.ko = 500;
	p.kc = 10;
	p.a = .8; % alpha, in the diffusible factor eq. (1/s)
	p.b = .6; % beta, in the diffusible factor eq. (1/s)
	p.tau = 1;
	p.V0 = -70;
	p.R = 1;
	p.Ec= 10;
	p.stim_scale = 1e-3;
end
if nargin < 4
	figure, hold on;
	for i = 1:8
		plothere(i)=subplot(4,2,i);
	end
end
if nargin == 4
	% check arguments OK
	if isempty(time)
		time = 0:.001:20;
	end
	if isempty(odor)
		odor = 0*time;
		odor(10000:11000) = 1;
	end
end

% scale stimulus
odor = p.stim_scale*odor;

% initial condition
ic= zeros(7,1);
ic(1) = .1; ic(2) = .1; ic(3) = .1; ic(4) = .1; 
ic(5) = 0;
ic(6) = .01; 
ic(7)=-50;

Tspan = [0 max(time)];
[T, Y] = ode23s(@(t,y) NagelWilsonModel(t,y,time,odor,p),Tspan,ic); % Solve ODE

% plot the output
l = {'R','R*','OR','OR*','C_o','D','V_m'};
for i = 1:7
	plot(plothere(i+1),T,Y(:,i),'k.')
	legend(plothere(i+1),l{i})
end
plot(plothere(1),time,odor)
