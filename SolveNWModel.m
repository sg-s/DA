% SolveNWModel.m
% numerically integrates the Nagel-Wilson Model
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,T2] = SolveNWModel(time,odor,p,plothere)

if nargin < 1
	% specify time
	time = 0:1e-3:4;
end
if nargin < 2
	% specify stimulus
	odor = 0*time;
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
	p.tau = .001;
	p.V0 = -70;
	p.R = 1;
	p.Ec= 10;
	p.stim_scale = 1;
	p.f0 = 10; % mean output firing rate
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
		load('NWtestdata.mat');
		time = t;
	end
	if isempty(odor)
		load('NWtestdata.mat');
		odor = o;
		time = [-10:mean(diff(time)):-mean(diff(time)) time];
		odor = [zeros(1,length(time)-length(odor)) odor];
		throw_away = -1;
	else
		% transient
		time = [-10:mean(diff(time)):-mean(diff(time)) time];
		odor = [zeros(1,length(time)-length(odor)) odor];
		throw_away = -1;
	end
end

% scale stimulus
odor = odor.*(10^p.stim_scale);
odor(isnan(odor)) = 0;


% initial condition
ic= zeros(7,1);
ic(1) = .1; ic(2) = .1; ic(3) = .1; ic(4) = .1; 
ic(5) = 0;
ic(6) = .01; 
ic(7)=-50;



Tspan = [min(time) max(time)];
[T, Y] = ode23s(@(t,y) NagelWilsonModel(t,y,time,odor,p),Tspan,ic); % Solve ODE

% throw away the transient
Y(T<throw_away,:) = [];
T(T<throw_away) = [];


% pass the voltage through the "spike filter"
t = 1:50; % in ms
K= normpdf(t,10,3).*(-8e-3*t+.0812);
K = K/max(K);
K = K*20;
% interpolate to get a uniform sampling, and to match the sampling of the "spike filter"
T2 = min(T):1e-3:max(T);
Y2 = NaN(length(T2),width(Y));
for i = 1:width(Y)
	Y2(:,i) = interp1(T,Y(:,i),T2);
end
f = filter(K,1,Y2(:,7)-mean(Y2(:,7))) + p.f0;


if ~isempty(plothere)
	% plot the output
	l = {'R','R*','OR','OR*','C_o','D','V_m'};
	for i = 1:7
		plot(plothere(i),T,Y(:,i),'k.')
		legend(plothere(i),l{i})
	end
	f(1:1000) = NaN;
	plot(plothere(8),Treal,R,T2,f)
	set(plothere(8),'YLim',[-5 200],'XLim',[min(T2)-5 max(T2)+5])
end

disp(p)