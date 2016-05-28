%% NagelWilsonIntegrate.m
% This function is a wrapper around NagelWilsonModel
% and solves that ode system using ode23t
% the formulation of this fuction is meant to be called by manipulate and other functions
function [C,D,R] = NagelWilsonIntegrate(S,p)

% list parameters for legibility
p.theta;
p.ka;
p.kb;
p.sa;
p.sb;
p.k0;
p.A;
p.B;


% define initial conditions for ode system
ic = zeros(6,1);
ic(1) = .1; ic(2) = .1; ic(3) = .1; ic(4) = .1; 
ic(5) = .4;
ic(6) = .4; 



time = 1e-3*(1:length(S));
Tspan = [min(time) max(time)];

[T, Y] = ode23t(@(t,y) NagelWilsonModel(t,y,time,S,p),Tspan,ic); % Solve ODE

% re-interpolate the solution to fit the stimulus
C = interp1(T,Y(:,5),time);
D = interp1(T,Y(:,6),time);

R = zeros(length(S),4);
for i = 1:4
	R(:,i) = interp1(T,Y(:,i),time);
end
