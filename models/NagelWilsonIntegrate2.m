%% NagelWilsonIntegrate2.m
% This function is a wrapper around NagelWilsonModelReduced
% and solves that ode system using ode23t
% the formulation of this fuction is meant to be called by manipulate and other functions
function [C,D,R] = NagelWilsonIntegrate2(S,p)

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
ic = zeros(5,1);
ic(1) = .1; ic(2) = .1; ic(3) = .1; 
ic(4) = 0.4;
ic(5) = .4; 


time = 1e-3*(1:length(S));
Tspan = [min(time) max(time)];

options = odeset('MaxStep',.1);
[T, Y] = ode23t(@(t,y) NagelWilsonModelReduced(t,y,time,S,p),Tspan,ic,options); % Solve ODE

% re-interpolate the solution to fit the stimulus
C = interp1(T,Y(:,4),time);
D = interp1(T,Y(:,5),time);

temp = zeros(length(S),3);
for i = 1:3
	temp(:,i) = interp1(T,Y(:,i),time);
end

% recover the OR species and insert it
R = zeros(length(S),4);
R(:,1:2) = temp(:,1:2);
R(:,4) = temp(:,3);
R(:,3) = 1- sum(temp,2);

