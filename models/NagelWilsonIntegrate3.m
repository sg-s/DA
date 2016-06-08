%% NagelWilsonIntegrate3.m
% This function is a wrapper around NagelWilsonModelReduced
% and solves that ode system using ode23t
% the formulation of this fuction is meant to be called by manipulate and other functions
function [C,D,R,Rx,ORx,OR] = NagelWilsonIntegrate3(S,p)

S = S(:);

% list parameters for legibility
p.theta; 
p.   ka; 
p.   kb; 
p.   sa;  
p.    A; 
p.    B; 
p.   ko; 
p.   kc; 


% define initial conditions for ode system
ic = .1+zeros(3,1);

ic(1) = 1/(1 +p.ka + S(1)*p.kb*(1+p.theta*p.ka));


time = 1e-3*(1:length(S));
Tspan = [min(time) max(time)];

options = odeset('MaxStep',.1);
[T, Y] = ode23t(@(t,y) NagelWilsonModelReduced3(t,y,time,S,p),Tspan,ic,options); % Solve ODE

% re-interpolate the solution to fit the stimulus
R = interp1(T,Y(:,1),time);
C = interp1(T,Y(:,2),time);
D = interp1(T,Y(:,3),time);

R = R(:);

% recover the other species
Rx = (1 - R - S.*R*p.kb)./(1 + p.theta*p.kb*S);
ORx = p.theta*p.kb*S.*Rx;
OR = 1 - R(:) - ORx(:) - Rx(:);

