% IFF.m
% incoherent feed-forward model from
% 10.7554/eLife.06694
% 
% created by Srinivas Gorur-Shandilya at 2:23 , 22 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [y,u] = IFF(S,p)


    %    b1: 9.9004
    %    b2: 6.5817
    %    b3: 0.0052
    %    b4: 0.2777
    %    b5: 0.0092
    %    a1: 5.4297
    %    a2: 0.0028
    % theta: 0.4964

% define all parameters
p.b1;
p.b2;
p.b3;
p.b4;
p.b5;

p.a1;
p.a2;

p.theta;

% lower bounds for all parameters is 0
lb.b1 = 1;
lb.b2 = 0;
lb.b3 = 0;
lb.b4 = 0;
lb.b5 = 0;

lb.a1 = 0;
lb.a2 = 0;

lb.theta = 0;


% some sensible upper bounds
ub.theta = 30;
ub.b2 = 10;
ub.b3 = 10;
ub.b4 = 10;




% define initial condition
y = zeros(length(S),1);
u = zeros(length(S),1);
x = zeros(length(S),1);

y(1) = 0;
u(1) = 0;

% strictly assuming time step is 1ms, and that all parameters are in ms units

for i = 2:length(S)
	fx = (p.b1*(S(i-1)))/(p.b2 + S(i-1) + u(i-1)*p.b3) - (p.b4*y(i-1)^2)/(y(i-1)^2 + p.theta^2) - p.b5*y(i-1);
	y(i) = y(i-1) + fx;
	if y(i) < 0
		y(i) = 0;
	end

	fx = p.a1*S(i-1) - p.a2*u(i-1);
	u(i) = u(i-1) + fx;
	if u(i) < 0
		u(i) = 0;
	end

end	