% created by Srinivas Gorur-Shandilya at 12:27 , 05 December 2013. Contact me at http://srinivas.gs/contact/
% DA_cost_function is a function that evaluates the response predicted by the DA model
% to the stimulus and and compares it to the actual response. it calcualtes the absolute
% error of the prediction, so has a minimum of 0 when the prediction is perfect. 
function [cost]  = DA_cost_function_dummy(x)
% convert the inputs into the parameter array that DA_integrate needs
p.A = x(1); % alpha in model
p.B = x(2); % beta in the model
p.C = x(3); % gamma in the model
p.tau_r = x(4); % these are the timescales
p.tau_y = x(5);
p.n_y = x(6);
p.tau_z = x(7);
p.n_z = x(8);

% create the stimulus
ints = 2.^[2:8];
S = reshape([zeros(length(ints),200),ints'*ones(1,1000),zeros(length(ints),500)]',1,[]);

% these are the parameters we want to find out (in reality, these will not be known)
x2(1) = 1; % alpha in model
x2(2) = 0.2; % beta in the model
x2(3)= 0.5; % gamma in the model
x2(4)= 30; % these are the timescales
x2(5) = 50;
x2(6) = 2;
x2(7) = 60;
x2(8) = 3;

Rguess = x.^2;
Ractual = x2.^2;

cost = sum(abs(Rguess-Ractual));