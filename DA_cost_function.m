% created by Srinivas Gorur-Shandilya at 12:27 , 05 December 2013. Contact me at http://srinivas.gs/contact/
% DA_cost_function is a function that evaluates the response predicted by the DA model
% to the stimulus and and compares it to the actual response. it calcualtes the absolute
% error of the prediction, so has a minimum of 0 when the prediction is perfect. 
function [cost]  = DA_cost_function(x)
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
p2.A = 1; % alpha in model
p2.B = 0.2; % beta in the model
p2.C = 0.5; % gamma in the model
p2.tau_r = 30; % these are the timescales
p2.tau_y = 50;
p2.n_y = 2;
p2.tau_z = 60;
p2.n_z = 3;

% now find the actual result
Ractual = DA_integrate(S,p2);

% now find the result from the guess
Rguess = DA_integrate(S,p);

cost = sum(abs(Rguess-Ractual));