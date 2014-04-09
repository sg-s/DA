% valdiates parameters for DA model. ensures somethings are +ve, some are bounded, etc.
function [p] = ValidateDAParameters(x)
p.A = x(1); % alpha in model
p.B = x(2); % beta in the model
p.C = (x(3)); % gamma in the model
p.tau_y = (x(4));
p.n_y = x(5);
p.tau_z = (x(6));
p.n_z = x(7);

% set tau to 0
p.tau_r = 0;