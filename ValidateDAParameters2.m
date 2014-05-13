% valdiates parameters for DA model (v2)
function [p] = ValidateDAParameters2(x)
p.A = x(1); % alpha in model
p.B = x(2); % beta in the model
p.C = (x(3)); % gamma in the model
p.tau_y = (x(4));
p.n_y = x(5);
p.tau_z = (x(6));
p.n_z = x(7);
p.r0 = x(8);

