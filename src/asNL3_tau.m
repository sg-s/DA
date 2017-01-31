% returns steady state timescale of asNL3 model
% doesn't solve the ODE -- so much faster

function [tau] =  asNL3_tau(S,p)

mS = mean(S);

F = p.kT*(log(mS*p.w_plus/p.w_minus));

H = p.E0 - p.E1*F;

k_plus = p.w_plus*exp(-H/p.kT);
k_minus = p.w_plus*exp(-(H+F)/p.kT);

tau = 1./(mS*k_plus + k_minus) + 0*S;