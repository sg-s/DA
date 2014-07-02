% valdiates parameters for XJW model (v2)
function p = ValidateXJWParameters(x)
if length(x) ~= 11
	error('Cannot validate parameters for the XJW model correctly')
end

% calcium dynamics
p.Cm = x(1); % membrance capacitance
p.C = x(2); % amount of calcium entering the neuron per spike
p.tau_Ca = x(3);

% conductances
p.gL = x(4);
p.gAHP = x(5);

% voltages
p.Vreset = x(6);
p.Vth = x(7);
p.Vrest = x(8);
p.Vk = x(9);

% stimulus
p.A = x(10);
p.B = x(11);

