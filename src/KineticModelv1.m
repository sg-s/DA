% non-adapting kinetic model, variant 1
% here, all the response kinetics come from receptor binding and unbinding. 
% 
function [R, a] = KineticModelv1(S,p)

% list parameters for clarity 
p.n; % co-operativity 
p.k_plus;
p.k_minus;

% output scale
p.A;
p.B;

% bounds
lb.n = 0.1;
lb.k_plus = 0;
lb.k_minus = 0;


ub.n = 10;
ub.k_plus = 1e3;
ub.k_minus = 1e3;


% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);


% interpolate 
vS = interp1(time,S,T); vS(isnan(vS)) = S(1);
va = 0*vS;

n = ceil(p.n);

% use a fixed-step Euler to solve this
for i = 2:length(vS)
	dydt = p.k_plus*(vS(i-1)^n)*(1-va(i-1)) - p.k_minus*va(i-1);

	va(i) = dydt*1e-4 + va(i-1);
	if va(i) < 0
		va(i) = 0;
	end

	if va(i) > 1
		va(i) = 1;
	end
end

% switch back to timestep of data
a = interp1(T,va,time);

% trivial scaling parameter
R = p.B +  a*p.A;






