%% adaptive NLN model with bacteria-like feedback of k_D to maintain constant activity at 1/2 
% implements a minimum k_D to prevent it exploding at low stimuli 
% the output nonlinearity is a simple threshold linear function
% uses a Euler solver to solve the differential equation 
% here, we allow the n of the input NL to vary
% 
function [R, a, b, kD, K] = aNLN5(S,p)

% list parameters for clarity 

% input NL parameters
p.k0; % minimum k_D 
p.B; % strength of adaptive changes 
p.Hill_n;

% filter
p.tau_y1;
p.tau_y2;
p.A; % ratio of lobes of bilobed filter 
p.n; % for the filter, not the nonlinearity 

% output scale
p.C;

% bounds

lb.Hill_n = 1;
ub.Hill_n = 8;

lb.k0 = 1e-3;
ub.k0 = .1;

lb.B = 0;
ub.B = 1e3;

lb.tau_y1 = 20;
ub.tau_y1 = 200;

lb.tau_y2 = 50;
ub.tau_y2 = 110;

lb.A = .5;
ub.A = 1;

lb.C = 1;
ub.C = 600;

lb.n = 1;
ub.n = 4;




% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);

R = 0*S;
a = 0*S;
kD = 0*S;

% interpolate 
vS = interp1(time,S,T); vS(isnan(vS)) = S(1);
vkD = 0*vS;
va = 0*vS;

hn = ceil(p.Hill_n);

% use a fixed-step Euler to solve this
for i = 2:length(vS)
	dydt = p.B*vkD(i-1)*(va(i-1) - 1/2);

	vkD(i) = dydt*1e-4 + vkD(i-1);
	if vkD(i) < p.k0
		vkD(i) = p.k0;
	end

	% update a
	KS = (vkD(i)/vS(i))^hn;
	va(i) = 1/(1 + KS);

end

% switch back to timestep of data
a = interp1(T,va,time);
kD = interp1(T,vkD,time);

% generate the filter
q.n = p.n;
q.tau1 = p.tau_y1;
q.tau2 = p.tau_y2;
q.A = p.A;
K = filter_gamma2(1:1e3,q);

b = filter(K,1,a);

b(b<0) = 0;
R = b*p.C;






