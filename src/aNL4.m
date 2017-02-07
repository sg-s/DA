% this is the LFP variant of aNLN4.m
% differences: 
% monolobed filter instead of a bilobed one
% adds a offset parameter because we don't know the actual absolute value of the LFP
% 
function [R, a, b, kD, K] = aNL4(S,p)

% list parameters for clarity 

% input NL parameters
p.k0; % minimum k_D 
p.B; % strength of adaptive changes 

% this uses pLFPfilter to generate the filter. 
p.n;
p.A;
p.tauA;
p.tauB;


% output scale
p.C;
p.D;

% bounds
lb.k0 = 1e-3;
ub.k0 = 1;

lb.B = 1;
ub.B = 100;

% filter
lb.tauA = 3;
lb.tauB = 10;
ub.tauA = 40;
ub.tauB = 140;

lb.A = 0.5;
ub.A = 1;

lb.C = 1e-3;
ub.C = 6;

lb.n = 1;
ub.n = 4;

lb.D = -100;
ub.D = 10;



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

% use a fixed-step Euler to solve this
for i = 2:length(vS)
	dydt = p.B*vkD(i-1)*(va(i-1) - 1/2);

	vkD(i) = dydt*1e-4 + vkD(i-1);
	if vkD(i) < p.k0
		vkD(i) = p.k0;
	end

	% update a
	va(i) = 1/(1 + (vkD(i)/vS(i)));

end

% switch back to timestep of data
a = interp1(T,va,time);
kD = interp1(T,vkD,time);

% make the filter
q.A = p.A;
q.tauA = p.tauA;
q.tauB = p.tauB;
q.n = p.n;
K = pLFPfilter(1:1e3,q);
b = filter(K,1,a);

% trivial scaling parameter
R = p.D +  b*p.C;






