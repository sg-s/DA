%% modified version of DA model designed to slow down
% Euler implementation of the ODE solver, because everything else is just too slow

function [R,y,z] = LFPmodel2_Euler(S,p)

% bounds
lb.s0 = -25;
ub.s0 = 20;

lb.A = 1;
lb.B = 10;
lb.tau_A = 1000;
lb.tau_r = 1;
lb.tau_y = 10;
lb.tau_z = 20;

ub.tau_z = 1e3;
ub.tau_y = 300;
ub.tau_r = 5;
ub.B = 30;



% list parameters for legibility
p.A;
p.B;
p.tau_A;
p.tau_y;
p.tau_z;
p.tau_r;
p.s0;

t = 0:3000; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,2,t);
Kz = generate_simple_filter(p.tau_z,2,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

% in the case that tau_r is 0, don't have to integrate; this approximation can
% also be used when tau_r/(1+p.B*z) << other time scales in y or z, for all or most z
if p.tau_r == 0
    R = p.A*y./(1+p.B*z);
    R = R(:);
    return;
end

R = 0*S;

for i = 2:length(S)
	tau_eff = ((1+p.tau_A*z(i))/(1+p.B*z(i)))*p.tau_r;
	dR = p.A*y(i) - (1 + p.B*z(i))*R(i-1);
	dR = (1/tau_eff)*dR;
	R(i) = R(i-1) + dR;
end

R = R + p.s0;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




