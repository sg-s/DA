%% adaptive model with adaptation occuring in channel opening, not receptors

function [R, b, d] = adaptingChannelModel(S,p)

% list parameters for clarity 

p.kD; % single, fixed kD, does not change
p.k1; % on rate for channel opening
p.k2; % off rate for channel closing 

p.n_channels; % number of channels 
p.C; % rate at which the diffusible factor is cleared out

% trivial parameters
p.A;
p.B;

% bounds
lb.kD = 1e-3;
lb.k1 = 1e-3;
lb.k2 = 1e-3;
lb.n_channels = 50;
lb.C = 1e-3;

ub.k1 = 10;
ub.k2 = 10;
ub.kD = 10;

ub.C = 1e3;
ub.n_channels = 1e3;

% compute the receptor bound fraction
a = S./(S + p.kD);

% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);

R = 0*S;
b = 0*S;
d = 0*S;


% interpolate 
va = interp1(time,a,T); va(isnan(va)) = a(1);
vb = 0*va;
vd = 0*va;


% use a fixed-step Euler to solve this
for i = 2:length(va)

	k1 = p.k1*exp(-vd(i-1));
	k2 = p.k2*exp(-vd(i-1));
	dydt = k1*va(i-1)*(1-vb(i-1)) - k2*vb(i-1);
	dydt = dydt/(1+ vd(i-1));

	vb(i) = dydt*1e-4 + vb(i-1);
	if vb(i) < 0
		vb(i) = 0;
	end
	if vb(i) > 1
		vb(i) = 1;
	end

	% update d
	dydt = p.n_channels*vb(i-1) - p.C*vd(i-1);

	vd(i) = dydt*1e-4 + vd(i-1);

	if vd(i) < 0
		vd(i) = 0;
	end


end

% switch back to timestep of data
b = interp1(T,vb,time);
d = interp1(T,vd,time);

R = p.A*b + p.B;





