%% model based on bacterial chemotaxis for the firing rate
% in this model, receptors are considered to be at steady state at all times 
% so the only ODE we are solving is the one for the k_D (or equivalent)
% similar to original model, but with no minimum k_D, and with a parameter that governs the fixed point of the activity. 
% 
function [R, a, b, e0] = bacteriaModelF2(S,p)

% work with matrices too
if size(S,2) > 1
	R = NaN*S; a = NaN*S;  b = NaN*S; e0 = NaN*S; 
	for i = 1:size(S,2)
		[R(:,i),a(:,i),b(:,i), e0(:,i)] = bacteriaModelF(S(:,i),p);
	end
	return
end

assert(isvector(S),'stimulus must be a vector')
assert(length(p)==1,'2nd argument should be a structure of length 1')


% list parameters for clarity 

% input NL parameters
p.B; % rate of adaptation 
p.a0; % fixed point of a

p.e_L; % lower bound on e0

p.K_1;
p.K_2;

% filter
p.tau_y1;
p.tau_y2;
p.A; % ratio of lobes of bilobed filter 
p.n; % for the filter, not the nonlinearity 

% output scale
p.C;

% bounds

lb.B = .1;
ub.B = 40;

lb.e_L = .1;

lb.a0 = 0.1;
ub.a0 = 0.5;

lb.K_1 = .01;
ub.K_1 = .01;


lb.K_2 = 10;
ub.K_2 = 1e3;

lb.tau_y1 = 17;
ub.tau_y1 = 127;

lb.tau_y2 = 56;
ub.tau_y2 = 296;

lb.A = .6;
ub.A = .8;

lb.C = 1;
ub.C = 1e4;

lb.n = 1;
ub.n = 3;




% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);

% interpolate 
S_ = interp1(time,S,T); S_(isnan(S_)) = S(1);
e0_ = 0*S_;
a_ = 0*S_;

% inital conditions 
a_(1) = p.a0;
Shat = (1 + S_(1)/p.K_2)/(1 + S_(1)/p.K_1);
e0_(1) = log((1-p.a0)/p.a0) - log(Shat);
if e0_(1) < p.e_L
	e0_(1) = p.e_L;
end

% use a fixed-step Euler to solve this
for i = 2:length(S_)
	dydt = p.B*(a_(i-1) - p.a0);

	e0_(i) = dydt*1e-4 + e0_(i-1);

	if e0_(i) < p.e_L
		e0_(i) = p.e_L;
	end

	% update a
	Shat = (1 + S_(i-1)/p.K_2)/(1 + S_(i-1)/p.K_1);
	E = exp(e0_(i-1) + log(Shat));
	a_(i) = 1/(1 + E);

end

% switch back to timestep of data
a = interp1(T,a_,time);
e0 = interp1(T,e0_,time);

% generate the filter
q.n = p.n;
q.tau1 = p.tau_y1;
q.tau2 = p.tau_y2;
q.A = p.A;
K = filter_gamma2(1:1e3,q);

b = filter(K,1,a);

b(b<0) = 0;
R = b*p.C;




