%% model based on bacterial chemotaxis for the LFP
% in this model, we consider the full model and keep alpha fixed, 
% and allow the barrier height to vary
% 

function [R, a, w_plus, w_minus, e0] = bacteriaModelX_fixedA_simple(S,p)

% work with matrices too
if size(S,2) > 1
	R = NaN*S; a = NaN*S;  w_plus = NaN*S; w_minus = NaN*S; e0 = NaN*S; 
	for i = 1:size(S,2)
		[R(:,i),a(:,i),w_plus(:,i), w_minus(:,i), e0(:,i)] = bacteriaModelX_fixedA_simple(S(:,i),p);
	end
	return
end

assert(isvector(S),'stimulus must be a vector')
assert(length(p)==1,'2nd argument should be a structure of length 1')


% list parameters for clarity 

% input NL parameters
p.B; % rate of adaptation 
p.e_L; % minimum k_D
p.A;

p.K_1;
p.K_2;

% filter
p.K_tau;
p.n; % for the filter, not the nonlinearity 



% bounds

lb.B = 2.9;
ub.B = 3;

lb.e_L = 0.86;


lb.K_1 = .01;
ub.K_1 = 1;

lb.K_2 = 400;
ub.K_2 = 400;

lb.K_tau = 2;
ub.K_tau = 6;

lb.n = 1;
ub.n = 5;

lb.A = 1;
ub.A = 120;



% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);

% interpolate 
S_ = interp1(time,S,T); S_(isnan(S_)) = S(1);
e0_ = 0*S_;
a_ = 0*S_;

w_minus = 0*S_;
w_plus = 0*S_;

% inital conditions 
a_(1) = .5;
Shat = (1 + S_(1)/p.K_2)/(1 + S_(1)/p.K_1);
e0_(1) =  - log(Shat);
if e0_(1) < p.e_L
	e0_(1) = p.e_L;
end

denom = 1 + S_(1)/p.K_2; 
num_num = S_(1)*(1 + exp(-e0_(1)));
num_denom = (p.K_1)*(1 + (p.K_2/p.K_1)*exp(-e0_(1)));
num = 1 + (num_num/num_denom);
w_plus(1) = p.A*(num/denom)/(1+exp(e0_(1)));

denom = 1 + S_(1)/p.K_1; 
w_minus(1) = p.A*(num/denom)/(1+exp(-e0_(1)));

% use a fixed-step Euler to solve this
for i = 2:length(S_)
	% update e0
	dydt = p.B*(a_(i-1) - 1/2);

	e0_(i) = dydt*1e-4 + e0_(i-1);
	if e0_(i) < p.e_L
		e0_(i) = p.e_L;
	end

	% compute rates 
	denom = 1 + S_(i-1)/p.K_2; 
	num_num = S_(i-1)*(1 + exp(-e0_(i-1)));
	num_denom = (p.K_1)*(1 + (p.K_2/p.K_1)*exp(-e0_(i-1)));
	num = 1 + (num_num/num_denom);
	w_plus(i) = p.A*(num/denom)/(1+exp(e0_(i-1)));

	denom = 1 + S_(i-1)/p.K_1; 
	w_minus(i) = p.A*(num/denom)/(1+exp(-e0_(i-1)));

	% update a
	dydt = w_plus(i)*(1-a_(i-1)) - w_minus(i)*a_(i-1);
	a_(i) = dydt*1e-4 + a_(i-1);
	if a_(i) < 0
		a_(i) = 0;
	end
	if a_(i) > 1
		a_(i) = 1;
	end

end

% switch back to timestep of data
a = interp1(T,a_,time);
e0 = interp1(T,e0_,time);
w_plus = interp1(T,w_plus,time);
w_minus = interp1(T,w_minus,time);



% % pass through a filter 
t = 0:(length(S)/10);
K = generate_simple_filter(p.K_tau,p.n,t);
K = K(1:find(K>1e-2*max(K),1,'last'));
R = filter(K,1,a);


R = R*p.output_scale;
R = R + p.output_offset;

	function f = generate_simple_filter(tau,n,t)
		f = t.^n.*exp(-t/tau); % functional form in paper
		f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
	end

end






