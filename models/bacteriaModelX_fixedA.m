%% model based on bacterial chemotaxis for the LFP
% in this model, we consider the full model and keep alpha fixed, 
% and allow the barrier height to vary
% 

function [R, a, w_plus, w_minus, e0] = bacteriaModelX_fixedA(S,p)

% work with matrices too
if size(S,2) > 1
	R = NaN*S; a = NaN*S;  w_plus = NaN*S; w_minus = NaN*S; e0 = NaN*S; 
	for i = 1:size(S,2)
		[R(:,i),a(:,i),w_plus(:,i), w_minus(:,i), e0(:,i)] = bacteriaModelX_fixedA(S(:,i),p);
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

p.Delta;
p.Gamma; 

p.K_1;
p.K_2;

% filter
p.K_tau;
p.n; % for the filter, not the nonlinearity 



% bounds

lb.B = 2;
ub.B = 1e3;

lb.e_L = 0;

lb.K_1 = .01;
ub.K_1 = .01;

lb.K_2 = 10;
ub.K_2 = 1e5;

lb.K_tau = 1;
ub.K_tau = 50;

lb.n = 1;
ub.n = 3;

% reduce the model a bit
lb.Delta = 0;
ub.Delta = 0;

lb.Gamma = 0;
ub.Gamma = 0;



% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);

% interpolate 
S_ = interp1(time,S,T); S_(isnan(S_)) = S(1);
e0_ = 0*S_;
a_ = 0*S_;

w_minus = 0*S_;
w_plus = 0*S_;

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
	num_denom = (p.K_1)*(1 + (p.K_2/p.K_1)*exp(-e0_(i-1))*exp(-p.Gamma));
	num = 1 + (num_num/num_denom)*exp(-p.Delta);
	w_plus(i) = p.A*(num/denom)*exp(-p.Gamma)/(1+exp(e0_(i-1)));

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






