%% model based on bacterial chemotaxis for the LFP
% 
% 

function [R, a, w_plus, w_minus, e0] = bacteriaModelX_coop(S,p)

% work with matrices too
if size(S,2) > 1
	R = NaN*S; a = NaN*S;  w_plus = NaN*S; w_minus = NaN*S; e0 = NaN*S; 
	for i = 1:size(S,2)
		[R(:,i),a(:,i),w_plus(:,i), w_minus(:,i), e0(:,i)] = bacteriaModelX_coop(S(:,i),p);
	end
	return
end

assert(isvector(S),'stimulus must be a vector')
assert(length(p)==1,'2nd argument should be a structure of length 1')


% list parameters for clarity 

% input NL parameters
p.B; % rate of adaptation 
p.e_L; % minimum k_D
p.w0; % attempt frequency

p.m; % cooperativity 

p.K_on;
p.K_off;

% filter
p.K_tau;
p.n; % for the filter, not the nonlinearity 



% bounds
lb.m = 10;
ub.m = 12;

lb.B = 10;
ub.B = 1e3;

lb.e_L = 0;
ub.e_L = 1e3;

lb.w0 = 1e2;
ub.w0 = Inf;

lb.K_tau = 17;
ub.K_tau = 17;

lb.n = 1;
ub.n = 1;

lb.K_on = .1;
ub.K_on = 1;

lb.K_off = 5;
ub.K_off = 10;


lb.output_scale = 0;

% solve the ODE
time = 1e-3*(1:length(S));
T = 1e-4:1e-4:max(time);

% interpolate 
S_ = interp1(time,S,T); S_(isnan(S_)) = S(1);
e0_ = 0*S_;
a_ = 0*S_;

w_minus = 0*S_;
w_plus = 0*S_;

m = ceil(p.m);

a_(1) = 1/2;
Shat = (1 + S_(1)/p.K_off)/(1 + S_(1)/p.K_on);
e0_(1) = - log(Shat);

% use a fixed-step Euler to solve this
for i = 2:length(S_)
	% update e0
	dydt = p.B*(a_(i-1) - 1/2);

	e0_(i) = dydt*1e-4 + e0_(i-1);
	if e0_(i) < p.e_L
		e0_(i) = p.e_L;
	end

	% compute rates 
	w_plus(i) = p.w0*((2*p.K_off)/(p.K_off + S_(i-1)))^m;
	w_plus(i) = w_plus(i)*exp(-e0_(i-1));

	w_minus(i) = p.w0*((2*p.K_on)/(p.K_on + S_(i-1)))^m;


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






