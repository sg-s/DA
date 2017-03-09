%% ORNmodel2
% a unified model for LFP and firing rate responses
% in this model, we consider the LFP-> firing rate transformation to be governmed by a ODE system that does something like a differentiating filter, whose timescale depends on kD

function [R,a,kD,tau,F,H] = ORNmodel2(S,p)


	% receptor binding, slowdown at receptors and adaptation
	p.w; % omega, attempt rate at jumping from bound to unbound states
	p.tau_adap; % receptor adaptation time 
	p.kD_min;

	% speed up in firing rate
	p.tau_firing; 
	p.C;

	% bounds
	lb.w = 0;
	lb.tau_adap = 1;
	lb.K_a = 0;
	lb.B = 0;
	lb.n_firing = 1;
	lb.K_f = 0;
	lb.tau_firing = 0;
	lb.kD_min = 1e-3;


	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S; kD = a; tau = a; F = a; H = a;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),kD(:,i),tau(:,i),F(:,i),H(:,i)] = ORNmodel2(S(:,i),p);
		end
		return
	end


	time = 1e-3*(1:length(S));
	T = 1e-4:1e-4:max(time);

	% interpolate 
	% variable_  stores the interpolated variable, that will be solved in the Euler integrator
	S_ = interp1(time,S,T); S_(isnan(S_)) = S(1);
	kD_ = 0*S_;
	a_ = 0*S_;
	F_ = 0*S_;
	H_ = 0*S_;


	% use a fixed-step Euler to solve this
	for i = 2:length(S_)


		% update k_D
		dydt = (1/p.tau_adap)*kD_(i-1)*(a_(i-1) - 1/2);

		kD_(i) = dydt*1e-4 + kD_(i-1);
		if kD_(i) < p.kD_min
			kD_(i) = p.kD_min;
		end

		% calculate the energies
		F_(i) = log(kD_(i));
		H_(i) = F_(i);

		k_minus = p.w*(exp(-H_(i)));
		k_plus = p.w*(exp(-H_(i) - F_(i)));

		% update a
		dydt = (1-a_(i-1))*k_plus*S_(i-1) - a_(i-1)*k_minus;
		a_(i) = dydt*1e-4 + a_(i-1);

		if a_(i) > 1
			a_(i) = 1;
		elseif a_(i) < 0
			a_(i) = 0;
		end
	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a_,time);
	F = interp1(T,F_,time);
	H = interp1(T,H_,time);
	kD = interp1(T,kD_,time);


	% now we need to filter it, with the timescale of the filter changing with time 
	tau = ceil(p.tau_firing*exp(-kD/p.C));

	% filter by solving a ODE that takes a derivative 
	R = 0*a;
	for i = 2:length(a)
		if i - tau(i) > 1
			dydt = a(i) - a(i-tau(i));
			dydt = dydt/tau(i);
		else
			dydt = 0;
		end
		R(i) = R(i-1) + dydt;

	end


end
