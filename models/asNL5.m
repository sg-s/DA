%% asNL5.m
% NL model that slows down as it adapts
% this is similar to asNL3, but has been completely re-written. Also, we solve for the evolving k_D, then calcualte the free energy F that then determines k+,k-. this is more numerically robust.

function [R,a,kD,F,H] = asNL5(S,p)

	% list parameters for legibility

	% receptors + adaptation
	p.k0; % minimum possible k_D
	p.w; % attempt frequency
	p.B; % amount of adaptation, also the rate of adaptation

	% filter + scale
	p.K_tau;
	p.output_offset;
	p.output_scale;
	p.n;

	% bounds 
	lb.n = 1;
	lb.K_tau = 10;
	lb.k0 = 1e-3;
	lb.w = 0;
	lb.B = 0;

	ub.K_tau = 100;
	ub.n = 1;
	ub.w = 2e3;
	ub.k0 = .1;
	ub.B = 100;




	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S;  kD = NaN*S; F = NaN*S; H = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),kD(:,i), F(:,i), H(:,i)] = asNL5(S(:,i),p);
		end
		return
	end

	assert(isvector(S),'stimulus must be a vector')
	assert(length(p)==1,'2nd argument should be a structure of length 1')

	time = 1e-3*(1:length(S));
	T = 1e-4:1e-4:max(time);

	% interpolate 
	% variable_  stores the interpolated variable, that will be solved in the Euler integrator
	S_ = interp1(time,S,T); S_(isnan(S_)) = S(1);
	kD_ = p.k0 + 0*S_;
	a_ = 0*S_;
	F_ = 0*S_;
	H_ = 0*S_;

	F_max = log(p.k0);

	% use a fixed-step Euler to solve this
	for i = 2:length(S_)
		% update k_D
		dydt = p.B*kD_(i-1)*(a_(i-1) - 1/2);

		kD_(i) = dydt*1e-4 + kD_(i-1);
		if kD_(i) < p.k0
			kD_(i) = p.k0;
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

	% % pass through a filter 
	t = 0:(length(S)/10);
	K = generate_simple_filter(p.K_tau,p.n,t);
	K = K(1:find(K>1e-2*max(K),1,'last'));
	R = filter(K,1,a);

	R = R + p.output_offset;
	R = R*p.output_scale;

		function f = generate_simple_filter(tau,n,t)
			f = t.^n.*exp(-t/tau); % functional form in paper
			f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
		end

end
