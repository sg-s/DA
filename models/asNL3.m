%% asNL3.m
% NL model that slows down as it adapts
% binding and unbinding rates are updated via an internal parameter, and do not directly "know" the stimulus, which is more realistic 

function [R,a,k_D,F,H] = asNL3(S,p)

	% list parameters for legibility
	p.E0;
	p.E1;
	p.w_minus;
	p.w_plus;

	p.tau_adap;
	p.kT; 
	

	% output filter
	p.K_tau;

	% bounds 
	lb.E0 = 0;
	lb.E1 = 0;
	lb.w_minus = 0;
	lb.w_plus = 0;
	lb.tau_adap = 100;
	lb.kT = 0;
	lb.K_tau = 0;


	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S;  k_D = NaN*S; F = NaN*S; H = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),k_D(:,i), F(:,i), H(:,i)] = asNL3(S(:,i),p);
		end
		return
	end


	time = 1e-3*(1:length(S));
	T = 1e-4:1e-4:max(time);

	% interpolate 
	vS = interp1(time,S,T); vS(isnan(vS)) = S(1);

	% make vectors for variables
	a = 0.5 + 0*vS;
	F0 = p.kT*log(p.w_minus/(p.w_plus*vS(1)));
	F = F0 + 0*vS;
	H = p.E0 + 0*F;
	k_minus = 0*F;
	k_plus = 0*F;

	for i = 2:length(vS)

		% compute barrier height
		H(i) = p.E0 - p.E1*F(i-1);

		% compute rate constants 
		k_plus(i) = p.w_plus*exp(-H(i-1)/p.kT);
		k_minus(i) = p.w_minus*exp(-(F(i-1) + H(i-1))/p.kT);

		% solve for a
		dydt = k_plus(i-1)*(1-a(i-1))*vS(i-1) - k_minus(i-1)*a(i-1);
		a(i) = dydt*1e-4 + a(i-1);
		if a(i) > 1
			a(i) = 1;
		elseif a(i) < 0
			a(i) = 0;
		end

		% now solve for F
		dydt = (1/p.tau_adap)*p.kT*(1/2 - a(i-1));
		F(i) = dydt*1e-4 + F(i-1);
		% if F(i) < 0
		% 	F(i) = 0;
		% end

	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a,time);
	F = interp1(T,F,time);
	H = interp1(T,H,time);
	k_minus = interp1(T,k_minus,time);
	k_plus = interp1(T,k_plus,time);
	k_D = k_minus./k_plus;

	% pass through a filter 
	t = 0:(length(S)/10);
	K = generate_simple_filter(p.K_tau,1,t);
	K = K(1:find(K>1e-2*max(K),1,'last'));
	R = filter(K,1,a);

		function f = generate_simple_filter(tau,n,t)
			f = t.^n.*exp(-t/tau); % functional form in paper
			f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
		end

end
