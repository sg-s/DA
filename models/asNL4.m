%% asNL.m
% NL model that slows down as it adapts
% LFP model built for fitting to real data

function [R,a,k1,k2,k_D,Shat] = asNL4(S,p)

	% bounds
	lb.A = 0;
	lb.B = 0;
	lb.adap_tau = 1;
	lb.K_tau = 1;
	lb.KD_min = eps;

	ub.adap_tau = 1e4;
	ub.K_tau = 1e3;


	% list parameters for legibility
	p.KD_min;

	% parameters for k1, k2 estimation
	p.A;
	p.B;

	% adaptation parameters
	p.adap_tau;

	% output filter
	p.K_tau;
	p.R_scale;


	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S; k1 = NaN*S; k2 = NaN*S; k_D = NaN*S; Shat = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),k1(:,i),k2(:,i),k_D(:,i),Shat(:,i)] = asNL_euler(S(:,i),p);
		end
		return
	end


	% compute the adaptation filter
	K_adap = ones(floor(1+p.adap_tau),1);
	

	% filter the stimulus with the adaptation filter
	Shat = filter(K_adap,sum(K_adap),S);

	% compute the rate constants 
	% k2 = p.A./Shat;
	% k1 = k2./Shat;

	% other form
	k1 = p.A*exp(-Shat./p.B);
	k2 = k1.*Shat;

	k_D = k2./k1;

	low_k_D = k_D < p.KD_min;
	k2(low_k_D) = k1(low_k_D)*p.KD_min;
	k_D(low_k_D) = p.KD_min;


	time = 1e-3*(1:length(S));
	T = 1e-4:1e-4:max(time);

	% interpolate 
	vS = interp1(time,S,T); vS(isnan(vS)) = S(1);
	vk1 = interp1(time,k1,T); vk1(isnan(vk1)) = k1(1);
	vk2 = interp1(time,k2,T); vk2(isnan(vk2)) = k2(1);

	% use a fixed-step Euler to solve this
	a = .1+0*vS;
	for i = 2:length(vS)
		dydt = vk1(i-1)*(1-a(i-1))*vS(i-1) - vk2(i-1)*a(i-1);
		a(i) = dydt*1e-4 + a(i-1);
		if a(i) > 1
			a(i) = 1;
		elseif a(i) < 0
			a(i) = 0;
		end

	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a,time);

	% pass through a filter 
	t = 0:(length(S)/10);
	K = generate_simple_filter(p.K_tau,1,t);
	K = K(1:find(K>1e-2*max(K),1,'last'));
	R = p.R_scale*(filter(K,1,a));


		function f = generate_simple_filter(tau,n,t)
			f = t.^n.*exp(-t/tau); % functional form in paper
			f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
		end

end
