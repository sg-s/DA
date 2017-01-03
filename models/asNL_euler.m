%% asNL.m
% NL model that slows down as it adapts
% euler implementation to avoid solving ODEs

function [R,a,k1,k2,k_D,Shat] = asNL_euler(S,p)

	% list parameters for legibility

	% parameters for k1, k2 estimation
	p.A;

	% adaptation parameters
	p.adap_tau;
	p.adap_lag;

	% output filter
	p.K_n;
	p.K_tau;

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
	Shat = circshift(Shat,floor(p.adap_lag));

	% compute the rate constants 
	k2 = p.A./Shat;
	k1 = k2./Shat;

	k_D = k2./k1;

	a = 0.5+0*S;

	for i = 2:length(S)
		dydt = k1(i-1)*(1-a(i-1))*S(i-1) - k2(i-1)*a(i-1);
		a(i) = dydt*1e-3 + a(i);
		if a(i) > 1
			a(i) = 1;
		end
		if a(i) < 0
			a(i) = 0;
		end
	end


	% pass through a filter 
	t = 0:(length(S)/10);
	K = generate_simple_filter(p.K_tau,p.K_n,t);
	K = K(1:find(K>1e-2*max(K),1,'last'));
	R = filter(K,1,a);

		function f = generate_simple_filter(tau,n,t)
			f = t.^n.*exp(-t/tau); % functional form in paper
			f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
		end

end
