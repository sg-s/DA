%% LFPmodelv6
% in this version, K_D adapts through a feed-forward mechanism
% both k1 and k2 decrease with k_D:
% k1 ~ 1/k_D^2, k2 ~ 1/k_D

function [R,a,k_D] = LFPmodelv6(S,p)

	% list parameters for legibility

	p.kappa1;
	p.kappa_max;
	p.kappa_min;
	p.adap_tau;
	p.k_min;

	% trivial parameters
	p.R_scale;
	p.R_offset;

	% bounds
	lb.kappa1 = 0;
	lb.kappa_max = 0;
	lb.kappa_min = 0;
	lb.adap_tau = 0;
	lb.k_min = 1e-6;

	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S; k1 = NaN*S; k2 = NaN*S; k_D = NaN*S; Shat = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),k1(:,i),k2(:,i),k_D(:,i),Shat(:,i)] = adaptingLFPmodel(S(:,i),p);
		end
		return
	end



	time = 1e-3*(1:length(S));
	T = 1e-4:1e-4:max(time);

	vS = interp1(time,S,T); vS(isnan(vS)) = S(1);

	% use a fixed-step Euler to solve this
	a = 0*vS;
	k_D = p.k_min + 0*vS;
	

	for i = 2:length(vS)

		k2 = p.kappa_max*(1 - ((k_D(i-1))/(k_D(i-1) + p.kappa1))) + p.kappa_min;
		k1 = k2/k_D(i-1);

		dydt = k1*(1-a(i-1))*vS(i-1) - k2*a(i-1);
		a(i) = dydt*1e-4 + a(i-1);
		if a(i) > 1
			a(i) = 1;
		elseif a(i) < 0
			a(i) = 0;
		end

		% also change k_D
		dydt = (1/p.adap_tau)*(vS(i-1) - k_D(i-1));
		k_D(i) = dydt*1e-4 + k_D(i-1);
		if k_D(i) < p.k_min
			k_D(i) = p.k_min;
		end


	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a,time);

	R = (p.R_offset + a)*p.R_scale;



end
