%% LFPmodelv4
% adapting model, k_D changes like in chemotaxis, k1 fixed
% we also impose a minimum k_D 

function [R,a,k_D] = LFPmodelv4(S,p)

	% list parameters for legibility

	p.k1;
	p.adap_tau;
	p.k_min;

	% trivial parameters
	p.R_scale;
	p.R_offset;

	% bounds
	lb.k1 = 0;
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
	

	k1 = p.k1;

	for i = 2:length(vS)

		k2 = k1*k_D(i-1);
		
		dydt = k1*(1-a(i-1))*vS(i-1) - k2*a(i-1);
		a(i) = dydt*1e-4 + a(i-1);
		if a(i) > 1
			a(i) = 1;
		elseif a(i) < 0
			a(i) = 0;
		end

		% also change k_D
		dydt = (1/p.adap_tau)*k_D(i-1)*(a(i-1) - 1/2);
		k_D(i) = dydt*1e-4 + k_D(i-1);
		if k_D(i) < p.k_min
			k_D(i) = p.k_min;
		end


	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a,time);

	R = (p.R_offset + a)*p.R_scale;



end