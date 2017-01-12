%% LFPmodelv1
% in this version, there is actually no adaptation
% k1, and k2 are fixed parameters 

function [R,a,k_D] = LFPmodelv1(S,p)

	% list parameters for legibility

	p.k1;
	p.k2;

	% trivial parameters
	p.R_scale;
	p.R_offset;


	% bounds
	lb.k1 = 0;
	lb.k2 = 0;



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
	

	k2 = p.k2;
	k1 = p.k1;

	k_D = k1./k2;

	for i = 2:length(vS)
		
		dydt = k1*(1-a(i-1))*vS(i-1) - k2*a(i-1);
		a(i) = dydt*1e-4 + a(i-1);
		if a(i) > 1
			a(i) = 1;
		elseif a(i) < 0
			a(i) = 0;
		end


	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a,time);

	R = (p.R_offset + a)*p.R_scale;



end
