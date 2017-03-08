%% ORNmodel
% this is the start of a new class of models, where we get k_D to adapt, but not perfectly, so that a(t) doesn't adapt perfectly, and we can use that to speed up responses later on

function [R,a,b,k_D,tau] = ORNmodel(S,p)

	% list parameters for legibility
	p.B;
	p.K_a;

	% LFP -> firing nonlinearity 
	p.n_firing;
	p.kD_firing;
	p.A_firing;

	p.tau_tau; % timescale over which the timescale of the filter is updated 
	p.tau_0; 

	% bounds
	lb.A_firing = 0;
	lb.kD_firing = 0;
	lb.n_firing = 1;



	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S; b = a; k_D = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),b(:,i),k_D(:,i)] = ORNmodel(S(:,i),p);
		end
		return
	end



	time = 1e-3*(1:length(S));
	T = 1e-4:1e-4:max(time);

	vS = interp1(time,S,T); vS(isnan(vS)) = S(1);

	% use a fixed-step Euler to solve this
	a = 0*vS;
	vkD = 0*vS + vS(1);

	
	for i = 2:length(vS)

		a(i) = vS(i-1)/(vS(i-1) + vkD(i-1));

		fkD = (1/2)*(vkD(i-1))/(vkD(i-1) + p.K_a);

		% also change vkD
		dydt = p.B*vkD(i-1)*(a(i-1) - fkD);
		vkD(i) = dydt*1e-4 + vkD(i-1);
		if vkD(i) < 0
			vkD(i) = 0;
		end
	end


	% re-interpolate the solution to fit the stimulus
	a = interp1(T,a,time);
	k_D = interp1(T,vkD,time);


	% pass this through a nonlinearity 
	b = a.^p.n_firing;
	b = b./(b + p.kD_firing.^p.n_firing);

	% now we need to filter it, with timescales 
	tau_tau = 1+floor(p.tau_tau);
	tau = floor(p.tau_0*exp(-filter(ones(tau_tau,1),1,a)));

	% filter by solving a ODE that takes a derivative 
	R = b;
	for i = 2:length(b)
		if i - tau(i) > 1
			dydt = b(i) - b(i-tau(i));
		else
			dydt = 0;
		end
		R(i) = R(i-1) + dydt;

	end


end
