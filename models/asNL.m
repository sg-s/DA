%% asNL.m
% NL model that slows down as it adapts

function [R,a,k1,k2,k_D,Shat] = asNL(S,p)

	% list parameters for legibility

	% parameters for k1, k2 estimation
	p.A;

	% adaptation parameters
	p.adap_tau;

	% output filter
	p.K_A;
	p.K_n;
	p.K_tau;

	% compute the adaptation filter
	K_adap = ones(floor(1+p.adap_tau),1);
	

	% filter the stimulus with the adaptation filter
	Shat = filter(K_adap,sum(K_adap),S);

	% compute the rate constants 
	k2 = p.A./Shat;
	k1 = k2./Shat;

	k_D = k2./k1;

	ic = 0;

	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	options = odeset('MaxStep',.01);
	[T, Y] = ode23t(@(t,y) asNL_ode(t,y,time,S,k1,k2),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	a = interp1(T,Y,time);

	% pass through a filter 
	t = 0:(length(S)/10);
	K = generate_simple_filter(p.K_tau,p.K_n,t);
	K = K(1:find(K>1e-2*max(K),1,'last'));
	R = filter(K,1,a);

		function f = generate_simple_filter(tau,n,t)
			f = t.^n.*exp(-t/tau); % functional form in paper
			f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
		end

		function dy = asNL_ode(t,y,time,S,k1,k2)
			% calculate the stimulus at the time point
			S_now = interp1(time,S,t); % Interpolate the data set (ft,f) at time t\
			k1_now = interp1(time,k1,t);
			k2_now = interp1(time,k2,t); 

			dy = k1_now*(1-y)*S_now - k2_now*y; 
		end
end
