%% asNL.m
% NL model that slows down as it adapts

function [R,a,k1,k2,k_D,Shat] = asNL(S,p)

	% list parameters for legibility

	% parameters for k1, k2 estimation
	p.A;

	% adaptation parameters
	p.adap_tau;

	% output filter
	p.K_n;
	p.K_tau;

	% work with matrices too
	if size(S,2) > 1
		R = NaN*S; a = NaN*S; k1 = NaN*S; k2 = NaN*S; k_D = NaN*S; Shat = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),a(:,i),k1(:,i),k2(:,i),k_D(:,i),Shat(:,i)] = asNL(S(:,i),p);
		end
		return
	end


	% compute the adaptation filter
	K_adap = ones(floor(1+p.adap_tau),1);
	

	% filter the stimulus with the adaptation filter
	Shat = filter(K_adap,sum(K_adap),S);

	% compute the rate constants 
	k2 = p.A./Shat;
	k1 = k2./Shat;

	k_D = k2./k1;

	ic = 0;

	time = linspace(0,length(S)*1e-3,length(S));
	Tspan = [min(time) max(time)];

	options = odeset('InitialStep',1e-3,'MaxStep',1,'RelTol',1e-3);
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
			S_now = interp1q(time,S,t); % Interpolate the data set (ft,f) at time t\
			k1_now = interp1q(time,k1,t);
			k2_now = interp1q(time,k2,t); 

			dy = k1_now*(1-y)*S_now - k2_now*y; 
		end

end
