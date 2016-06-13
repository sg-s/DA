%% LFPmodel2.m
% a simple LFP model that slows down as well as Weber-Fechner gain scaling
% this is modified from LFPmodel as follows:
% instead of a logarithmic activation function, we use a DA-style activation model.

function R = LFPmodel2(S,p)

	% list parameters for legibility
	p.A;
	p.B;
	p.ko;
	p.tau;

	% filter parameters
	p.tau_y;

	ic = .1;

	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	% compute the filter
	t = 0:3000;
	K = t.^n.*exp(-t/tau); 
	K = K/tau^(n+1)/gamma(n+1); % normalize appropriately

	options = odeset('MaxStep',.1);
	[T, Y] = ode23t(@(t,y) LFPmodel2_ode(t,y,time,S,p),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	R = interp1(T,Y,time);

		function dy = LFPmodel2_ode(t,y,time,odor,p)
			% calculate the odor at the time point
			O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

			% cal

			eff_tau = p.tau*(1 + p.A*O)/(1 + p.B*O);
			dy = p.ko + log(O) - y;
			dy = dy/eff_tau;
		end
end