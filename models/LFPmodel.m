%% LFPmodel.m
% a simple LFP model that slows down as well as Weber-Fechner gain scaling

function R = LFPmodel(S,p)

	% list parameters for legibility
	p.A;
	p.B;
	p.C;
	p.ko;
	p.tau;

	ic = .1;

	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	options = odeset('MaxStep',.1);
	[T, Y] = ode23t(@(t,y) LFPmodel_ode(t,y,time,S,p),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	R = interp1(T,Y,time);

		function dy = LFPmodel_ode(t,y,time,odor,p)
			% calculate the odor at the time point
			O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t
			eff_tau = p.tau*(p.A + p.B*O)/(1 + p.C*O);
			dy = p.ko + log(O) - y;
			dy = dy/eff_tau;
		end
end