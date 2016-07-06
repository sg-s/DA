%% LFPmodel.m
% a simple LFP model that slows down as well as Weber-Fechner gain scaling

function R = LFPModelForFitting(S,p)

	% list parameters for legibility
	p.A;
	p.B;
	p.ko;
	p.tau;
	p.scale1;
	p.scale0;

	ic = .1;

	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	options = odeset('MaxStep',.1);
	[T, Y] = ode23t(@(t,y) LFPmodel_ode(t,y,time,S,p),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	R = p.scale0 + p.scale1*interp1(T,Y,time);


		function dy = LFPmodel_ode(t,y,time,odor,p)
			% calculate the odor at the time point
			O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t
			eff_tau = p.tau*(1 + p.A*O)/(1 + p.B*O);
			dy = p.ko + log(O) - y;
			dy = dy/eff_tau;
		end
end