%% SchulzeLouisModel.m
% model specified in this paper:
% https://www.ncbi.nlm.nih.gov/pubmed/26077825
% http://dx.doi.org/10.7554/eLife.06694
% https://elifesciences.org/content/4/e06694



function [R,U] = SchulzeLouisModel(S,p)


	if size(S,2) > 1
		R = NaN*S;
		U = NaN*S;
		for i = 1:size(S,2)
			[R(:,i),U(:,i)] = SchulzeLouisModel(S(:,i),p);
		end

		return
	end

	% list parameters for legibility
	p.a1;
	p.a2;
	p.a3;
	p.b1;
	p.b2;
	p.b3;
	p.b4;
	p.theta;
	p.b5;
	p.n;

	% some bounds
	ub.n = 2;
	lb.n = 2;

	lb.a1 = 0;
	lb.a2 = 0;
	lb.a3 = 0;
	lb.b1 = 0;
	lb.b2 = 0;
	lb.b3 = 0;
	lb.b4 = 0;
	lb.theta = 0;
	lb.b5 = 0;

	ic = [.1 .1];

	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	options = odeset('MaxStep',.1);
	[T, Y] = ode23t(@(t,y) SchulzeLouisModel_ode(t,y,time,S,p),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	R = interp1(T,Y(:,2),time);
	U = interp1(T,Y(:,1),time);

		function dY = SchulzeLouisModel_ode(t,Y,time,odor,p)
			% calculate the odor at the time point
			O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

			u = Y(1);
			y = Y(2);

			dY = 0*Y;
			dY(1) = p.a1*O - p.a2*u + p.a3*y;
			dY(2) = p.b1*O/(p.b2 + O + p.b3*u) - (p.b4*(y^p.n))/(p.theta^p.n + y^p.n) - p.b5*y;
		end
end