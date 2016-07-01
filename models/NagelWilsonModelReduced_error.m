%% NagelWilsonReduced.m
% reduced version of the model in Nagel and Wilson
% this uses a conservation equation to eliminate one of the reaction ODEs
% 
function [C,D,R] = NagelWilsonModelReduced_error(S,p)

	% list parameters for legibility
	p.theta;
	p.ka;
	p.kb;
	p.sa;
	p.sb;
	p.k0;
	p.A;
	p.B;


	% define initial conditions for ode system
	ic = zeros(5,1);
	ic(1:3) = .25;
	ic(4) = 0.4;
	ic(5) = .4; 


	time = 1e-3*(1:length(S));
	Tspan = [min(time) max(time)];

	options = odeset('MaxStep',.1);
	[T, Y] = ode23t(@(t,y) NagelWilsonModelReduced_error_ode(t,y,time,S,p),Tspan,ic,options); % Solve ODE

	% re-interpolate the solution to fit the stimulus
	C = interp1(T,Y(:,4),time);
	D = interp1(T,Y(:,5),time);

	temp = zeros(length(S),3);
	for i = 1:3
		temp(:,i) = interp1(T,Y(:,i),time);
	end

	% recover the OR species and insert it
	R = zeros(length(S),4);
	R(:,1:2) = temp(:,1:2);
	R(:,4) = temp(:,3);
	R(:,3) = 1- sum(temp,2);

	function dy = NagelWilsonModelReduced_error_ode(t,y,time,odor,p)

		dy = NaN*y;

		% calculate the odor at the time point
		O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

		% receptor-ligand binding

		R = y(1);
		RX = y(2);
		ORX = y(3);

		dR = -R*(p.ka*p.sa + O*p.kb*p.sb) + RX*p.sa + p.sb*(1 - R - RX - ORX);
		dRX = -R*(p.ka*p.sa) - RX*(p.sa + O*p.theta*p.kb*p.sb) + ORX*p.sb;
		dORX = RX*O*p.theta*p.kb*p.sb + (1 - R - RX - ORX)*p.theta*p.ka*p.sa - ORX*(p.sa+p.sb);

		dy(1:3) = [dR dRX dORX];
		  
		% diffusible factor and channel opening
		dy(4) = (RX+ORX)*p.ko*(1/y(5))*(1-y(4)) - p.kc*y(4); % according to Nagel and Wilson
		dy(5) = p.A*y(4) - p.B*y(5);
	end

end



