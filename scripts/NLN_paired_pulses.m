
pHeader;

%% Can NLN models fool paired pulse protocols?
% In this document, I check if NLN models can fool paired pulse protocols, like the one used by Cao et al., that are used to check for the timescale of adaptation in neurons. 


%% To do so, I generate responses using a NLN model, using paired pulses as my stimulus, while varying the spacing between the pulse pairs. 

% stimulus
pulse_sep = ceil(logspace(2.5,4,10)); % milliseconds
T = 12e3;
R = NaN(T,length(pulse_sep));
first_pulse = 1e3;

% parameters for input nonlinearity 
hill_param = [1 .5 1];


% construct a filter
clear p
p.n = 2; p.tau1 = 20; p.tau2 = 100; p.A = 0.3;
K = filter_gamma2(1:1e3,p);

for i = 1:length(pulse_sep)
	% construct stimulus
	S = zeros(T,1);
	S(first_pulse) = 1;
	S(first_pulse+pulse_sep(i)) = 1;
	S = filter(ones(50,1),1,S);

	% pass it through the model
	x = hill(hill_param,S);
	R(:,i) = filter(K,1,x);

end


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


