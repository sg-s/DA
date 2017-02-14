% determines the cutoff in the response when the model deviations from the response are uncorrelated with the instantaneous stimulus
% meant to be used with plot_tau_gain_nat_stim

function rc = findResponseCutoffTauGainPlot(S,R,P,use_spear)


S = S(:);
R = R(:);
P = P(:);

minrc = min(R);
maxrc = max(R);

allrc = linspace(minrc,maxrc,100);
rho = NaN*allrc;

for i = 1:100
	D = R - P;
	rm_this = isnan(D) | isnan(S) | R < allrc(i);
	x = S(~rm_this);
	y = D(~rm_this);
	if use_spear
		rho(i) = corr(x(1:10:end),y(1:10:end),'type','Spearman');
	else
		rho(i) = corr(x(1:10:end),y(1:10:end),'type','Pearson');
	end
end

rc = allrc(find(rho>0,1,'first'));