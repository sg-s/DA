% computes contrast and gains for the variance switching experiment
% this is a helper script for modular analysis

% for each trial, compute the relative gain change and the input contrast change
input_contrast_change_K1p = NaN(width(reshaped_dLFP),1);
input_contrast_change_K2K1p = NaN(width(reshaped_dLFP),1);
input_contrast_change_K2p = NaN(width(reshaped_dLFP),1);
input_contrast_change_K3p = NaN(width(reshaped_dLFP),1);
gain_change_K1p = NaN(width(reshaped_dLFP),1);
gain_change_K2K1p = NaN(width(reshaped_dLFP),1);
gain_change_K2p = NaN(width(reshaped_dLFP),1);
gain_change_K3p = NaN(width(reshaped_dLFP),1);

for i = 1:width(reshaped_dLFP)
	if r2_K1p(i) > min_r2
		y = reshaped_dLFP(1e3:4e3,i);
		x = K1p(1e3:4e3,i);
		ff = fit(x(:),y(:),'poly1');
		hi_gain = ff.p1;
		y = reshaped_dLFP(6e3:9e3,i);
		x = K1p(6e3:9e3,i);
		ff = fit(x(:),y(:),'poly1');
		lo_gain = ff.p1;
		gain_change_K1p(i) = lo_gain./hi_gain;
		clear lo_gain hi_gain

		input_contrast_change_K1p(i) = std(reshaped_PID(6e3:9e3,i))/std(reshaped_PID(1e3:4e3,i));
	end

	if r2_K2K1p(i) > min_r2
		y = reshaped_fA(1e3:4e3,i);
		x = K2K1p(1e3:4e3,i);
		r = max(y) - min(y);
		l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
		ff = fit(x(l),y(l),'poly1');
		hi_gain = ff.p1;
		y = reshaped_fA(6e3:9e3,i);
		x = K2K1p(6e3:9e3,i);
		r = max(y) - min(y);
		l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
		ff = fit(x(l),y(l),'poly1');
		lo_gain = ff.p1;
		gain_change_K2K1p(i) = lo_gain./hi_gain;
		clear lo_gain hi_gain

		input_contrast_change_K2K1p(i) = std(K1p(6e3:9e3,i))/std(K1p(1e3:4e3,i));
	end

	if r2_K2p(i) > min_r2
		y = reshaped_fA(1e3:4e3,i);
		x = K2p(1e3:4e3,i);
		r = max(y) - min(y);
		l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
		ff = fit(x(l),y(l),'poly1');
		hi_gain = ff.p1;
		y = reshaped_fA(6e3:9e3,i);
		x = K2p(6e3:9e3,i);
		r = max(y) - min(y);
		l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
		ff = fit(x(l),y(l),'poly1');
		lo_gain = ff.p1;
		gain_change_K2p(i) = lo_gain./hi_gain;
		clear lo_gain hi_gain

		input_contrast_change_K2p(i) = std(reshaped_dLFP(6e3:9e3,i))/std(reshaped_dLFP(1e3:4e3,i));
	end


	if r2_K3p(i) > min_r2
		y = reshaped_fA(1e3:4e3,i);
		x = K3p(1e3:4e3,i);
		r = max(y) - min(y);
		l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
		ff = fit(x(l),y(l),'poly1');
		hi_gain = ff.p1;
		y = reshaped_fA(6e3:9e3,i);
		x = K3p(6e3:9e3,i);
		r = max(y) - min(y);
		l = (y > (r/3 + min(y))) & (y < (max(y) - r/3));
		ff = fit(x(l),y(l),'poly1');
		lo_gain = ff.p1;
		gain_change_K3p(i) = lo_gain./hi_gain;
		clear lo_gain hi_gain

		input_contrast_change_K3p(i) = std(reshaped_PID(6e3:9e3,i))/std(reshaped_PID(1e3:4e3,i));
	end
end