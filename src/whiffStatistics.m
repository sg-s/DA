function [ws] = whiffStatistics(S,X,R,history_length)


Shat = filter(ones(300,1),history_length,S);


[ws.stim_peaks, ws.stim_peak_loc] = findpeaks(S,'MinPeakProminence',.25,'MinPeakDistance',100);
ws.min_stim_before_whiff = NaN*ws.stim_peak_loc;
ws.min_stim_before_whiff_locs = NaN*ws.stim_peak_loc;
ws.peak_firing_rate = NaN*ws.stim_peaks;
ws.peak_LFP = NaN*ws.stim_peaks;
ws.peak_firing_locs = NaN*ws.stim_peaks;
ws.peak_LFP_locs = NaN*ws.stim_peaks;

for i = 2:(length(ws.stim_peak_loc)-1)
	a = max([ws.stim_peak_loc(i) - 90 ws.stim_peak_loc(i-1)]);
	[ws.min_stim_before_whiff(i), loc] = min(S(a:ws.stim_peak_loc(i)));
	ws.min_stim_before_whiff_locs(i) = loc + a;
end
for i = 2:(length(ws.stim_peak_loc)-1)
	z = min([ws.stim_peak_loc(i)+200 ws.min_stim_before_whiff_locs(i+1)]);
	[ws.peak_firing_rate(i),loc] = max(R(ws.stim_peak_loc(i):z));
	ws.peak_firing_locs(i) = loc + ws.stim_peak_loc(i);

	z = min([ws.stim_peak_loc(i)+200 ws.min_stim_before_whiff_locs(i+1)]);
	[ws.peak_LFP(i),loc] = min(X(ws.stim_peak_loc(i):z));
	ws.peak_LFP_locs(i) = loc + ws.stim_peak_loc(i);
end

ws.stim_history_length_before_whiff =  [NaN; Shat(nonnans(ws.min_stim_before_whiff_locs)); NaN];