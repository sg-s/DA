function [ws] = whiffStatistics(S,X,R,history_length,varargin)


% options and defaults
options.MinPeakProminence = .25;
options.MinPeakDistance = 100;
options.debug = false;

if nargout && ~nargin 
	varargout{1} = options;
    return
end

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options.(temp) = varargin{ii+1};
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

Shat = filter(ones(300,1),history_length,S);


[ws.stim_peaks, ws.stim_peak_loc] = findpeaks(S,'MinPeakProminence',options.MinPeakProminence,'MinPeakDistance',options.MinPeakDistance);

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

if options.debug
	figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	ax1 = subplot(2,1,1); hold on
	plot(S,'k')
	plot(ws.stim_peak_loc,ws.stim_peaks,'ro')

	ax2 = subplot(2,1,2); hold on
	plot(R,'k')
	plot(ws.peak_firing_locs,ws.peak_firing_rate,'ro')
	linkaxes([ax1 ax2],'x')
end

ws.stim_history_length_before_whiff =  [NaN; Shat(nonnans(ws.min_stim_before_whiff_locs)); NaN];