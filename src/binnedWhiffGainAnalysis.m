% this performs an analysis on time series of stimulus, lfp and response from nat. stim data
% the idea is to show fast gain control in a model-free manner

function [s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(S,X,R,varargin)


% options and defaults
options.history_length = 300;
options.bin_size = .1; % [90% - 110%]
options.min_n = 4; 
options.min_shat_range = .9;
options.whiff_height_frac = .01;

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

for j = 1:size(S,2)
	ws(j) = whiffStatistics(S(:,j),X(:,j),R(:,j),options.history_length,'MinPeakProminence',max(S(:,j)*options.whiff_height_frac),'debug',false);
end

pS = vertcat(ws.stim_peaks);
last_stim_peak = circshift(pS,1);
time_to_last_stim_peak = vertcat(ws.stim_peak_loc) - circshift(vertcat(ws.stim_peak_loc),1);
pX = vertcat(ws.peak_LFP);
pR = vertcat(ws.peak_firing_rate);
Shat = vertcat(ws.stim_history_length_before_whiff);


n = NaN*vertcat(pS);
for i = 1:length(pS)
	n(i) = sum(pS > pS(i)*(1-options.bin_size) & pS < pS(i)*(1+options.bin_size));
end

% make a plot of rank order of responses vs. mean stim in preceding X ms

s = []; x = []; r = []; whiff_s = []; group = [];
c =1;
for i = 1:length(pS)
	if n(i) > options.min_n - 1

		these = (pS > pS(i)*(1-options.bin_size) & pS < pS(i)*(1+options.bin_size));
		stim = Shat(these);
		whiffs = pS(these);
		lfp = -pX(these);
		firing = pR(these);

		if max(stim) - min(stim) > options.min_shat_range*(max(Shat)- min(Shat))
			group = [group; ones(length(stim),1)*c];
			s = [s; stim(:)];
			x = [x; lfp/nanmean(lfp)];
			r = [r; firing/nanmean(firing)];
			whiff_s = [whiff_s; whiffs(:)];
			c = c + 1;
		end
	end
end
rm_this = isnan(s) | isnan(x);
s(rm_this) = []; x(rm_this) = []; r(rm_this) = []; whiff_s(rm_this) = []; group(rm_this) = [];