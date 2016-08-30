%% cleanMSGdata
% helper function for webers_generally_observed

function cdata = cleanMSGdata(cdata,varargin)

% options and defaults
options.extract_filter = true;
options.use_cache = true;

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
    		options = setfield(options,temp,varargin{ii+1});
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

% unpack
v2struct(cdata);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces
% bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
% PID(:,bad_trials) = [];
% LFP(:,bad_trials) = [];
% fA(:,bad_trials) = [];
% fB(:,bad_trials) = [];
% paradigm(bad_trials) = [];
% fly(bad_trials) = [];
% orn(bad_trials) = [];

% save the original LFP
raw_LFP = 10*LFP;
% filter raw LFP to remove spikes, and remove baseline
for i = 1:width(raw_LFP)
	raw_LFP(:,i) = filtfilt(ones(30,1),30,raw_LFP(:,i));
	raw_LFP(:,i) = raw_LFP(:,i) - nanmean(raw_LFP(1:5e3,i));
end

% filter the LFP
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,LFP(:,i));
	LFP(:,i) = LFP(:,i)*10; % to get the units right, now in mV
end

if ~exist('a','var')
	a = 35e3; z = 55e3;
end

% extract filters 
K1 = []; K2 = []; fA_pred = []; LFP_pred = []; LFP_gain = []; fA_gain = [];
if options.extract_filter
	[K2,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',options.use_cache,'a',a,'z',z);
	[K1,LFP_pred,LFP_gain] = extractFilters(PID,LFP,'use_cache',options.use_cache,'a',a,'z',z);
end



% repack everything
clear cdata
cdata = v2struct(a,z,K1,K2,PID,LFP,raw_LFP,fA,fB,LFP_pred,fA_pred,fA_gain,LFP_gain,orn,fly,paradigm,AllControlParadigms);
