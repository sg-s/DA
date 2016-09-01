% consolidateData2.m
% consolidate data accepts a path, and merges all the .mat files there into one consolidate blob
% it makes many assumptions on how the data should be: it should be a Kontroller format
% and should have PID, voltage.
% this is not very general purpose, but hopefully will become more hardened to edge cases as time progresses
% 
% created by Srinivas Gorur-Shandilya at 11:09 , 06 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function consolidated_data = consolidateData2(pathname,varargin)

if ~nargin
	pathname = pwd;
end


% options and defaults
options.use_cache = true;
options.fA = true;
options.PID = true;
options.fB = true;
options.LFP = true;
options.fly = true;
options.orn = true;
options.AllControlParadigms = true;
options.paradigm_hashes = true;
options.max_length = 70e3; % for firing rate, etc.

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


fA = [];
fB = [];
PID = [];
fly = [];
paradigm = [];
A_spikes = [];
B_spikes = [];
A_spike_amp = [];
B_spike_amp = [];
orn = [];
AllControlParadigms = struct;
AllControlParadigms.Name = '';
AllControlParadigms.Outputs = [];
AllControlParadigms(1) = [];
paradigm_hashes = {};
sequence = []; % for each file. 

if ~strcmp(pathname(end),oss)
	pathname = [pathname oss];
end

allfiles = dir([pathname '*.mat']);
% remove the consolidated data from this
rm_this = [find(strcmp('cached_log.mat',{allfiles.name})) find(strcmp('consolidated_data.mat',{allfiles.name})) find(strcmp('cached.mat',{allfiles.name})) find(strcmp('template.mat',{allfiles.name}))];
if ~isempty(rm_this)
	allfiles(rm_this) = [];
end

% always be caching. 
hash = dataHash(allfiles);
if options.use_cache
	consolidated_data = [];
	try
		consolidated_data = load([pathname 'consolidated_data.mat'],'consolidated_data');
		consolidated_data = consolidated_data.consolidated_data;
	catch
	end
	if ~isempty(consolidated_data)
		return
	end
end

% find the longest trial in all the data. this will be the length of the combined data
disp('Determining longest data length...')
ll = 0;
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	for j = 1:length(data)
		if ~isempty(data(j).voltage)
			ll = max([ll length(data(j).voltage)]);
		end
	end
end

if ll > options.max_length*10
	ll = options.max_length*10;
end

disp('Longest trial observed is:')
disp(ll)
ll = ceil(ll/10);

% make placeholders
LFP = zeros(ll,0);
PID = zeros(ll,0);
fA = zeros(ll,0);
fB = zeros(ll,0);
A_spikes = sparse(ll*10,0);
B_spikes = sparse(ll*10,0);
A_spike_amp = sparse(ll*10,0);
B_spike_amp = sparse(ll*10,0);


for i = 1:length(allfiles)
	clear spikes data metadata
	load(strcat(pathname,allfiles(i).name));
	disp(strcat(pathname,allfiles(i).name));
	for j = 1:length(data)
		clear this_PID this_LFP this_A_spikes this_B_spikes this_hash this_paradigm 
		if ~isempty(data(j).voltage)

			disp(['Paradigm : ' mat2str(j)]) 
			% figure out which control paradigm this is
			this_hash = dataHash(ControlParadigm(j));
			if isempty(find(strcmp(this_hash,paradigm_hashes)))
				AllControlParadigms(end+1) = orderfields(ControlParadigm(j));
				this_paradigm = length(AllControlParadigms);
				paradigm_hashes{end+1} = this_hash;
			else
				this_paradigm = find(strcmp(this_hash,paradigm_hashes));
			end

			% intercept early to get max length correctly
			if length(data(j).voltage) > ll*10

				data(j).voltage = data(j).voltage(:,1:ll*10);
				data(j).PID = data(j).PID(:,1:ll*10);
				try
					spikes(j).A = spikes(j).A(:,1:ll*10);
					spikes(j).B = spikes(j).B(:,1:ll*10);
					spikes(j).amplitudes_A = spikes(j).A(:,1:ll*10);
					spikes(j).amplitudes_B = spikes(j).B(:,1:ll*10);
					
				catch
					
				end
			end

			disp('Loading PID and raw voltage...')
			this_LFP = data(j).voltage;
			this_PID = data(j).PID;

			this_fA = NaN*this_LFP(:,1:10:end)';
			this_fB = NaN*this_LFP(:,1:10:end)';

			disp(['There appear to be ' mat2str(width(this_PID)) ' trials for this paradigm.'])

			if exist('spikes','var')
				disp('spikes variable exists for this file.')
				if length(spikes) < j
					disp('No spike info for this paradigm.')
					this_A_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
					this_B_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
					this_A_spike_amp = sparse(size(this_LFP,1),size(this_LFP,2));
					this_B_spike_amp = sparse(size(this_LFP,1),size(this_LFP,2));
				else
					if length(spikes(j).A) > 10 && max(max(spikes(j).A)) > 0
						disp('A Spikes exist for this paradigm, loading...')
						this_A_spikes = spikes(j).A;
						this_A_spike_amp = spikes(j).amplitudes_A;

						if width(this_A_spikes) < width(this_PID)
							disp('Spike width less than data width, padding...')
							this_A_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
							this_A_spike_amp = sparse(size(this_LFP,1),size(this_LFP,2));
							for k = 1:width(this_A_spikes)
								try
									this_A_spikes(k,:) = spikes(j).A(k,:);
									this_A_spike_amp(k,:) = spikes(j).amplitudes_A(k,:);
								end
							end
						end
					else
						disp('No A spikes found!')
						% create dummy this_A_spikes
						this_A_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
						this_A_spike_amp = sparse(size(this_LFP,1),size(this_LFP,2));
					end

					% grab B spikes
					if length(spikes(j).B) > 10 && max(max(spikes(j).B)) > 0
						disp('B Spikes exist for this paradigm, loading...')
						this_B_spikes = spikes(j).B;
						this_B_spike_amp = spikes(j).amplitudes_B;
						if width(this_B_spikes) < width(this_PID)
							disp('Spike width less than data width, padding...')
							this_B_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
							this_B_spike_amp = sparse(size(this_LFP,1),size(this_LFP,2));
							for k = 1:width(this_B_spikes)
								try
									this_B_spikes(k,:) = spikes(j).B(k,:);
									this_B_spike_amp(k,:) = spikes(j).amplitudes_B(k,:);
								catch
								end
							end
						end
					else
						disp('No B spikes found!')
						% create dummy this_B_spikes
						this_B_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
						this_B_spike_amp = sparse(size(this_LFP,1),size(this_LFP,2));
					end
				


					% censor trace use use_trace_fragment
					use_trace_fragment = [];
					try
						use_trace_fragment = spikes(j).use_trace_fragment;
					catch
					end

					if ~isempty(use_trace_fragment)
						disp('Censoring parts of trial based on annotation...')
						if width(use_trace_fragment) == width(this_PID)
							this_LFP(~logical(use_trace_fragment)) = NaN;
						else
							for k = 1:width(use_trace_fragment)
								this_LFP(k,~logical(use_trace_fragment(k,:))) = NaN;
							end
						end
					end
				end
			end

			% convert to firing rate
			this_fA = spiketimes2f(this_A_spikes,1e-4*(1:length(this_A_spikes)),1e-3,3e-2);
			this_fB = spiketimes2f(this_B_spikes,1e-4*(1:length(this_B_spikes)),1e-3,3e-2);

			rm_this = [];
			try
				rm_this = find(spikes(j).discard);
			catch
			end
			if ~isempty(rm_this)
				disp('Discarding some trials acc. to annotation:')
				this_PID(rm_this,:) = [];
				this_LFP(rm_this,:) = [];
				this_A_spikes(rm_this,:) = [];
				this_B_spikes(rm_this,:) = [];
				this_A_spike_amp(rm_this,:) = [];
				this_B_spike_amp(rm_this,:) = [];
				this_fA(:,rm_this) = [];
				this_fB(:,rm_this) = [];

			else
				disp('No trials to discard')
			end

			disp('Downsampling data...')
			this_PID = this_PID(:,1:10:end)';
			this_LFP = this_LFP(:,1:10:end)';


			if length(this_LFP) < size(LFP,1) && length(this_LFP) > 1% this is to account for shorter paradigms, and pad them
				disp('Short paradgim, need to pad...')
				padding = NaN(size(LFP,1)-length(this_LFP),width(this_LFP));
				this_LFP = [this_LFP; padding];
				padding = NaN(size(PID,1)-length(this_PID),width(this_PID));
				this_PID = [this_PID; padding];
				padding = NaN(size(fA,1)-length(this_fA),width(this_fA));
				try
					this_fA = [this_fA; padding];
					this_fB = [this_fB; padding];
				catch
					this_fA = [this_fA(:); padding];
					this_fB = [this_fB(:); padding];
				end

				% pad spikes too
				disp('Padding spikes...')
				
				padding = NaN(size(A_spikes,1)-length(this_A_spikes),width(this_A_spikes));
				this_A_spikes = [this_A_spikes'; padding];
				this_A_spike_amp = [this_A_spike_amp'; padding];

				padding = NaN(size(B_spikes,1)-length(this_B_spikes),width(this_B_spikes));
				this_B_spikes = [this_B_spikes'; padding];
				this_B_spike_amp = [this_B_spike_amp'; padding];
			else
				disp('No padding, rotating spike matrix...')
				% still need to pad spikes because of stupidity (the spike lengths may be slightly off because we divide the longest length by 10, then multiply by 10 again)
				padding = NaN(size(A_spikes,1)-length(this_A_spikes),width(this_A_spikes));
				this_A_spikes = [this_A_spikes'; padding];
				this_A_spike_amp = [this_A_spike_amp'; padding];

				padding = NaN(size(B_spikes,1)-length(this_B_spikes),width(this_B_spikes));
				this_B_spikes = [this_B_spikes'; padding];
				this_B_spike_amp = [this_B_spike_amp'; padding];
			end

			% check that PID and fA match
			if length(unique([width(this_PID) width(this_fA) width(this_LFP) width(this_A_spikes)])) > 1
				disp('Non matching PID and fA widths')
				keyboard
			end

			if length(unique([width(PID) width(fA) width(LFP) width(A_spikes)])) > 1
				disp('Non matching all_PID and fA widths')
				keyboard
			end

			% consolidate
			try
				LFP = [LFP this_LFP];
				PID = [PID this_PID];
				fA =  [fA  this_fA ];
				fB =  [fB  this_fB ];
				A_spikes = [A_spikes this_A_spikes];
				B_spikes = [B_spikes this_B_spikes];
				A_spike_amp = [A_spike_amp this_A_spike_amp];
				B_spike_amp = [B_spike_amp this_B_spike_amp];
			catch er
				disp(er)
				keyboard
				
			end
			

			paradigm = [paradigm  this_paradigm*ones(1,width(this_PID))];
			orn = [orn  i*ones(1,width(this_PID))];

			% figure out the fly # from the file name
			try
				us = strfind(allfiles(i).name,'_'); % underscores
				this_fly = str2double(allfiles(i).name(strfind(allfiles(i).name,'_F')+2:us(find(us>strfind(allfiles(i).name,'_F')+2,1,'first'))-1));
				this_fly = 100*str2double(strrep(allfiles(i).name(1:us(3)-1),'_','')) + this_fly;
				fly = [fly; this_fly*ones(width(this_PID),1)];
			catch
				warning('Error in determining fly ID. Fly IDs may not match data.')
			end

			% also add the sequence 
			try
				this_sequence = find(timestamps(1,:)==j);
				if length(rm_this) 
					this_sequence(rm_this) = [];	
				end
				sequence = [sequence this_sequence];
			catch
				warning('Error in assembling sequence info.')
			end
			

		end
	end
end


% cache all of this
consolidated_data = struct;
consolidated_data.PID = PID;
consolidated_data.LFP = LFP;
consolidated_data.fA = fA;
consolidated_data.fB = fB;
consolidated_data.paradigm = paradigm;
consolidated_data.AllControlParadigms = AllControlParadigms;
consolidated_data.paradigm_hashes = paradigm_hashes;
consolidated_data.orn = orn;
consolidated_data.fly = fly;
consolidated_data.A_spikes  = A_spikes;
consolidated_data.B_spikes  = B_spikes;
consolidated_data.A_spike_amp  = A_spike_amp;
consolidated_data.B_spike_amp  = B_spike_amp;
consolidated_data.sequence = sequence;
save([pathname 'consolidated_data.mat'],'consolidated_data');

