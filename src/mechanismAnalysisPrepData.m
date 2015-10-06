% MechanismAnalysis_PrepData
% this function prepares data from the experiments investigating the mechanism of adaptation
% for a use example, look at Mechanism3.m
% 
% created by Srinivas Gorur-Shandilya at 10:22 , 21 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [stim,resp,ParadigmNames,paradigm] = mechanismAnalysisPrepData(datapath,haz_data,stim,resp,ParadigmNames,paradigm,use_light)

h = dataHash(whos);
cached_data = cache(h);
if isempty(cached_data)
else
	stim = cached_data.stim;
	resp = cached_data.resp;
	ParadigmNames = cached_data.ParadigmNames;
	paradigm = cached_data.paradigm;
	return
end

load(datapath);
if isempty(paradigm)
	p_offset = 0;
else
	p_offset = max(paradigm);
end

if isempty(haz_data)
	% need to figure out which paradigms to use on our own 
	haz_data = find(Kontroller_ntrials(data));
	rm_this = [];
	for i = haz_data
		if length(spikes) < i
			rm_this = [rm_this i];
		else
			if length(spikes(i).A) < 2e5 % 10 seconds
				rm_this = [rm_this i];
			end

			try
				if isempty(setdiff(1:width(spikes(i).A),find(spikes(i).discard)))
					rm_this = [rm_this i];
				end
			catch
			end

			if use_light
				if strfind(ControlParadigm(i).Name,'OdourFlicker')
					rm_this = [rm_this i];
				end
			else
				if strfind(ControlParadigm(i).Name,'LightFlicker') 
					rm_this = [rm_this i];
				end
				if strfind(ControlParadigm(i).Name,'LongPulse')
					rm_this = [rm_this i];
				end

				if strfind(ControlParadigm(i).Name,'Microscope')
					rm_this = [rm_this i];
				end

			end

		end
	end
	haz_data = setdiff(haz_data,rm_this);
end

for i = 1:length(haz_data)
	time = 1e-4*(1:length(data(haz_data(i)).PID));
	try
		[temp, tA] = spiketimes2f(spikes(haz_data(i)).A,time,1e-3);
	catch
		keyboard
	end
	rm_this = sum(temp) == 0;

	% add PID
	if ~use_light
		temp_stim = temp;
		for j = 1:width(temp)
			temp_stim(:,j) = interp1(time,data(haz_data(i)).PID(j,:),tA);
		end
	else
		temp_stim = interp1(time,ControlParadigm(haz_data(i)).Outputs(1,:),tA)';
		temp_stim = repmat(temp_stim,1,width(temp));
	end

	temp_stim(:,rm_this) = [];
	temp(:,rm_this) = [];

	% censor parts of the trace, if need be. 
	if isfield(spikes,'use_trace_fragment')
		if ~isempty(spikes(haz_data(i)).use_trace_fragment)
			if width(spikes(haz_data(i)).use_trace_fragment) < width(temp)
				for ww = 1:width(spikes(haz_data(i)).use_trace_fragment)
					censor_this = ~spikes(haz_data(i)).use_trace_fragment(ww,:);
					censor_this = censor_this(:,1:10:end);
					temp(censor_this,ww) = NaN;
				end
			else
				if any(rm_this)
					try
						censor_this = ~spikes(haz_data(i)).use_trace_fragment(~rm_this,:);
					catch
						disp('93')
						keyboard
					end
					censor_this = censor_this(:,1:10:end);
				else
					censor_this = ~spikes(haz_data(i)).use_trace_fragment;
					censor_this = censor_this(:,1:10:end);
				end
				if length(censor_this) ~= length(temp) || width(censor_this) ~= width(temp)
					censor_this = censor_this(:);
					try
						temp(censor_this) = NaN;
					catch
						disp('78')
						keyboard
					end
				else
					temp(censor_this) = NaN;
				end
			end
		end
	end
	
	% add to fA
	try
		resp = [resp temp];
	catch
		disp(datapath)
		disp('56')
		keyboard
	end
	paradigm = [paradigm (i + p_offset)*(ones(1,width(temp)))];

	stim = [stim temp_stim];

	ParadigmNames = [ParadigmNames strrep(ControlParadigm(haz_data(i)).Name,'_','-')];
end

% cache data
cached_data.stim = stim;
cached_data.resp = resp;
cached_data.paradigm = paradigm;
cached_data.ParadigmNames = ParadigmNames;

cache(h,cached_data);
