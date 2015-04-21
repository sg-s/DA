% MechanismAnalysis_PrepData
% this function prepares data from the experiments investigating the mechanism of adaptation
% for a use example, look at Mechanism3.m
% 
% created by Srinivas Gorur-Shandilya at 10:22 , 21 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,stim,resp,ParadigmNames,paradigm,use_light)

load(datapath);
if isempty(paradigm)
	p_offset = 0;
else
	p_offset = max(paradigm);
end

for i = 1:length(haz_data)
	time = 1e-4*(1:length(data(haz_data(i)).PID));
	[temp, tA] = spiketimes2f(spikes(haz_data(i)).A,time,1e-3);
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
	
	% add to fA
	resp = [resp temp];
	paradigm = [paradigm (i + p_offset)*(ones(1,width(temp)))];

	stim = [stim temp_stim];

	ParadigmNames = [ParadigmNames strrep(ControlParadigm(haz_data(i)).Name,'_','-')];
end