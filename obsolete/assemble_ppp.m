% Assemble ppp data
% small script to assemble ppp data from raw into something usable
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

allfiles = dir('2014_10_09_CSF*.mat');
c = 1;
add_this=1;

for i = 1:length(allfiles)
	textbar(i,length(allfiles))
	load(allfiles(i).name)

	for j = 1:length(spikes)
		add_this = 1;
		if length(spikes(j).A) > 1 && isvector(spikes(j).A) 
		else
			add_this = 0;
		end
		if isfield(spikes,'discard')
			if spikes(j).discard == 1
				add_this = 0;
			end
		end 
		if add_this
			% has data

			% make a time vector
			time = (1:length(spikes(j).A))/SamplingRate;

			% build a firing rate estimate
			[f,t] = spiketimes2f(spikes(j).A,time,3e-3,3e-2);

			ppp_data(c).time = t;
			ppp_data(c).f = f;

			% load the PID
			ppp_data(c).PID=interp1(time,data(j).PID,t);


			us = strfind(ControlParadigm(j).Name,'_');

			% check if it makes sense 
			if max(ControlParadigm(j).Outputs(6,:))
				% there is a probe

				% load the valve
				ppp_data(c).c_valve=interp1(time,ControlParadigm(j).Outputs(5,:),t);
				ppp_data(c).p_valve=interp1(time,ControlParadigm(j).Outputs(6,:),t);


				% pull out the conditioning pulse height
				a = strfind(ControlParadigm(j).Name,'c_');
				z = us(find(us>a+1,1,'first'));
				ppp_data(c).c_height = str2double(ControlParadigm(j).Name(a+2:z-1));

				% pull out the probe pulse height
				a = strfind(ControlParadigm(j).Name,'p_');
				z = us(find(us>a+1,1,'first'));
				ppp_data(c).p_height = str2double(ControlParadigm(j).Name(a+2:z-1));


				
			else

				% there is no probe, which is silly. swap the probe and conditioning valves around
				ppp_data(c).c_valve=interp1(time,ControlParadigm(j).Outputs(6,:),t);
				ppp_data(c).p_valve=interp1(time,ControlParadigm(j).Outputs(5,:),t);

				% pull out the conditioning pulse height
				a = strfind(ControlParadigm(j).Name,'p_');
				z = us(find(us>a+1,1,'first'));
				ppp_data(c).c_height = str2double(ControlParadigm(j).Name(a+2:z-1));

				% pull out the probe pulse height
				a = strfind(ControlParadigm(j).Name,'c_');
				z = us(find(us>a+1,1,'first'));
				ppp_data(c).p_height = str2double(ControlParadigm(j).Name(a+2:z-1));


			end

			% add some metadata
			ppp_data(c).original_name = allfiles(i).name;

			ppp_data(c).paradigm_name = ControlParadigm(j).Name;


			% pull out the lag
			a = strfind(ControlParadigm(j).Name,'l_');
			z = us(find(us>a+1,1,'first'));
			ppp_data(c).lag = str2double(ControlParadigm(j).Name(a+2:z-1));

			% pull out the width
			a = strfind(ControlParadigm(j).Name,'w_');
			z = us(find(us>a+1,1,'first'));
			if isempty(z)
				z = length(ControlParadigm(j).Name);
			end
			ppp_data(c).width = str2double(ControlParadigm(j).Name(a+2:z));

			% neuron info
			if any(strfind(allfiles(i).name,'ab2'))
				ppp_data(c).neuron = 'ab2';
			else
				ppp_data(c).neuron = 'ab3';
			end

			
			c = c+ 1;

		end
	end
	clear j
end
clear i
