% KineticsAnalysis.m
% similar to Gain Analysis, performs a kinetics analysis on flickering stimulus data
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = KineticsAnalysis(data,history_lengths,pred,plothere)

switch nargin
case 0
	help KineticsAnalysis
	return
case 1
	history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
	example_history_length = 0.135;
case 2
	if isempty(history_lengths)
		history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
		example_history_length = 0.135;
	end
case 3
	if isempty(history_lengths)
		history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
		example_history_length = 0.135;
	end
end

% unpack data
stim = data.PID(:);
resp = data.ORN(:);
time = data.time(:);
valve = data.Valve(:);
dt = mean(diff(time));
hl = round(history_lengths/dt);


% compute shat(the smoothed stimulus)
shat = ComputeSmoothedStimulus(stim,hl);

% fix the delay b/w the valve and the orn
d = 24;
time(1:d) = [];
resp(1:d) = [];
valve(end-d+1:end) = [];
pred(1:d) = [];


% find the time to half max
[ValveOns,ValveOffs]=ComputeOnsOffs(valve);
t_half_max = NaN*ValveOns;
t_half_max_loc = NaN*ValveOns;

% debug
figure, hold on
plot(time,resp)



for i = 1:length(ValveOns)-1
	% determine how far to look in the future
	a = ValveOns(i)-5;
	z = ValveOffs(i)+5;
	this_resp = resp(a:z);
	this_resp(1:min_loc) = 0;

	[min_val,min_loc]=min(resp(a:a+10));
	[max_val,max_loc]=max(this_resp);

	
	scatter(time(a+min_loc),min_val,'k','filled')
	scatter(time(a+max_loc),max_val,'r','filled')

	% find when resp exceeds 1/2 this first
	
	t_half_max(i) = find(this_resp>min_val+(max_val-min_val)/2,1,'first');
	t_half_max_loc(i) = a+t_half_max(i)-1;

	scatter(time(t_half_max_loc(i)),resp(t_half_max_loc(i)),'b','filled')

	

end



t_half_max_loc(end) = [];
t_half_max(end) = [];


% calculate gain
gain = resp./pred;
