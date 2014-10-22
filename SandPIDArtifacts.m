% SandPIDArtifacts.m
% sands down spikes on the end of every PID pulse, which are believed to be artifacts from the electrical noise of valve switching
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function data = SandPIDArtifacts(data,td)

[~,offs] = ComputeOnsOffs(data(td).Valve);
data(td).PID2 = data(td).PID;

for i = 1:length(offs)-1
	x = offs(i)-10:offs(i)+30;
	y = data(td).PID(x);
	xx = x(22:28);
	x(22:28) = [];
	y(22:28) = [];
	yy = interp1(x,y,xx,'pchip');
	data(td).PID2(offs(i)+12:offs(i)+18) = yy;

end