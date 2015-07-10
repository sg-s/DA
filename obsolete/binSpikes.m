% binSpikes.m
% bins spikes from a raw spike sequence into 1ms bins
% 
% created by Srinivas Gorur-Shandilya at 4:20 , 25 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function bs = binSpikes(spikes)

if size(spikes,2) > size(spikes,1)
	spikes = spikes';
end

bs = zeros(floor(length(spikes)/10),width(spikes));


for j = 1:width(spikes)
	for i = 1:length(bs)-1
		a = 10*(i-1) + 1;
		z = a + 10;
		if any(spikes(a:z,j))
			bs(i,j) = 1;
		end
	end
end