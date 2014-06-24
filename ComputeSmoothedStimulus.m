function [shat] = ComputeSmoothedStimulus(stimulus,hl)
shat = NaN(length(hl),length(stimulus));
for i = 1:length(hl)
	if hl(i) == 0
		shat(i,:) = stimulus;
	else
		shat(i,:) = filter(ones(1,hl(i))/hl(i),1,stimulus);
		shat(i,1:hl(i)) = NaN;
	end
 	
end
