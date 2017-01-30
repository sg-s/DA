% helper function called by script
% this exists solely to clean up other scripts

function [cdata, data] =  assembleScaledNatStim(cdata)


% delete junk 
rm_this =  isnan(sum(cdata.PID)) | isnan(sum(cdata.LFP));
fn = fieldnames(cdata);
for i = 1:length(fn)
	temp = cdata.(fn{i});
	if size(temp,1) == length(rm_this)
		cdata.(fn{i})(rm_this,:) = [];
	elseif size(temp,2) == length(rm_this)
		cdata.(fn{i})(:,rm_this) = [];
	end
end

% only retain data where we have many paradigms on one neuron
useless_neurons = [];
for i = 1:max(cdata.orn)
	if length(unique(cdata.paradigm(cdata.orn== i))) > 1
	else
		useless_neurons = [useless_neurons; i];
	end
end
rm_this = ismember(cdata.orn,useless_neurons);
fn = fieldnames(cdata);
for i = 1:length(fn)
	temp = cdata.(fn{i});
	if size(temp,1) == length(rm_this)
		cdata.(fn{i})(rm_this,:) = [];
	elseif size(temp,2) == length(rm_this)
		cdata.(fn{i})(:,rm_this) = [];
	end
end

clear data
all_orns = unique(cdata.orn);
for i = 1:length(all_orns)

	data(i).R = [];
	data(i).S = [];
	data(i).X = [];
	
	this_orn = all_orns(i);
	

	these_paradigms = unique(cdata.paradigm(cdata.orn == this_orn));
	c = lines(length(these_paradigms));
	time = 1e-3*(1:length(cdata.PID));
	for j = 1:length(these_paradigms)
		plot_this = cdata.paradigm == these_paradigms(j) & cdata.orn == this_orn;
		S = cdata.PID(:,plot_this); 
		X = 10*cdata.LFP(:,plot_this); 
		R = cdata.fA(:,plot_this);

		% remove baseline from each trial
		for k = 1:size(X,2)
			X(:,k) = X(:,k) - mean(X(1:5e3,k));
			S(:,k) = S(:,k) - min(S(1:5e3,k));
		end

		% average
		S = mean(S,2); X = mean(X,2);
		R(:,sum(R)==0) = [];
		R = mean(R,2);

		data(i).R = [data(i).R R];
		data(i).X = [data(i).X X];
		data(i).S = [data(i).S S];


	end

end