% temp script to fit aNLN2 model to nat. stim. data

while 1

	% ab2A 2 butanone

	cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
	[~, data] =  assembleScaledNatStim(cdata);

	for this_orn = 1:length(data)
		clear fd
		for i = 1:size(data(this_orn).X,2)
			S = data(this_orn).S(:,i); S = S - min(S);
			R = data(this_orn).R(:,i);
			fd(i).stimulus = S;
			fd(i).response = R;
			fd(i).response(1:5e3) = NaN;
		end
		if this_orn ~= 2
			fd(1) = [];
		end


		% fit 
		fitModel2Data(@aNLN2,fd,'nsteps',200,'make_plot',false);

	end

	% ab3A 2-butanone 
	clear data fd
	% get all data 
	cdata = consolidateData2(getPath(dataManager,'c2bce18a6b0a7e89e9c6832dcc27e39b'));
	[~, data] =  assembleScaledNatStim(cdata);


	S = data(this_orn).S; S = S - min(S);
	R = data(this_orn).R;
	fd.stimulus = S;
	fd.response = R;
	fd.response(1:5e3) = NaN;

	% fit 
	fitModel2Data(@aNLN2,fd,'nsteps',200,'make_plot',false);

	end




end % end infinite loop