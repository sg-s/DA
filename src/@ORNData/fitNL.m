% ORNData/fitNL.m
% fits a nonlinearity to the data
% 
% created by Srinivas Gorur-Shandilya at 1:34 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [o,ff] = fitNL(o)


% recursively call this to work on arrays of objects
if length(o) > 1
	original_dimensions = size(o);
	o = o(:);
	for i = 1:length(o)
		o(i) = fitNL(o(i));
	end
	o = reshape(o,original_dimensions);
	return
end

% respect constrains on where to calculate this 
if isempty(o.use_these_trials)
	utt = 1:o.n_trials;
else
	utt = o.use_these_trials;		
end
if isempty(o.use_this_segment)
	uts = 1:length(o.stimulus);
else
	uts = o.use_this_segment;
end	




% check if data is there
if ~isempty(o.firing_rate) && ~isempty(o.firing_projected) && ~isempty(o.stimulus)
	r = o.firing_rate(uts,:);
	p = o.firing_projected(uts,:);

	for i = 1:width(r)
		if ismember(i,find(utt))
			temp1 = p(:,i); temp2 = r(:,i);
			rm_this = isnan(temp1) | isnan(temp2);
			temp1(rm_this) = []; temp2(rm_this) = [];
			temp2 = temp2/max(temp2);

			ft = fittype('hill2(x,k,n,x_offset)');
			ff = fit(temp1,temp2,ft,'StartPoint',[.5 2 nanmean(temp1)],'Lower',[0 1 -1],'Upper',[max(temp1)*100 10 1],'MaxIter',1e4);

			% figure, hold on
			% plot(temp1,temp2,'k.')
			% plot(temp1,ff(temp1),'rx')

			temp1 = o.firing_projected(:,i);
			% temp1 = temp1 - min_offset; temp1 = temp1/max_offset;
			o.firing_projected(:,i) = ff(temp1)*max(o.firing_rate(uts,i));


		end
	end

end
if ~isempty(o.LFP) && ~isempty(o.LFP_projected) && ~isempty(o.stimulus)
	r = o.LFP(uts,utt);
	p = o.LFP_projected(uts,utt);

	for i = 1:width(r)
		if ismember(i,find(utt))
			temp1 = p(:); temp2 = r(:);
			rm_this = isnan(temp1) | isnan(temp2);
			temp1(rm_this) = []; temp2(rm_this) = [];

			temp1 = -temp1; temp2 = - temp2;
			temp1 = temp1 - min(temp1);
			temp2 = temp2 - min(temp2);
			temp1 = temp1/max(temp1);
			temp2 = temp2/max(temp2);
			ft = fittype('hill2(x,k,n)');
			ff = fit(temp1,temp2,ft,'StartPoint',[.5 2],'Lower',[0 1],'Upper',[1 4],'MaxIter',1e3);

			% ff = fit(temp1(1:100:end),temp2(1:100:end),'rat24','StartPoint',[-2 2 1 -1 -2 .1 1]);

			temp1 = -o.LFP_projected(:,i); 
			temp1 = temp1 - min(temp1); temp1 = temp1/max(temp1);
			temp1 = -ff(temp1);

			ff2 = fit(temp1(~isnan(temp1)),o.LFP(~isnan(temp1),i),'poly1');
			temp1 = ff2(temp1);


			o.LFP_projected(:,i) = temp1;


		end
	end
end