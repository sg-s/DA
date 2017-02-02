function plotScaledNatStimFilters(data,filter_type)

w = min([1300 500*length(data)]);
figure('outerposition',[0 0 w 901],'PaperUnits','points','PaperSize',[1300 901]); hold on
time = 1e-3*(1:length(data(1).S));

ss = 2;
history_length = 500;
scatter_size = 20;



for i = 1:length(data)
	c = lines(size(data(i).S,2));

	% compute overall sclae
	shat = computeSmoothedStimulus(data(i).S(:),history_length);
	MS = max(shat);
	mS = min(shat);

	for j = 1:size(data(i).S,2)
		S = data(i).S(:,j);
		X = data(i).X(:,j);
		R = data(i).R(:,j);
		switch filter_type
		case 'LFP' 
			K1 = fitFilter2Data(S,X,'offset',200);
		otherwise
			if max(R) > 0 
				K1 = fitFilter2Data(S,R,'offset',200);
			else
				continue
			end
		end
		K1 = K1(100:end-100);
		filtertime = 1e-3*(1:length(K1)) - .1;
		fp = convolve(time,S,K1,filtertime);
		K1 = K1/(nanstd(fp)/nanstd(S)); % normalise correctly 
		fp = convolve(time,S,K1,filtertime);

		subplot(2,length(data),i); hold on
		plot(filtertime,K1,'MarkerSize',20,'Color',c(j,:))

		ylabel('Filter (norm)')
		xlabel('Lag (s)')


		subplot(2,length(data),length(data)+i); hold on
		switch filter_type
		case 'LFP' 
			%plot(fp(1:ss:end),X(1:ss:end),'.-','MarkerSize',20,'Color',c(j,:))
			% colour by mean stim in preceding 500ms 
			shat = computeSmoothedStimulus(S,history_length);
			shat = shat-mS;
			shat = shat/MS;
			shat = 1+floor(shat*99);
			shat(isnan(shat)) = 1;

			% make the output analysis plot
			cc = parula(100);
			c = cc(shat,:);
			scatter(fp(1:ss:end),X(1:ss:end),scatter_size,c(1:ss:end,:),'filled')
			ylabel('LFP  (mV)')
		otherwise
			% colour by mean stim in preceding 500ms 
			shat = computeSmoothedStimulus(S,history_length);
			shat = shat-mS;
			shat = shat/MS;
			shat = 1+floor(shat*99);
			shat(isnan(shat)) = 1;

			% make the output analysis plot
			cc = parula(100);
			c = cc(shat,:);
			scatter(fp(1:ss:end),R(1:ss:end),scatter_size,c(1:ss:end,:),'filled')
			%plot(fp(1:ss:end),R(1:ss:end),'.-','MarkerSize',20,'Color',c(j,:))
			ylabel('Firing rate  (Hz)')
		end

		xlabel('Proj. Stim (V)')
	end
end

switch filter_type
case 'LFP'
	for i = 1:length(data)
		subplot(2,length(data),i); hold on
		set(gca,'YDir','reverse')
		subplot(2,length(data),length(data)+i); hold on
		set(gca,'YDir','reverse','XDir','reverse')
	end
otherwise 
end