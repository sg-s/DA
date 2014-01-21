function [NR_final, lambda_final, a_final, D_final,fitN,a_notnormalised] = build_LNmodel_reg(time, stim, f, fsd, trials, memory, deltat, reg, fitfun, figon, regL2)

Ltot = length(stim);% length of stimulus/response vector

L2 = floor(2*Ltot/3);  % 2/3 length of whole data, for building filter
Ltest = floor(Ltot/3);  % 1/3 length, for validation
i2 = (1:L2)'; % indices of data used for filter computation
itest = (Ltot-Ltest+1:Ltot)';  % indicies of data used for check
NR = [];
for jj = 1%:2
    stim2 = stim(i2);    
    f2 = f(i2);
    fsd2 = fsd(i2);
    time2 = time(i2);
        
    f2 = f2(memory+1:end);
    fsd2 = fsd2(memory+1:end);
    time2 = time2(memory+1:end);    
    s = zeros(L2-memory, memory+1);
    for i=1:L2-memory
        s(i,:) = stim2(memory+i:-1:i);
    end
    
    % L2 regularized linear regression
    if regL2 
        C = s'*s + reg*eye(memory+1);    % C = cov(s) + reg*eye(memory+1);
        b = C\(s'*f2);%./sum(f); % b is the bloody filter
        
    else
%         [s_ s_mean s_var] = normalize(s);
%         f2_ = center(f2);
        [b info] = elasticnet(s, f2, reg, 4*10^4, 0);
        
%         [b info] = elasticnet(s, f2, reg);
%         b(:,1)=[]; info.lambda(:,1)=[];
%         b_ind = find(info.lambda<500,1);
%         b = b(:,b_ind);
%         info.lambda = info.lambda(:,b_ind);

%         if reg > find_lambdamax_l1_ls(s', f2)
%             display('ERROR: reg > reg_MAX');
%             quit;
%         else
%             [a, status] = l1_ls(s, f2, reg);
%         end
    end
    
    %%
    Dl = [];
    for lambda_i = 1:size(b,2)
        a = b(:,lambda_i);
        a = a/max(a);
        fSTA = s*a;
        
        % static non linear part
        nbins = 20;
        if fitfun==1; mfSTA = min(fSTA);
        else  mfSTA = 0; end
        fSTA = fSTA - mfSTA;
        ds = (max(fSTA)-min(fSTA))/nbins;
        N=[];
        for i = 1:nbins
            ind1 = fSTA <= min(fSTA)+i*ds;
            ftemp = f2(ind1);
            ind2 = fSTA(ind1) > min(fSTA)+(i-1)*ds;
            N(i) = mean(ftemp(ind2));
        end
        bfSTA = min(fSTA):ds:max(fSTA);
        if fitfun==1 % fit hill function
            fun_fit = fittype('A*(1/(1+(K/x)^n))', 'coefficients',{'A' 'K', 'n'}, 'independent', 'x');
            options = fitoptions('Method', 'NonlinearLeastSquares');
            %         options.StartPoint = [100 abs(mfSTA) 4]; % mofified on 10/19/2011
            options.StartPoint = [100 mean(bfSTA(1:nbins)) 3];
%             fitN = fit(bfSTA(1:nbins)', N', fun_fit, options);
            fitN = fit(fSTA, f2, fun_fit, options);
            D = coeffvalues(fitN);          
            D = [D mfSTA min(f2)];
            fSTA = fSTA+mfSTA;
            fLN = D(1)*(1./(1+(D(2)./fSTA).^D(3)));
        else % fit line/rectifier
            %     th =0;
            %     b_ini = find(bfSTA>th, 1);
            b_ini=1;
            q_fit = fittype('poly1');
%             fitN = fit(bfSTA(b_ini:nbins)', N(b_ini:nbins)', q_fit);
            fitN = fit(fSTA, f2, q_fit);
            D = coeffvalues(fitN);
            fLN = D(1)*fSTA + D(2);
            zeroresp = min(f2);
            fLN = fLN.*(fLN>zeroresp) + zeroresp*(fLN<zeroresp);
        end
%         Pr = mean(abs(f2-fLN).^2);
%         Pn = mean(fsd2.^2);
%         NR(lambda_i) = sqrt(Pr/Pn);
%     
        if figon==3
            figure;
            subplot(2,2,1)
            plot((1:memory+1)*deltat, a, 'r', 'linewidth', 2)
            subplot(2,2,2); hold on
            plot(fSTA, f2, '.', 'color', [.85 .85 .85])
%             plot(bfSTA(1:nbins)+mfSTA, N, 'k', 'linewidth', 2)
            MfSTA = max(fSTA);
            if fitfun==1; ffitted = @(x) D(1)*(1/(1+(D(2)/(x-D(4)))^D(3)));
            else ffitted = @(x) D(1)*x + D(2); end
            fplot(ffitted, [min(fSTA) MfSTA+mfSTA], 'r');
            
            subplot(2,2,3:4); hold on
            errorbar_cazzo(time2, f2, fsd2/sqrt(trials))
            plot(time2, fLN, 'r')
        end
        
    
        % cross-validation
        memory2 = memory;
        
        stimtest = stim(itest);
        ftest= f(itest);
        fsdtest = fsd(itest);
        timetest = time(itest);
        
        ftest = ftest(memory2+1:end);
        fsdtest = fsdtest(memory2+1:end);
        timetest = timetest(memory2+1:end);
        
        stest = zeros(Ltest-memory2, memory2+1);
        for i=1:Ltest-memory2
            stest(i,:) = stimtest(memory2+i:-1:i);
        end
        fSTAtest = stest*a(1:memory2+1);
        fSTAtest = fSTAtest - mfSTA;
        if fitfun==1
            fLNtest = D(1)*(1./(1+(D(2)./fSTAtest).^D(3)));
            fLNtest(fSTAtest<0)=0; % need this to prevent complex numbers coming out from the hill function
        else
            fLNtest = D(1)*fSTAtest + D(2);
            %         fLN = fLN.*(fLN>zeroresp) + zeroresp*(fLN<zeroresp);
        end
        Pr = mean(abs(ftest-fLNtest).^2);
        Pn = mean(fsdtest.^2);
        NR(lambda_i) = sqrt(Pr/Pn);
        Dl(lambda_i, :) = D;
        
        if figon==2
            figure;
            subplot(2,2,1)
            plot((1:memory2+1)*deltat, a(1:memory2+1), 'r', 'linewidth', 2)
            subplot(2,2,2); hold on
            plot(fSTAtest+mfSTA, ftest, '.', 'color', [.85 .85 .85])
            plot(bfSTA(1:nbins)+mfSTA, N, 'k', 'linewidth', 2)
            MfSTA = max(fSTAtest);
            if fitfun==1
                ffitted = @(x) D(1)*(1/(1+(D(2)/(x-D(4)))^D(3)));
                fplot(ffitted, [mfSTA MfSTA+mfSTA], 'r');
            else
                ffitted = @(x) D(1)*x + D(2);
                fplot(ffitted, [min(fSTA) MfSTA+mfSTA], 'r');
            end
            
            subplot(2,2,3:4); hold on
            errorbar_cazzo(timetest, ftest, fsdtest/sqrt(trials))
            plot(timetest, fLNtest, 'r')
            %         title(['NR=' num2str(NR,2) ' CR=' num2str(CR,2)], 'FontSize',24)
        end

    end
    %%
    if regL2        
        id = 1;
        lambda_final = 0;
    else
        [min_res id] = min(NR);
        lambda_final = info.lambda(id);    
    end
    NR_final = NR(id);
    a_final = b(:,id);
    % here carlotta is normalising the filter...
    % comment this out if you dont want to normalise the filter
    a_notnormalised = a_final;
    a_final = a_final/max(a_final);

    D_final = Dl(id,:);
    
    if figon==1
        figure;
        subplot(2,2,1)
        plot((1:memory2+1)*deltat, a_final(1:memory2+1), 'r', 'linewidth', 2)

        subplot(2,2,2); hold on    
        fSTAtest = stest*a_final(1:memory2+1);
        if fitfun==1
            fLNtest = D_final(1)*(1./(1+(D_final(2)./(fSTAtest-D_final(4))).^D_final(3)));
            fLNtest(fSTAtest<D_final(4))=0; % need this to prevent complex numbers coming out from the hill function
        else
            fLNtest = D_final(1)*fSTAtest + D_final(2);
            %         fLN = fLN.*(fLN>zeroresp) + zeroresp*(fLN<zeroresp);
        end
        plot(fSTAtest, ftest, '.', 'color', [.85 .85 .85])
        if fitfun==1
            ffitted = @(x) D_final(1)*(1/(1+(D_final(2)/(x-D_final(4)))^D_final(3)));
        else
            ffitted = @(x) D_final(1)*x + D_final(2);
        end
        fplot(ffitted, [min(fSTAtest) max(fSTAtest)], 'r');
        
        subplot(2,2,3); hold on
        errorbar_cazzo(timetest, ftest, fsdtest/sqrt(trials))
        plot(timetest, fLNtest, 'r')
        
        subplot(2,2,4); hold on
        if length(NR)>1
            plot(log10(info.lambda), NR, '-ok')
            plot(log10(lambda_final), NR_final, 'or')
        end
    end
end



end

