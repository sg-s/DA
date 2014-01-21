function [model_values] = build_LNmodel_optreg(time, stim, f, fsd, trials, memory, shift_input, deltat, reg_vec, fitfun, figon, regL2)

    disp('build_LNmodel_optreg starting...')
    if shift_input>0        
        time = time(1:end-shift_input+1);
        stim = stim(shift_input:end);
        f = f(1:end-shift_input+1);
        fsd = fsd(1:end-shift_input+1);
    end

    n = length(reg_vec);
    NR = zeros(1,n);
    lambda = zeros(1,n);
    a = zeros(n,memory+1);
    D = [];
    
    wb = waitbar(0, 'regularizing' );
    for i=1:n       
        waitbar(i/n, wb, ['regularizing ' num2str(reg_vec(i))]);

        [NR(i), lambda(i), a(i,:), D(i,:),~,a_notnormalised] = build_LNmodel_reg(time, stim, f, fsd, trials, memory, deltat, reg_vec(i), fitfun, 0, regL2);
    end
    close(wb)
    [temp ireg_opt] = min(NR);
    
    NR_final = NR(ireg_opt);
    reg_final = reg_vec(ireg_opt);
    lambda_final = lambda(ireg_opt);
    a_final = a(ireg_opt, :)';
    D_final = D(ireg_opt, :);
    
    if n>1
        figure; hold on
        plot(log10(reg_vec), NR, 'ko-');
        plot(log10(reg_final), NR_final, '>');
    end

    time2 = time(memory+1:end);
    stim2 = stim(memory+1:end);
    f2 = f(memory+1:end);
    fsd2 = fsd(memory+1:end);

    [fSTA, fLN,model_values] = build_LNmodel_allvalues(time2, stim, f2, fsd2, trials, memory, deltat, fitfun, figon, a_final, D_final);
    % add the not-normalised filter back
    model_values.a_notnormalised = a_notnormalised;
    disp('build_LNmodel_optreg DONE')
end

