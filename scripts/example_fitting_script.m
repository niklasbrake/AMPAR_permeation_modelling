% Parameter format
%    pars 1-7: A1_A2 (k1x,k1rx,k2x,k2rx,k3x,k4rx,k4x),
%    pars 8-14: gamma2 (k1x,k1rx,k2x,k2rx,k3x,k4rx,k4x),
%    pars 15-21: gamma2+CNIH3 (k1x,k1rx,k2x,k2rx,k3x,k4rx,k4x),
%    pars 22-28: R607E_gamma2 (k1x,k1rx,k2x,k2rx,k3x,k4rx,k4x),
%    par 29: delta 1 (A1_A2, gamma2, gamma2+CNIH3, R607E_gamma2)
%    par 30: delta 3 (A1_A2, gamma2, gamma2+CNIH3)
%    par 31: delta 3 (R607E_gamma2)

% Load data with individual patches
addpath('../models')
addpath('../models/model_fits')

load('data_all_patches.mat')

% Use best fit as starting point in parameter space
load('model4.2_fit.mat','P');
init_pars = zeros(1,7*4+2);
fn = fieldnames(P);
for i = 1:4
    temp = P.(fn{i});
    init_pars(7*(i-1)+1:7*i) = temp(1:7);
end
init_pars(29:31) = [0.1,0.5,0.5];

% Runn 100 boostrap fits of the data
for i = 1:100
    fprintf('Bootstrap fitting %d/100...',i);
    temp = bootstrap(data,init_pars);
    fn = fieldnames(temp);
    for j = 1:length(fn)
        P_bootstrap.(fn{j})(i,:) = temp.(fn{j});
    end
    fprintf('complete.\n',i);
end

function P = bootstrap(data,init_pars)

    % Bootstrap sample of data
    fn = fieldnames(data);
    for i = 1:4
        n = length(data.(fn{i}).Ca);
        idcs = randi(n,n,1);
        data.(fn{i}).Ca = data.(fn{i}).Ca(idcs);
        data.(fn{i}).gnorm = data.(fn{i}).gnorm(idcs,:);
    end

    % Options for fitting algorithm
    myfun = @(optimValues,options) options.InitialTemperature.*0.96.^optimValues.k;
    options = optimoptions(@simulannealbnd,'PlotFcn','saplotbestf','MaxStallIterations',300,'ReannealInterval',70,'InitialTemperature',.3,'FunctionTolerance',1e-6,'TemperatureFcn',myfun,'PlotFcns',{},'Display','off');
    F = @(pars) problem_fun(data,pars);
    ub = ones(1,7*4+2);
    lb = zeros(1,7*4+2);

    % Run optimization algorithm on model
    pars = simulannealbnd(F,init_pars,lb,ub,options);

    % Format parameters into structure for each AMPAR
    for i = 1:4
        field = fn{i};
        P.(field) = pars([7*(i-1)+1:7*i,29,30]);
    end
    P.(field) = pars([7*(i-1)+1:7*i,29,31]);

    % model.fun = @model4;
    % model.P = P;
    % figureNB;
    % for i = 1:length(fn)
    %     plot_results(data,model,fn{i},i);
    % end
end

function F = problem_fun(data,pars)

    % Assign parameters to each AMPAR
    pars0 = zeros(7,4);
    pars0([8,9],:) = repmat(pars(29:30)',[1,4]);
    pars0(9,4) = pars(31);
    for i = 1:4
        pars0(1:7,i) = pars(7*(i-1)+1:7*i);
    end

    % Loop through each AMPAR and track cummulative error
    F = 0;
    fn = fieldnames(data);
    for i = 1:length(fn)
        % Split data based on calcium concentration
        [I,Ca] = findgroups(data.(fn{i}).Ca);
        % Simulate model for each calcium concentration
        G = model4(pars0(:,i),Ca,data.(fn{i}).V*1e-3);
        % Compute L2 error
        err = (data.(fn{i}).gnorm-G(I,:)).^2;
        % Take mean to weight each AMPAR the same (scaling factor
        % arbitrary but related to optimiation algorithm options)
        F = F + 100*nanmean(err(:));
    end
end
