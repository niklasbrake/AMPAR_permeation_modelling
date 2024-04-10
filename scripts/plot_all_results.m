addpath('../models')
addpath('../models/model_fits')

data = load('data.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot all model fits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model1 = load('model1_fit.mat','P');
model1.fun = @model1;
model1.name = 'Model 1';
model1.BIC = plot_all(data,model1);

model2 = load('model2_fit.mat','P');
model2.fun = @model2;
model2.name = 'Model 2';
model2.BIC = plot_all(data,model2);

model3 = load('model3_fit.mat','P');
model3.fun = @model3;
model3.name = 'Model 3';
model3.BIC = plot_all(data,model3);

model4 = load('model4_fit.mat','P');
model4.fun = @model4;
model4.name = 'Model 4';
[model4.BIC,model4.SSE] = plot_all(data,model4);

model5 = load('model5_fit.mat','P');
model5.fun = @model5;
model5.name = 'Model 5';
model5.BIC = plot_all(data,model5);

model4_1 = load('model4.1_fit.mat','P');
model4_1.fun = @model4;
model4_1.name = 'Model 4.1';
[~,model4_1.SSE] = plot_all(data,model4_1);

model4_2 = load('model4.2_fit.mat','P');
model4_2.fun = @model4;
model4_2.name = 'Model 4.2';
[~,model4_2.SSE] = plot_all(data,model4_2);

model4_3 = load('model4.3_fit.mat','P');
model4_3.fun = @model4;
model4_3.name = 'Model 4.3';
[~,model4_3.SSE] = plot_all(data,model4_3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot BIC result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian information criterion (BIC) for each model
BIC = [model1.BIC;model2.BIC;model3.BIC;model4.BIC;model5.BIC];
clrs = zeros(5,3);
clrs(1,:) = hsv2rgb([340/360,   50/100,     60/100]);
clrs(2,:) = hsv2rgb([32/360,    50/100,     60/100]);
clrs(3,:) = hsv2rgb([22/360,    25/100,     70/100]);
clrs(4,:) = hsv2rgb([217/360,   50/100,     60/100]);
clrs(5,:) = hsv2rgb([207/360,   25/100,     70/100]);
figureNB;
    axes('ColorOrder',clrs);
    hold on;
    bar(BIC');
    xticks(1:4)
    xticklabels({'A1+A2','A1+A2\gamma2','A1+A2\gamma2+CNIH3','A1+A2_{R607E}\gamma2'});
    ylabel('BIC');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Nested model result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum of squared error (SSE) for nested models 4.0 to 4.3
SSE = [model4.SSE,model4_1.SSE,model4_2.SSE,model4_3.SSE];

% Number of parameters in each nested model
p = [36,33,31,30];

% Calculate total number of data points used for fitting
fn = fieldnames(data);
for i = 1:4
    gnorm(:,:,i) = data.(fn{i}).gnorm;
end
N = sum(~isnan(gnorm(:)));

% Compute F-statistic and p values
for i = 1:4
    for j = i+1:4
        F(i,j) = (SSE(j)-SSE(i))/SSE(i)*(N-p(i))/(p(i)-p(j));
        pvalue(i,j) = 1-fcdf(F(i,j),p(i)-p(j),N-p(i));
    end
end

figureNB;
subplot(2,2,1)
    bar(p,'k')
    ylim([25,40]);
    xticks(1:4)
    xticklabels({'4.0','4.1','4.2','4.3'});
    xlabel('Model');
    ylabel('Parameter count')
subplot(2,2,2)
    bar(SSE,'k')
    xticks(1:4)
    xticklabels({'4.0','4.1','4.2','4.3'});
    xlabel('Model');
    ylabel('Total error (SSE)')
subplot(2,1,2);
    imagesc(fliplr(log10(F(:,2:end))));
    xticks(1:3);
    xticklabels({'4.3','4.2','4.1'});
    yticks(1:3);
    yticklabels({'4.0','4.1','4.2','4.3'});
    C = colorbar;
    C.Label.String = 'log_{10} F-stat';
    colormap(clrsSequential(1e3))
    set(gca,'CLim',[0,2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameters of Model 4.2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot bootstrap uncertainty
model4_2_boot.fun = @model4;
model4_2_boot.P = load('model4.2_bootstrap.mat');
plot_bootstrap_uncertainty(model4_2_boot);

% Extract parameters from boostrap structure
p = [];
ids = [];
fn = fieldnames(model4_2_boot.P);
for j = 1:length(fn)
    p = [p;model4_2_boot.P.(fn{j})];
    ids = [ids;zeros(size(model4_2_boot.P.(fn{j}),1),1)+j];
end
k1x = p(:,1);
k1rx = p(:,2);
k2x = p(:,3);
k2rx = p(:,4);
k3x = p(:,5);
k3rx = k2rx.*k1rx.*k3x./(k2x.*k1x);
d1 = p(:,8);
d3 = p(:,9);

% Compute dissociation constants
p0(:,1) = p(:,2)./p(:,1);
p0(:,2) = p(:,4)./p(:,3);
p0(:,3) = k3rx./k3x;
p0(:,4) = p(:,7)./p(:,6);
p0(:,5) = p(:,8);
p0(:,6) = p(:,9);
p1 = splitapply(@mean,p0,ids);
p2 = splitapply(@std,p0,ids);

% Plot histogram of delta values
figureNB(7,4);
    histogram([d1;d3],'BinWidth',0.005,'FaceColor','k','FaceAlpha',0.6,'EdgeColor','none','Normalization','pdf');
    xlim([0,1]);
    yticks([]);
    xlabel('Fraction of mem. potential');
    ylabel('Probability');
    gcaformat;

% Plot dissociation constants associated with each transition in model
AMPARs = {'A1+A2','A1+A2\gamma2','A1+A2\gamma2+CNIH3','A1+A2_{R607E}\gamma2'};
dh = [-0.15,-0.05,0.05,0.15]*2;
clrs = lines(4)*0.8;
figureNB(11,4);
axes('Position',[0.10, 0.20, 0.45, 0.75])
for i = 1:4
    for j = 1:4
        bar(i+dh(j),p1(j,i),'EdgeColor','none','BarWidth',0.085*2,'FaceColor',clrs(j,:));
        line(i+dh(j)*[1,1],p1(j,i)+p2(j,i).*[-1,1],'color','k','LineWidth',1);
        hold on;
    end
    text(0.6,2.3-0.2*i,AMPARs{i},'color',clrs(i,:),'FontSize',6)
end
xlim([0.5,4.5])
xticks(1:4)
ylim([-0.1,2.2])
xtl = {'Ca_O\leftrightarrowB_1','B_1\leftrightarrowB_2', ...
    'Ca_O\leftrightarrowB_2', 'Ca_i\leftrightarrowB_2'};
xticklabels(xtl)
xlabel('Reaction')
ylabel('Dissociation constant')
gcaformat;


clrs = lines(4);
clrs(2,:) = clrs(3,:);
clrs(3,:) = clrs(4,:);
clrs(4,:) = clrs(1,:);
clrs(1,:) = clrs(2,:);
mrks = {'s','^','<','o'};
Q = [0.06,0.42,0.96,4.65];
axes('Position',[0.66, 0.20, 0.24, 0.75])
for i = [4,1,3]
    E = plot(Q,p1(:,i),'.-','color',clrs(i,:),'MarkerSize',4,'Marker',mrks{i},'LineWidth',1);
    E.MarkerEdgeColor = 'none';
    E.MarkerFaceColor = clrs(i,:);
    for j = 1:4
        y = p0(ids==j,i);
        line(Q(j)*[1,1],mean(y)+std(y)*[-1,1],'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    ylim([-0.1,2.2])
    xlim([1/32,8])

    xlabel('P_{Ca}/P_{Na}')
    ylabel('Ca^{2+} dissociation (mM)')
    text(8,p1(end,i),xtl{i},'FontSize',6,'color',clrs(i,:));

    [rho ,p] = corr(Q(ids)',p0(:,i),'type','spearman');
    fprintf('%s: \\rho=%f (p=%d)\n',xtl{i},rho, p)
end
set(gca,'xscale','log')
xticks([0.25/4,0.25,1,4])
xticklabels({'1/8','1/4','1','4'})
gcaformat

function [BIC,SSE] = plot_all(data,model)
% PLOT_ALL  loops through each AMPAR and plots the fit of model. Wrapper
%   for PLOT_single.
%   [BIC,SSE] = plot_all(data,model) takes data and model structures and
%       returns the Bayesian information criterion (BIC) and sum of
%       squared error (SSE).

    fn = fieldnames(data);
    fig = figureNB(23,11);
    fig.Name = model.name;
    for i = 1:length(fn)
        [BIC(i),SSE(i)] = plot_single(data,model,fn{i},i);
    end
    SSE = sum(SSE(:));
end


function [BIC,SSE] = plot_single(data,model,field,i)
% PLOT_SINGLE  Plots plots of fit of model to AMPAR data and returns
%   goodness of fit.
%   [BIC,SSE] = plot_single(data,model,field,i) takes data and model structures,
%   as well as the name of the AMPAR (field) and its ordinal number (i). Returns
%   the Bayesian information criterion (BIC) and sum of squared error (SSE).

    % Load data from data structure and simulate model
    data = data.(field);
    V = data.V;
    gnorm = data.gnorm;
    gnorm_SE = data.gnorm_SE;
    model_gnorm = model.fun(model.P.(field),data.Ca,1e-3*data.V);

    % Compute SSE
    err = (gnorm-model_gnorm).^2;
    SSE = nansum(err(:));

    % Compute BIC
    n = sum(~isnan(err(:)));
    k = numel(model.P.(field));
    BIC = n*log(SSE) + k*log(n);

    % PLot model output with higher resolution in calcium space.
    Ca = 10.^linspace(-3,3,1e3);
    model_gnorm = model.fun(model.P.(field),Ca,1e-3*data.V);

    for j = 1:length(V)
        subplot(4,8,8*(i-1)+j);
        for iC = 1:length(data.Ca)
            plot(data.Ca(iC),gnorm(iC,j),'.k','MarkerSize',10);
            hold on;
            plot(data.Ca(iC)*[1,1],gnorm(iC,j)+gnorm_SE(iC,j)*[-2.33,2.33],'color','k','LineWidth',1);
        end

        plot(Ca,model_gnorm(:,j),'color','r','LineWidth',1.5);
        set(gca,'xscale','log');
        xlim([0.02,200]);
        drawnow;
        if(i==1)
            title(sprintf('V = %d mV',V(j)));
        end
        if(j==1)
            ylabel(strrep(field,'_','\_'));
        end
        xlabel('[Ca^{2+}] (mM)')
        xticks([0.1,1,10,100]);
        xticklabels([0.1,1,10,100]);
        gcaformat;
        if(i==1 & j==1)
            text(0.03,0.2,model.name,'color','r','FontSize',7)
            text(0.03,0.35,'Data','color','k','FontSize',7)
        end
    end
end
function plot_bootstrap_uncertainty(model)
    load('data_all_patches.mat')
    fn = fieldnames(data);
    xtl = {'A1+A2','A1+A2\gamma2','A1+A2\gamma2+CNIH3','A1+A2_{R607E}\gamma2'};
    figureNB(18.3,8);
    for i = 1:4
        field = fn{i};
        data0 = data.(field);
        clrs = clrsSequential(12); clrs = clrs(5:end,:);
        V = data0.V;
        gnorm = data0.gnorm;

        Ca = 10.^linspace(-3,3,100);

        p = model.P.(field);
        for k = 1:100
            model_gnorm(:,:,k) = model.fun(p(k,:),Ca,1e-3*data0.V);
        end
        model.gnorm.(field) = model_gnorm;

        clrs = lines(4)*0;
        count = 0;
        for j = 1:8
            count = count+1;
            subplot(4,8,8*(i-1)+count)
            plot(data0.Ca,gnorm(:,j),'.','MarkerSize',5,'color','k');
            hold on;
            plot(Ca,mean(model_gnorm(:,j,:),3),'color','r','LineWidth',1);
            y0 = min(model_gnorm(:,j,:),[],3);
            y1 = max(model_gnorm(:,j,:),[],3);
            fill([Ca,fliplr(Ca)],[y0',fliplr(y1')],'r','FaceAlpha',0.2,'EdgeColor','none');
            set(gca,'xscale','log');
            xlim([0.02,200]);
            drawnow;
            xticks([0.1,1,10,100]);
            yticks([0,0.5,1]);
            xticklabels([0.1,1,10,100]);
            gcaformat;
            if(i==1)
                title(sprintf('V = %d mV',V(j)));
            end
            if(i==4)
                xlabel('[Ca^{2+}] (mM)')
            end
        end
        subplot(4,8,8*(i-1)+1)
        ylabel(xtl{i});
        drawnow;
    end
end
function varargout = figureNB(x,y)
    if(nargin==0)
        x=9;
        y=9;
    end
    set(0,'units','centimeters');
    SS = get(0,'screensize');

    fig = figure('color','w','units','centimeters','ToolBar','none');
    fig.Position = [SS(3)/2-x/2,SS(4)/2-y/2,x,y];
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[x,y],'Renderer','Painters');
    set(fig, 'DefaultAxesFontName', 'Arial');
    set(fig, 'DefaultTextFontName', 'Arial');

    if(nargout==1)
        varargout{1} = fig;
    else
        varargout = {};
    end
end
function gcaformat
    ax = gca;
    set(ax,'box','off');
    set(ax,'tickdir','out');
    set(ax,'linewidth',0.75);
    set(ax,'fontsize',7);
    set(ax.Title,'FontSize',7);
    xax = get(ax,'xaxis');
    xax.Label.FontSize = 7;
    xax.TickLabelRotation = 0;
    yax = get(ax,'yaxis');
    for i = 1:length(yax)
        yax(i).Label.FontSize = 7;
        yax(i).TickLabelRotation = 0;
    end
    tickLength = 0.05; % 1/2 mm
    U = ax.Units;
    set(ax,'Units','centimeters');
    L = max(ax.Position(3:4));
    t = tickLength/L;
    set(ax,'TickLength',[t,t]);
    set(ax,'Units',U);
end
function clrs = clrsSequential(N)
    sequential_CM = [255,255,255; ...
                    255,247,188; ...
                    254,227,145; ...
                    254,196,79; ...
                    251,154,41; ...
                    236,112,20; ...
                    204,76,2; ...
                    153,52,3; ...
                    102,37,6]/255;
    I1 = linspace(0,1,9);
    I2 = linspace(0,1,N);
    clrs = interp1(I1,sequential_CM,I2);
end