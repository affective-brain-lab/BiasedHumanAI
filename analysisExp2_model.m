clc; clear; close all;
rng('default')

%% Load data
allData = readtable('Exp2.csv');
nSubject = max(allData.subject);
nTrial = 120;
subNum = repmat([1:nSubject]', 1, nTrial)';
allData.subject = subNum(:);

%% Variables
allB0 = [];
allC = [];
allR = [];
allA = [];
allConfi = [];
allColl = [];
allBias = [];
allAcc = [];
allSubject = [];
allCond = [];
nSubject = max(allData.subject);

%% Extract data
for iSubject = 1:nSubject
    
    % Variables
    data = allData(allData.subject == iSubject, :);
    nTrial = height(data);
    
    b0 = mean(data.response(1:30) - data.evidence(1:30));
    allB0 = [allB0; repmat(b0, 90, 1)];
    
    c = data.evidence(31:nTrial);
    allC = [allC; c];
    
    r = data.response(31:nTrial);
    allR = [allR; r];
    
    a = data.responseAI(31:nTrial);
    allA = [allA; a];
    
    confi = data.confidence(31:nTrial);
    confi = (confi - min(confi)) ./ (max(confi) - min(confi));
    confi = max(confi) - confi;
    allConfi = [allConfi; confi];
    
    coll = data.collaboration(31:nTrial);
    allColl = [allColl; coll];
    
    bias = data.response(31:nTrial) - data.evidence(31:nTrial);
    allBias = [allBias; bias];
    
    accuracy = abs(data.response(31:nTrial) - data.evidence(31:nTrial));
    allAcc = [allAcc; accuracy];
    
    allSubject = [allSubject; repmat(iSubject, 90, 1)];
    
    cond = data.condition(31:nTrial);
    allCond = [allCond; cond];
    
    %% Bias per condition
    % 1 2 3
    % 2 3 1
    % 3 1 2
    if cond(1) == 1
        biasA(iSubject, 1) = mean(bias(1:30));
        biasB(iSubject, 1) = mean(bias(31:60));
        biasN(iSubject, 1) = mean(bias(61:90));
        order(iSubject, 1) = 1;
    elseif cond(1) == 2
        biasA(iSubject, 1) = mean(bias(61:90));
        biasB(iSubject, 1) = mean(bias(1:30));
        biasN(iSubject, 1) = mean(bias(31:60));
        order(iSubject, 1) = 2;
    else
        biasA(iSubject, 1) = mean(bias(31:60));
        biasB(iSubject, 1) = mean(bias(61:90));
        biasN(iSubject, 1) = mean(bias(1:30));
        order(iSubject, 1) = 3;
    end
    
end

data = [allB0/100, allC/100, allA/100, allR/100, allConfi, allColl, allBias, allAcc, allCond];
choice = [allR/100];
subject = allSubject;
meanB0 = mean(allB0)/100;

%% Difficulty ~ Collaboration
difficulty = 0.5 - abs(allC/100 - .5);
dataTable = table(difficulty, allColl, subject);
fitlme(dataTable, 'allColl ~  1 + difficulty + (difficulty|subject)')

%% Baseline model
tic
model_Base = @(PHI,data) model_Base_func(PHI,data);
beta0 = [0, 0];
options = statset('TolX',1e-4, 'MaxIter', 1000, 'Display', 'Iter');
[beta_Base, PSI_Base, stats_Base, B_Base] = nlmefitsa(data,choice,subject,[],model_Base,beta0, 'Options', options, 'LogLikMethod','lin', 'NIterations', [100, 100, 100]);
results(1,1) = stats_Base.logl;
results(1,2) = stats_Base.aic;
results(1,3) = stats_Base.bic;
disp('Model 1');
toc

%% Model RL_Base
tic
model_RL_Base = @(PHI,data) model_RL_Base_func(PHI,data);
beta0 = [beta_Base', 0];
options = statset('TolX',1e-4, 'MaxIter', 1000, 'Display', 'Iter');
[beta_RL_Base, PSI_RL_Base, stats_RL_Base, B_RL_Base] = nlmefitsa(data,choice,subject,[],model_RL_Base,beta0, 'Options', options, 'LogLikMethod','lin', 'NIterations', [100, 100, 100]);
results(2,1) = stats_RL_Base.logl;
results(2,2) = stats_RL_Base.aic;
results(2,3) = stats_RL_Base.bic;
CI1 = beta_RL_Base + tinv(.975, 120)*std(B_RL_Base')'/sqrt(nSubject);
CI2 = beta_RL_Base - tinv(.975, 120)*std(B_RL_Base')'/sqrt(nSubject);
disp('Model 2');
toc

%% Model PH
tic
model_PH = @(PHI,data)  model_RL_PH_func(PHI,data);
beta0 = [beta_Base', 0, 0];
options = statset('TolX',1e-4, 'MaxIter', 1000, 'Display', 'Iter');
[beta_PH, PSI_PH, stats_PH, B_PH] = nlmefitsa(data, choice,subject,[],model_PH,beta0, 'Options', options, 'LogLikMethod','lin', 'NIterations', [100, 100, 100]);
results(3,1) = stats_PH.logl;
results(3,2) = stats_PH.aic;
results(3,3) = stats_PH.bic;
disp('Model 3');
toc

%% Figure - AIC
figure
AIC = round(results(:, 2));
baseline = 8380;
b = barh(AIC + baseline);
xlabel('AIC')
yticklabels({'Baseline model';'Learning model';'Learning model PH'})
b.FaceColor = 'flat';
b.CData(1,:) = [0.8500 0.3250 0.0980];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.4940 0.1840 0.5560];

minX = -5;
maxX = 90;
xlim([-5, 85])
ax = gca; 
ax.XTick = linspace(minX + 15, maxX, 5);
labels = -baseline + minX:20:-baseline + maxX;
ax.XTickLabel = strsplit(num2str(labels + 15)); 
set(gca, 'FontSize', 16)
text(15, 2, ['\leftarrow' '   Best model'], 'FontSize', 20)

%% Figure - Simulated Vs. Real bias
beta_Model = beta_RL_Base;
B_Model = B_RL_Base;

for i = 1:nSubject
    subjectData = data(allSubject == i, :);
    predSubject1 = model_Base_func_pred(beta_Base' + B_Base(:, i)', subjectData);
    predSubject2 = model_RL_Base_func_pred(beta_RL_Base' + B_RL_Base(:, i)', subjectData);
    if subjectData(1, 9) == 1
        biasFinal(i, 1) = mean(subjectData(1:30,4) - subjectData(1:30,2));
        biasFinal(i, 2) = mean(subjectData(31:60,4) - subjectData(31:60,2));
        biasFinal(i, 3) = mean(subjectData(61:90,4) - subjectData(61:90,2));
        biasFinal(i, 4) = mean(predSubject2(1:30) - subjectData(1:30,2));
        biasFinal(i, 5) = mean(predSubject2(31:60) - subjectData(31:60,2));
        biasFinal(i, 6) = mean(predSubject2(61:90) - subjectData(61:90,2));
    elseif subjectData(1, 9) == 2
        biasFinal(i,1) = mean(subjectData(61:90,4) - subjectData(61:90,2));
        biasFinal(i,2) = mean(subjectData(1:30,4) - subjectData(1:30,2));
        biasFinal(i,3) = mean(subjectData(31:60,4) - subjectData(31:60,2));
        biasFinal(i,4) = mean(predSubject2(61:90) - subjectData(61:90,2));
        biasFinal(i,5) = mean(predSubject2(1:30) - subjectData(1:30,2));
        biasFinal(i,6) = mean(predSubject2(31:60) - subjectData(31:60,2));
    elseif subjectData(1, 9) == 3
        biasFinal(i,1) = mean(subjectData(31:60,4) - subjectData(31:60,2));
        biasFinal(i,2) = mean(subjectData(61:90,4) - subjectData(61:90,2));
        biasFinal(i,3) = mean(subjectData(1:30,4) - subjectData(1:30,2));
        biasFinal(i,4) = mean(predSubject2(31:60) - subjectData(31:60,2));
        biasFinal(i,5) = mean(predSubject2(61:90) - subjectData(61:90,2));
        biasFinal(i,6) = mean(predSubject2(1:30) - subjectData(1:30,2));
    end
end

biasFinal = 100*biasFinal;

%% Correlation Actual bias - Predicted bias
figure
plotCorr(biasFinal(:, 1), biasFinal(:, 4), 'Bias (data)', 'Bias (model)', [0 0.4470 0.7410])
figure
plotCorr(biasFinal(:, 2), biasFinal(:, 5), 'Bias (data)', 'Bias (model)', [0.8500 0.3250 0.0980])
figure
plotCorr(biasFinal(:, 3), biasFinal(:, 6), 'Bias (data)', 'Bias (model)', [0.9290 0.6940 0.1250])

%% Figure - Pred by condition
allPredBase = [];
allPred = [];
for i = 1:nSubject
    subjectData = data(allSubject == i, :);
    pred = model_RL_Base_func_pred(beta_RL_Base' + B_RL_Base(:, i)', subjectData);
    allPred = [allPred; pred];
end

%% Figure all trials predictions
figure
hold on
responses = reshape(allR, 90, nSubject);
meanResponse = mean(responses, 2);
predictions = reshape(allPred, 90, nSubject);
meanPredictions = mean(predictions, 2)*100;
p1 = plot(1:90, meanResponse, 'b', 'LineWidth',2)
patch([[1:90] flip([1:90])], [[meanResponse' + 2*std(responses')/sqrt(nSubject)] flip(meanResponse' - 2*std(responses')/sqrt(nSubject))], 'b', 'FaceAlpha',0.25, 'EdgeColor','none')
p2 = plot(1:90,meanPredictions, 'sk', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2)

legend([p1,p2],{'Data'; 'Model'})
xlabel('Trial')
ylabel('Mean response')
set(gca, 'FontSize', 18)
xlim([1,90])


%%%%%%%%%%%%%%%%
%% Functions %%%
%%%%%%%%%%%%%%%%

%% model_Base_func
function y = model_Base_func(PHI,data)

a = PHI(:, 1);
b = PHI(:, 2);
y = a + b.*data(:, 2);

end

%% model_Base_func_pred
function y = model_Base_func_pred(PHI,data)

a = PHI(:, 1);
b = PHI(:, 2);
y = a + b.*data(:, 2);

end

%% model_RL_Base_func
function y = model_RL_Base_func(PHI, data)

b0 = PHI(:, 1);
b1 = PHI(:, 2);
alpha = PHI(:, 3);

c = data(:, 2);
a = data(:, 3);
r = data(:, 4);

p(1, 1) = 0;

for iTrial = 2:90
    PE(iTrial) = a(iTrial - 1) - r(iTrial - 1);
    p(iTrial, 1) = p(iTrial-1, 1) + alpha*PE(iTrial);
end

y = b0 + b1.*data(:, 2) + p;

end

%% model_RL_Base_func_pred
function [y, p] = model_RL_Base_func_pred(PHI, data)

b0 = PHI(:, 1);
b1 = PHI(:, 2);
alpha = PHI(:, 3);

a = data(:, 3);
nTrial = size(data, 1);

p(1, 1) = 0;
y(1, 1) = b0 + b1.*data(1, 2) + p(1,1);

for iTrial = 2:nTrial
    PE(iTrial) = a(iTrial - 1) - y(iTrial - 1);
    p(iTrial, 1) = p(iTrial-1, 1) + alpha*PE(iTrial);
    y(iTrial, 1) = b0 + b1.*data(iTrial, 2) + p(iTrial,1);
end

end

%% model_RL_PH_func
function y = model_RL_PH_func(PHI, data)

b0 = PHI(:, 1);
b1 = PHI(:, 2);
alpha = PHI(:, 3);
eta =  PHI(:, 4);

c = data(:, 2);
a = data(:, 3);
r = data(:, 4);

p(1, 1) = 0;
alphaD(1) = alpha;

for iTrial = 2:90
    PE(iTrial) = a(iTrial - 1) - r(iTrial - 1);
    alphaD(iTrial) = eta*abs(PE(iTrial)) + (1 - eta)*(alphaD(iTrial - 1));
    p(iTrial, 1) = p(iTrial - 1, 1) + alphaD(iTrial)*PE(iTrial);
end
y = b0 + b1.*data(:, 2) + p;
end

%% plotCorr
function plotCorr(X,Y, xLabel, yLabel, colorFig, titleFig)

hold on
[r, p] = corr(X, Y);
axis([-20,50,-20,50])

scatter(X, Y ,100, 'filled' ,'MarkerEdgeColor',colorFig,...
    'MarkerFaceColor',colorFig,...
    'LineWidth',0.5,...
    'MarkerFaceAlpha',.6,...
    'MarkerEdgeAlpha',1)
xlabel(xLabel, 'FontSize', 40);
ylabel(yLabel, 'FontSize', 40);
ax = gca;
ax.XAxis.FontSize = 17;
ax.YAxis.FontSize = 17;


h1 = lsline
set(h1, 'LineWidth', 2.5, 'Color', colorFig)


if p < .001
    title(['\it r\rm\bf = ' num2str(r,2) ',\it p\rm\bf < .001'], 'FontSize', 20);
else
    title(['\it r\rm\bf = ' num2str(r,2) ',\it  p\rm\bf = ' num2str(p,1)], 'FontSize', 20);
end
box on

mdl = fitlm(X, Y);
xl = xlim;
[~, yci] = predict(mdl,[xl(1):.5:xl(2)]');
ciUpper = yci(:,2);
ciLower = yci(:,1);
fill([[xl(1):.5:xl(2)], fliplr([xl(1):.5:xl(2)])], [ciUpper; flip(ciLower)]', colorFig, 'facealpha', .25,'EdgeColor','none');

end