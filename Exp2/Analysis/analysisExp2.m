clc; clear all; close all;
rng('default')

%% Load data
allData = readtable('Exp2.csv');
nSubject = max(allData.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allData(allData.subject == iSubject, :);
    nTrial = height(data);
    
    %% Evidence
    evidence = data.evidence;
    
    %% Response
    response = data.response;
    
    %% Response AI
    responseAI = data.responseAI;
    
    %% Confidence
    confidence = data.confidence;
    confidence = (confidence - min(confidence)) ./ (max(confidence) - min(confidence));
    
    %% Collaboration
    collaboration = data.collaboration;
    
    %% Condition
    condition = data.condition;
    
    %% Bias baseline
    biasBase(iSubject, 1) = mean(response(condition == 0) - evidence(condition == 0));
    biasA(iSubject, 1) = mean(response(condition == 1) - evidence(condition == 1));
    biasB(iSubject, 1) = mean(response(condition == 2) - evidence(condition == 2));
    biasN(iSubject, 1) = mean(response(condition == 3) - evidence(condition == 3));
    
    %% Accuracy
    errorBase(iSubject, 1) = mean(abs(evidence(condition == 0) - response(condition == 0)));
    errorA(iSubject, 1) = mean(abs(evidence(condition == 1) - response(condition == 1)));
    errorB(iSubject, 1) = mean(abs(evidence(condition == 2) - response(condition == 2)));
    errorN(iSubject, 1) = mean(abs(evidence(condition == 3) - response(condition == 3)));
    
    %% Collaboration
    collaborationBase(iSubject, 1) = 0;
    collaborationA(iSubject, 1) = mean(collaboration(condition == 1));
    collaborationB(iSubject, 1) = mean(collaboration(condition == 2));
    collaborationN(iSubject, 1) = mean(collaboration(condition == 3));
    
    %% Confidence
    confidenceBase(iSubject, 1) = mean(confidence(condition == 0));
    confidenceA(iSubject, 1) = mean(confidence(condition == 1));
    confidenceB(iSubject, 1) = mean(confidence(condition == 2));
    confidenceN(iSubject, 1) = mean(confidence(condition == 3));
    confidenceAll(iSubject, 1) = mean(confidence);
    
    %% Order
    order(iSubject, 1) = data.order(1);
    
    %% Influence
    zBiasTemp = zscore([biasA(iSubject, 1), biasB(iSubject, 1), biasN(iSubject, 1)]);
    zAccTemp = zscore([-errorA(iSubject, 1), -errorB(iSubject, 1), -errorN(iSubject, 1)]);
    zBias(iSubject, 1) = zBiasTemp(2);
    zAcc(iSubject, 1) = zAccTemp(1);
    if order(iSubject, 1) == 1
        perceivedInfluenceA(iSubject, 1) = mean(data.perceivedInfluenceA);
        perceivedInfluenceB(iSubject, 1)= mean(data.perceivedInfluenceB);
        perceivedInfluenceN(iSubject, 1)= mean(data.perceivedInfluenceN);
    elseif order(iSubject, 1) == 2
        perceivedInfluenceA(iSubject, 1) = mean(data.perceivedInfluenceN);
        perceivedInfluenceB(iSubject, 1)= mean(data.perceivedInfluenceA);
        perceivedInfluenceN(iSubject, 1)= mean(data.perceivedInfluenceB);
    else
        perceivedInfluenceA(iSubject, 1) = mean(data.perceivedInfluenceB);
        perceivedInfluenceB(iSubject, 1)= mean(data.perceivedInfluenceN);
        perceivedInfluenceN(iSubject, 1)= mean(data.perceivedInfluenceA);
    end
    
end

%% Accuracy
createFigure(errorBase, errorA, errorB, errorN, [-15, 15], 'AI Induced accuracy')

%% Bias
createFigure(biasBase, biasA, biasB, biasN, [-20, 20], 'AI Induced bias')

%% Actual influence
figure
plot([zBias, zAcc]', '-o', 'Color', [.8, .8, .8])
hold on
plot(mean([zAcc, zBias]),'-ok', 'LineWidth', 3, 'MarkerFaceColor', 'k')
xticks([1,2])
xticklabels({'Accurate AI', 'Biased AI'})
xlabel({'Interaction with'})
ylabel({'Actual Influence (Z-scored)'})
xlim([0,3])
ylim([-1.5,1.5])
set(gca, 'FontSize', 16)
text(1.5, .925, '\it p\rm = 0.90', 'FontSize', 14, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
line([1.1,1.9], [.75, .75], 'Color','k', 'LineWidth', 2, 'LineStyle', '--')

%% Actual influence 2
actualBias = (biasBase - biasB)./abs(biasBase);
actualAccuracy = -(errorB - errorBase)./abs(errorB);

%% Perceived influence
figure
influence = [perceivedInfluenceA, perceivedInfluenceB, perceivedInfluenceN];
plot(influence', '-o', 'Color', [.8, .8, .8])
hold on
plot(mean(influence),'-ok', 'LineWidth', 3, 'MarkerFaceColor', 'k')
xticks([1,2,3])
xticklabels({'Accurate AI', 'Biased AI', 'Noisy AI'})
xlabel({'Interaction with'})
ylabel({'Perceived influence'; '(participants rating)'})
xlim([0,4])
set(gca, 'FontSize', 16)
text(1.5, 5.325, '\it p\rm <0.001', 'FontSize', 14, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
text(2, 6.325, '\it p\rm <0.001', 'FontSize', 14, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
text(2.5, 5.325, '\it p\rm  = 0.11', 'FontSize', 14, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
line([1.1,1.9], [5, 5], 'Color','k', 'LineWidth', 1.5, 'LineStyle', '--')
line([1.1,2.9], [6, 6], 'Color','k', 'LineWidth', 1.5, 'LineStyle', '--')
line([2.1,2.9], [5, 5], 'Color','k', 'LineWidth', 1.5, 'LineStyle', '--')

%% Comparisons
[p1, d1, CI1] = permutationTest(allData.evidence(1:30), 50);
[p2, d2, CI2] = permutationTest(biasBase, 0);
[p3, d3, CI3] = permutationTest(biasA, biasBase);
[p4, d4, CI4] = permutationTest(biasB, biasBase);
[p5, d5, CI5] = permutationTest(biasN, biasBase);
[p6, d6, CI6] = permutationTest(biasB, biasA);
[p7, d7, CI7] = permutationTest(biasB, biasN);
[p8, d8, CI8] = permutationTest(biasA, biasN);
[p9, d9, CI9] = permutationTest(errorA, errorBase);
[p10, d10, CI10] = permutationTest(errorB, errorBase);
[p11, d11, CI11] = permutationTest(errorN, errorBase);
[p12, d12, CI12] = permutationTest(errorA, errorN);
[p13, d13, CI13] = permutationTest(errorB, errorN);
[p14, d14, CI14] = permutationTest(errorA, errorB);
[p15, d15, CI15] = permutationTest(collaborationA, collaborationBase);
[p16, d16, CI16] = permutationTest(collaborationB, collaborationBase);
[p17, d17, CI17] = permutationTest(collaborationN, collaborationBase);
[p18, d18, CI18] = permutationTest(collaborationA, collaborationN);
[p19, d19, CI19] = permutationTest(collaborationB, collaborationN);
[p20, d20, CI20] = permutationTest(collaborationN, collaborationB);
[p21, d21, CI21] = permutationTest(confidenceA, confidenceBase);
[p22, d22, CI22] = permutationTest(confidenceB, confidenceBase);
[p23, d23, CI23] = permutationTest(confidenceN, confidenceBase);
[p24, d24, CI24] = permutationTest(confidenceA, confidenceN);
[p25, d25, CI25] = permutationTest(confidenceB, confidenceN);
[p26, d26, CI26] = permutationTest(confidenceN, confidenceB);
[p27, d27, CI27] = permutationTest(influence(:, 1), influence(:, 2));
[p28, d28, CI28] = permutationTest(influence(:, 2), influence(:, 3));
[p29, d29, CI29] = permutationTest(influence(:, 1), influence(:, 3));
[p30, d30, CI30] = permutationTest(zAcc, zBias);
[p31, d31, CI31] = permutationTest(actualAccuracy, actualBias);

%% Regressions
allData.bias = allData.response - allData.evidence;
allData.acc = abs(allData.response - allData.evidence);
allData.difficuly = abs(allData.evidence - 50);
allData.time = repmat([1:30]', 4*120, 1);
allData(allData.condition == 0, :) = [];

% Bias coding
allData.conditionAvsB = zeros(length(allData.condition), 1);
allData.conditionNvsB = zeros(length(allData.condition), 1);
allData.conditionAvsB(allData.condition == 1, :) = 1;
allData.conditionNvsB(allData.condition == 3, :) = 1;

% Accuracy coding
allData.conditionBvsA = zeros(length(allData.condition), 1);
allData.conditionNvsA = zeros(length(allData.condition), 1);
allData.conditionBvsA(allData.condition == 2, :) = 1;
allData.conditionNvsA(allData.condition == 3, :) = 1;

% Resgression bias
mdlBias = fitlme(allData, 'bias ~ conditionAvsB + conditionNvsB + evidence + time + (conditionAvsB|subject) + (conditionNvsB|subject) + (evidence|subject) + (time|subject)')
[~,~,stats1] = fixedEffects(mdlBias,'DFMethod','satterthwaite','alpha',0.05)

% Regression acc
mdlAcc = fitlme(allData, 'acc ~ conditionBvsA + conditionNvsA + evidence + time + (conditionBvsA|subject) + (conditionNvsA|subject) + (evidence|subject) + (time|subject)')
[~,~,stats2] = fixedEffects(mdlAcc,'DFMethod','satterthwaite','alpha',0.05)

% Resgression difficulty
mdlDifficulty = fitlme(allData, 'collaboration ~ difficuly + (difficuly|subject)')
[~,~,stats3] = fixedEffects(mdlDifficulty,'DFMethod','satterthwaite','alpha',0.05)

%%%%%%%%%%%%%%%
%% Functions %%
%%%%%%%%%%%%%%%

%% Figure
function createFigure(b, x, y, z, limY, labelY)

figure
hold on
data = [x, y, z] - b;

switch labelY
    case 'AI Induced accuracy'
        data = b - [x, y, z];
end
%% Colors
accurateColor = [0, 0.5313, 0.8633];
biasedColor = [0.9961, 0.5430, 0.0039];
noisyColor = [.7608, 0.2314, 0.1294];

%% Error
plot([1, 2, 3]' + -.025 + 0.05*repmat(rand(1, size(data, 1)), 3, 1), data','-o' ,'Color', [.7, .7, .7, .2], 'LineWidth', .1, 'MarkerSize', 3);
errorBar = std(data)/sqrt(size(data,1));
errorbar([1], mean(data(:,1)), errorBar(1), 'o', 'LineWidth', 2, 'Color', accurateColor, 'MarkerFaceColor', accurateColor)
errorbar([2], mean(data(:,2)), errorBar(2), 'o', 'LineWidth', 2, 'Color', biasedColor , 'MarkerFaceColor', biasedColor )
errorbar([3], mean(data(:,3)), errorBar(3), 'o', 'LineWidth', 2, 'Color', noisyColor, 'MarkerFaceColor', noisyColor)

%% Quantiles
coordLineStyle = '';
h = boxplot([data(:, 1), data(:, 2), data(:,3)], 'Symbol', coordLineStyle, 'Colors' , [accurateColor; biasedColor; noisyColor], 'Whisker' , 1, 'Widths',0.15);
set(h,'LineWidth',1)
xlim([0, 4])
ylim([limY])

ylabel(labelY)
set(gca, 'FontSize', 16)

xticks([1, 2, 3])
xticklabels({'Accurate AI','Biased AI', 'Noisy AI'})
xlabel('Interaction with')

%% Violin
rangeViolin = linspace(limY(1), limY(2), 100);
[c1,v]=ksdensity(data(:, 1), rangeViolin, 'Bandwidth', (limY(2) - limY(1))/60);
[c2,v]=ksdensity(data(:, 2), rangeViolin, 'Bandwidth',(limY(2) - limY(1))/60);
[c3,v]=ksdensity(data(:, 3), rangeViolin, 'Bandwidth',(limY(2) - limY(1))/60);
scale = 0.15/max([c1,c2,c3]);
plot(1+scale*((c1)),rangeViolin,'Color', [accurateColor, 0.5])
plot(1-scale*((c1)), rangeViolin,'Color', [accurateColor, 0.5])
plot(2+scale*((c2)),rangeViolin,'Color', [biasedColor, 0.5])
plot(2-scale*((c2)), rangeViolin,'Color', [biasedColor, 0.5])
plot(3+scale*((c3)),rangeViolin,'Color', [noisyColor, 0.5])
plot(3-scale*((c3)), rangeViolin,'Color', [noisyColor, 0.5])
fill([1+scale*((c1)), fliplr(1-scale*((c1)))], [rangeViolin, fliplr(rangeViolin)],accurateColor,'FaceAlpha',.1,'EdgeColor','none');
fill([2+scale*((c2)), fliplr(2-scale*((c2)))], [rangeViolin, fliplr(rangeViolin)],biasedColor,'FaceAlpha',.1,'EdgeColor','none');
fill([3+scale*((c3)), fliplr(3-scale*((c3)))], [rangeViolin, fliplr(rangeViolin)],noisyColor,'FaceAlpha',.1,'EdgeColor','none');

%% Significance
yl = ylim;
box on
space = abs(yl(1) - yl(2))/7;
yLine1 = max(mean(data)) + space;
line([1.1 1.9], [yLine1, yLine1], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '--')
[sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,y);
text(1.5, yLine1 + yHeight, ['\it p \rm' sigSymbol], 'FontSize', sizeFont, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')

yLine1 = max(mean(data)) + space;
line([2.1 2.9], [yLine1, yLine1], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '--')
[sigSymbol, sizeFont, yHeight] = sigSymbolFunc(y,z);
text(2.5, yLine1 + yHeight, ['\it p \rm' sigSymbol], 'FontSize', sizeFont, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')

yLine2 = max(mean(data)) + 2*space;
line([1.1 2.9], [yLine2, yLine2], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '--')
[sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,z);
text(2, yLine2 + yHeight, ['\it p \rm' sigSymbol], 'FontSize', sizeFont, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')

if nargin > 5
    line([0, 4], [0, 0], 'Color', 'r', 'LineStyle', '--')
    [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,b);
    text(1.425, yHeight, ['\it p \rm' sigSymbolFunc(x,b)], 'FontSize', sizeFont, 'LineStyle', '--', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r')
    [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(y,b);
    text(2.425, yHeight, ['\it p \rm' sigSymbolFunc(y,b)], 'FontSize', sizeFont, 'LineStyle', '--', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r')
    [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(z,b);
    text(3.425, yHeight,  ['\it p \rm' sigSymbolFunc(z,b)], 'FontSize', sizeFont, 'LineStyle', '--', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r')
end
    function [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,y)
        [p] = permutationTest(x,y);
        sizeFont = 32;
        yHeight = 0.25*space;
        if p > .05
            sigSymbol = 'n.s.';
            sizeFont = 12;
            yHeight = 0.25*space;
        elseif p > .01 & p<= .05
            sigSymbol = '*';
        elseif p > .001 & p<= .01
            sigSymbol = '**';
        else
            sigSymbol = '***';
        end
        if p < 0.001
         sigSymbol = ['< 0.001']; 
        elseif p < 0.01
              sigSymbol = [ '= ' num2str(p, 1)];
        else
        sigSymbol = [ '= ' num2str(p, 2)];
        end
        sizeFont = 12;
    end

end

function [p, d, CI] = permutationTest(x, y)

%% Permutation test
sizeData = size(x, 1);
nPerm = 1e5;
meanDiff = zeros(sizeData, 1);

if size(y, 1) == 1
    y = repmat(y, sizeData, 1);
end

for iPerm = 1:nPerm
    signs = rand(sizeData,1) > .5;
    meanDiff(iPerm) = nanmean(signs .* x + -1 .* signs .* y);
end

actualMean = nanmean(x - y);
pLeft = nanmean(actualMean <= sort(meanDiff));
pRight = nanmean(actualMean >= sort(meanDiff));
p = 2*min([pLeft,pRight]);

%% Effect size
d = round(nanmean(x - y)/ nanstd(x - y), 2);

%% CI
diff = x - y;
meanSample = zeros(sizeData, 1);
for iPerm = 1:nPerm
    meanSample(iPerm, 1) = mean(diff(randi(sizeData, sizeData, 1)));
end
CI = round(quantile(meanSample, [.025, .975]), 2);

end