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
ylabel('Actual Influence (z-scored)')
xlim([0,3])
ylim([-1.5,1.5])
set(gca, 'FontSize', 16)
text(1.5, .9, 'n.s.', 'FontSize', 16, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
line([1.1,1.9], [.75, .75], 'Color','k', 'LineWidth', 2)

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
ylabel('Perceived influence (participant rating)')
xlim([0,4])
set(gca, 'FontSize', 16)
text(1.5, 5.3, '***', 'FontSize', 32, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
text(2, 6.3, '***', 'FontSize', 32, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
text(2.5, 5.3, 'n.s.', 'FontSize', 16, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
line([1.1,1.9], [5, 5], 'Color','k', 'LineWidth', 1.5)
line([1.1,2.9], [6, 6], 'Color','k', 'LineWidth', 1.5)
line([2.1,2.9], [5, 5], 'Color','k', 'LineWidth', 1.5)

%% Comparisons
[p, d, CI] = permutationTest(allData.evidence(1:30), 50);
[p, d, CI] = permutationTest(biasBase, 0);
[p, d, CI] = permutationTest(biasA, biasBase);
[p, d, CI] = permutationTest(biasB, biasBase);
[p, d, CI] = permutationTest(biasN, biasBase);
[p, d, CI] = permutationTest(biasB, biasA);
[p, d, CI] = permutationTest(biasB, biasN);
[p, d, CI] = permutationTest(biasA, biasN);
[p, d, CI] = permutationTest(errorA, errorBase);
[p, d, CI] = permutationTest(errorB, errorBase);
[p, d, CI] = permutationTest(errorN, errorBase);
[p, d, CI] = permutationTest(errorA, errorN);
[p, d, CI] = permutationTest(errorB, errorN);
[p, d, CI] = permutationTest(errorA, errorB);
[p, d, CI] = permutationTest(collaborationA, collaborationBase);
[p, d, CI] = permutationTest(collaborationB, collaborationBase);
[p, d, CI] = permutationTest(collaborationN, collaborationBase);
[p, d, CI] = permutationTest(collaborationA, collaborationN);
[p, d, CI] = permutationTest(collaborationB, collaborationN);
[p, d, CI] = permutationTest(collaborationN, collaborationB);
[p, d, CI] = permutationTest(confidenceA, confidenceBase);
[p, d, CI] = permutationTest(confidenceB, confidenceBase);
[p, d, CI] = permutationTest(confidenceN, confidenceBase);
[p, d, CI] = permutationTest(confidenceA, confidenceN);
[p, d, CI] = permutationTest(confidenceB, confidenceN);
[p, d, CI] = permutationTest(confidenceN, confidenceB);
[p, d, CI] = permutationTest(influence(:, 1), influence(:, 2));
[p, d, CI] = permutationTest(influence(:, 2), influence(:, 3));
[p, d, CI] = permutationTest(influence(:, 1), influence(:, 3));
[p, d, CI] = permutationTest(zAcc, zBias);
[p, d, CI] = permutationTest(actualAccuracy, actualBias);

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
line([1.1 1.9], [yLine1, yLine1], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '-')
[sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,y);
text(1.5, yLine1 + yHeight, sigSymbol, 'FontSize', sizeFont, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')

yLine1 = max(mean(data)) + space;
line([2.1 2.9], [yLine1, yLine1], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '-')
[sigSymbol, sizeFont, yHeight] = sigSymbolFunc(y,z);
text(2.5, yLine1 + yHeight, sigSymbol, 'FontSize', sizeFont, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')

yLine2 = max(mean(data)) + 2*space;
line([1.1 2.9], [yLine2, yLine2], 'LineWidth', .5, 'Color', 'k', 'LineStyle', '-')
[sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,z);
text(2, yLine2 + yHeight, sigSymbol, 'FontSize', sizeFont, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')

if nargin > 5
    line([0, 4], [0, 0], 'Color', 'r', 'LineStyle', '--')
    [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,b);
    text(1.3, yHeight, sigSymbolFunc(x,b), 'FontSize', sizeFont, 'LineStyle', '--', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r')
    [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(y,b);
    text(2.3, yHeight, sigSymbolFunc(y,b), 'FontSize', sizeFont, 'LineStyle', '--', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r')
    [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(z,b);
    text(3.3, yHeight, sigSymbolFunc(z,b), 'FontSize', sizeFont, 'LineStyle', '--', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r')
end
    function [sigSymbol, sizeFont, yHeight] = sigSymbolFunc(x,y)
        [~, p] = ttest(x,y);
        sizeFont = 32;
        yHeight = 0.17*space;
        if p > .05
            sigSymbol = 'n.s.';
            sizeFont = 12;
            yHeight = 0.22*space;
        elseif p > .01 & p<= .05
            sigSymbol = '*';
        elseif p > .001 & p<= .01
            sigSymbol = '**';
        else
            sigSymbol = '***';
        end
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
