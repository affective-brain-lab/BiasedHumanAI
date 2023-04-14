clc; clear all; close all;
rng('default')

%% Iterate through the experiments
for iExp = 1:2
    
    %% Load data
    if iExp == 1
        allData = readtable('Exp2-Accurate.csv');
    else
        allData = readtable('Exp2-Biased.csv');
    end
    nSubject = max(allData.subject);
    
    %% Iterate through the subjects
    for iSubject = 1:nSubject
        
        data = allData(allData.subject == iSubject, :);
        
        %% Bias & Accuracy
        bias = data.response - data.evidence;
        accuracy = abs(data.response - data.evidence);
        
        %% Blocks
        for iBlock = 1:6
            if iExp == 1
                accuracyBlock(iSubject, iBlock) = mean(accuracy(data.block == iBlock));
                
            elseif iExp == 2
                biasBlock(iSubject, iBlock) = mean(bias(data.block == iBlock));
                
            end
        end
    end
end

%% Measures
meanBiasB = biasBlock(:,2:6) - biasBlock(:,1);
meanAccuracyA = -accuracyBlock(:,2:6) + accuracyBlock(:,1);

%% Bias
figure
hold on
color = [0.9961, 0.5430, 0.0039];
color2 = [0.9961, 0.9492, 0.8945];
accurateColor = [0, 0.5313, 0.8633];
biasedColor = [0.9961, 0.5430, 0.0039];
noisyColor = [.7608, 0.2314, 0.1294];
x = 1:5;
y = mean(meanBiasB);
xconf = [x x(end:-1:1)] ;
yconf = [y - std(meanBiasB)/sqrt(size(meanBiasB, 1)) y(end:-1:1) + fliplr(std(meanBiasB)/sqrt(size(meanBiasB, 1)))];
p1 = errorbar(1:5, y, std(meanBiasB)/sqrt(size(meanBiasB, 1)),'-o','MarkerSize', 15,'LineWidth', 1.5, 'Color', color, 'MarkerEdgeColor', color, 'MarkerFaceColor', color2, 'CapSize', 0)
line([0.5, 6.5], [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')
axis([0.75 5.25 -1 8.5])
xlabel('Block of interaction with a Biased AI')
ylabel({'Induced bias'})
set(gca, 'FontSize', 16)
[p1, d1, CI1] = permutationTest(mean(meanBiasB,2), 0);

%% Trend analysis bias
nSubject = 50;
nBlock = 5;
nItem = nSubject*nBlock;
response = reshape(meanBiasB', nItem, 1);
block = repmat([1:nBlock]', nSubject, 1);
subject = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTable = table(subject, block, response);
trendAnalysis = fitlme(trendTable,...
    'response ~  1 + block + (block|subject)')

%% Accuracy
figure
hold on
color = [0, 0.5313, 0.8633];
color2 = [0.8945    0.9492    0.9805];
x = 1:5;
y = mean(meanAccuracyA);
xconf = [x x(end:-1:1)] ;
yconf = [y - std(meanBiasB)/sqrt(size(meanAccuracyA, 1)) y(end:-1:1) + fliplr(std(meanAccuracyA)/sqrt(size(meanAccuracyA, 1)))];
errorbar(1:5, y, std(meanAccuracyA)/sqrt(size(meanAccuracyA, 1)),'-o','MarkerSize', 15,'LineWidth', 1.5, 'Color', color, 'MarkerEdgeColor', color, 'MarkerFaceColor', color2, 'CapSize', 0)
line([0.5, 6.5], [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')
axis([0.75 5.25 -1 6])
xlabel('Block of interaction with an Accurate AI')
ylabel({'Induced bias'})
set(gca, 'FontSize', 16)
[p2, d2, CI2] = permutationTest(mean(meanAccuracyA,2), 0);

%% Trend analysis bias
nSubject = 50;
nBlock = 5;
nItem = nSubject*nBlock;
response = reshape(meanAccuracyA', nItem, 1);
block = repmat([1:nBlock]', nSubject, 1);
subject = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTable = table(subject, block, response);
trendAnalysis = fitlme(trendTable,...
    'response ~  1 + block + (block|subject)')

%%%%%%%%%%%%%%%
%% Functions %%
%%%%%%%%%%%%%%%

%% permutationTest
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