clc; clear all;  close all;
rng('default')

%%%%%%%%%%%%%%
%% Relative %%
%%%%%%%%%%%%%%

%% Load data
allDataRelative = readtable('allDataRelative.csv');
nSubject = max(allDataRelative.subjectRelative);
N = nSubject;

%% Iterate through the subjects
for iSubject = 1:nSubject
 
    %% Subject data
    data = allDataRelative(allDataRelative.subjectRelative == iSubject, :);
      
    %% Recoding and normalizing
    data.response(data.type == 1) = 8 - data.response(data.type == 1);
    data.norm(data.type == 1) = 8 - data.norm(data.type == 1);
    
    %% Bias relative
    block = [ones(45,1); 2*ones(45,1); 3*ones(45,1); 4*ones(45,1)];
    for iBlock = 1:4
        bias(iSubject, iBlock) = mean(data.response(block == iBlock) - data.norm(block == iBlock));        
    end
         
end

%% Bias relative
figure
hold on
value = bias(:,2:4) - bias(:,1);
meanBias = mean(bias(:,2:4) - bias(:,1),2);
[p1, d1, CI1] = permutationTest(meanBias, 0);
errorbar(1:3, mean(bias(:,2:4)) - mean(bias(:,1)), std(value/sqrt(N)),'-o', 'LineWidth', 2, 'Color', [0 112 192]/256) % , 'MarkerFaceColor', 'b'
xticklabels({'1'; '2';'3'})
xlabel({'Interaction block'})
xticks([1,2,3])
ylabel({'AI induced bias';'\fontsize{18} (Relative 2 - Relative 1)'})                               
xlim([.5 3.5])
ylim([ -.5, 0.5])
title('Relative Evaluation')
line([0 3.5], [0,0], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2)
text(3.58, 0, ' \rightarrow   Bias toward men', 'FontSize', 15, 'HorizontalAlignment','left', 'Color', 'k', 'Rotation', 90)
text(3.58, 0, ' Bias toward women  \leftarrow', 'FontSize', 15, 'HorizontalAlignment','right', 'Color', 'k', 'Rotation', 90)
set(gca,'FontSize', 20)

%%%%%%%%%%%%%%
%% Separate %%
%%%%%%%%%%%%%%

%% Load data
allDataSeparate = readtable('allDataSeparate.csv');
nSubject = max(allDataSeparate.subjectSeparate);
N = nSubject;

%% Iterate through the subjects
for iSubject = 1:nSubject
 
    %% Subject data
    data = allDataSeparate(allDataSeparate.subjectSeparate == iSubject, :);
      
    %% Bias separate
    block = [ones(10,1); 2*ones(10,1)];
    biasSeparate(iSubject, 1) = mean(data.rating2(block == 1 & data.rating2_Person < 6));
    biasSeparate(iSubject, 3) = mean(data.rating2(block == 1 & data.rating2_Person > 5));
    biasSeparate(iSubject, 2) = mean(data.rating2(block == 2 & data.rating2_Person < 6));
    biasSeparate(iSubject, 4) = mean(data.rating2(block == 2 & data.rating2_Person > 5));
    biasSeparate(iSubject, 5) = biasSeparate(iSubject, 2) - biasSeparate(iSubject, 1);
    biasSeparate(iSubject, 6) = biasSeparate(iSubject, 4) - biasSeparate(iSubject, 3);
    biasSeparate(iSubject, 7) = biasSeparate(iSubject, 2) - biasSeparate(iSubject, 1) - biasSeparate(iSubject, 4) + biasSeparate(iSubject, 3);
       
end

%% Bias separate
figure
hold on

[p2, d2, CI2] = permutationTest(biasSeparate(:,7), 0);
[p3, d3, CI3] = permutationTest(biasSeparate(:,6), 0);
[p4, d4, CI4] = permutationTest(biasSeparate(:,5), 0);
[~,~,~,stats] = ttest(biasSeparate(:,7));
errorbar(1:2, [mean(biasSeparate(:,5)), mean(biasSeparate(:,6))], [stats.sd/sqrt(N), stats.sd/sqrt(N)],'.-', 'LineWidth', 2, 'Color', [0 176 80]/256)
xticklabels({'Men'; 'Women'})
xticks([1,2])
ylabel({'AI induced bias'; '\fontsize{18} (Separate 2 - Separate 1)'})
xlabel({'Evaluation of'})
xlim([.5 2.5])
ylim([-.5, .5])
title('Separate Evaluation')
text(1.5, .35, '**','FontSize', 40,'HorizontalAlignment', 'center')
line([.5 2.5], [0,0], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2)
text(2.58, 0, ' \rightarrow   Positive bias', 'FontSize', 15, 'HorizontalAlignment','left', 'Color', 'k', 'Rotation', 90)
text(2.58, 0, ' Negative bias  \leftarrow', 'FontSize', 15, 'HorizontalAlignment','right', 'Color', 'k', 'Rotation', 90)
set(gca,'FontSize', 20)

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
