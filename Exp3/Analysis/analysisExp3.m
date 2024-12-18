clc; close all; clear all;
rng('default')

% Value Labels for Demographic Categories: 
% 1 = Black Woman
% 2 = Black Man
% 3 = Asian Woman
% 4 = Asian Man
% 5 = White Woman
% 6 = White Man

data = readtable('Exp3.csv');
nSubject = 50;

%% AI
dataAI = data(data.control == 0, :);
for iSubject = 1:nSubject
dataSubject = dataAI(dataAI.subject == iSubject, :);   
for iChoice = 1:6
summaryTableBefore(iSubject, iChoice) = sum(dataSubject.choice == iChoice & dataSubject.condition == 0)/100;
summaryTableAfter(iSubject, iChoice) = sum(dataSubject.choice == iChoice & dataSubject.condition == 1)/100;
end
end
summaryTableBefore = fliplr(summaryTableBefore);
summaryTableAfter = fliplr(summaryTableAfter);
differenceAI = 100*(summaryTableAfter - summaryTableBefore);

createFigure(differenceAI)

%% Control
dataControl = data(data.control == 1, :);
dataControl.subject = dataControl.subject - 100;
for iSubject = 1:nSubject
dataSubject =  dataControl(dataControl.subject == iSubject, :);   
for iChoice = 1:6
summaryTableBefore(iSubject, iChoice) = sum(dataSubject.choice == iChoice & dataSubject.condition == 0)/100;
summaryTableAfter(iSubject, iChoice) = sum(dataSubject.choice == iChoice & dataSubject.condition == 1)/100;
end
end
summaryTableBefore = fliplr(summaryTableBefore);
summaryTableAfter = fliplr(summaryTableAfter);
differenceControl = 100*(summaryTableAfter - summaryTableBefore);

createFigure(differenceControl)

%% Between group comparison for white men
[p1, d1, CI1] = permutationTest2(differenceAI(:, 1), differenceControl(:, 1));

% See SPSS analysis for the Multinomial logistic regression results

%%%%%%%%%%%%%%%
%% Functions %%
%%%%%%%%%%%%%%%

function [p, d, CI] = permutationTest2(x, y)

%% Permutation test
sizeX = size(x, 1);
sizeY = size(y, 1);

nPerm = 1e5;
meanDiff = zeros(sizeX + sizeY, 1);
allData = [x; y];

for iPerm = 1:nPerm
    Perm = randperm(sizeX + sizeY);
    meanDiff(iPerm) = mean(allData(Perm(1:sizeX))) - mean(allData(Perm(sizeX + 1:end)));
end

actualMean = mean(x) - mean(y);
pLeft = nanmean(actualMean <= sort(meanDiff));
pRight = nanmean(actualMean >= sort(meanDiff));
p = 2*min([pLeft,pRight]);

%% Effect size
d = round((mean(x) - mean(y)) / sqrt(((sizeX - 1)*std(x).^2 + sizeY*std(y).^2) / (sizeX + sizeY - 2)), 2);

%% CI
meanSample = zeros(nPerm, 1);
for iPerm = 1:nPerm
    meanSample(iPerm, 1) = mean(x(randi(sizeX, sizeX, 1))) - mean(y(randi(sizeY, sizeY, 1)));
end
CI = round(quantile(meanSample, [.025, .975]), 2);

end

function createFigure(difference)
figure
hold on

% Colors
color1 = [0, 70, 120]/256;
color2 = [173, 216, 230]/256;
color3 = [200, 53, 48]/256;
color4 = [244, 129, 97]/256;
color5 = [58, 58, 58]/256;
color6 = [190, 190, 190]/256;
width = 1.5;
nSubject = length(difference);

% Define the colors for each condition
colors = {color1, color2, color3, color4, color5, color6};

% Loop through each condition (1 to 6)
for i = 1:6
    % Plot error bars without markers ('none' marker type suppresses the marker)
    errorbar(i, mean(difference(:,i)), std(difference(:,i))/sqrt(nSubject), 'LineWidth', width, 'Color', colors{i}, 'CapSize', 0, 'Marker', 'none')
    hold on
end

% Overlay transparent markers after plotting error bars
for i = 1:6
    scatter(i, mean(difference(:,i)), 275, 'w', 'filled', 'MarkerEdgeColor', 'w')
    scatter(i, mean(difference(:,i)), 275, colors{i}, 'filled', 'MarkerEdgeColor', colors{i}, 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 1, 'LineWidth', width - 0.5)
end

% Release hold
hold off

b.FaceColor = 'flat';
b.CData(1,:) = color1;
b.CData(2,:) = color2;
b.CData(3,:) = color3;
b.CData(4,:) = color4;
b.CData(5,:) = color5;
b.CData(6,:) = color6;
b.FaceAlpha = 0.675;

line([0, 7], [0, 0], 'LineStyle', '--', 'Color', 'r')

xticks(1:6)
xticklabels({'White Man','White Woman','Asian Man', 'Asian Woman', 'Black Man', 'Black Woman'})
set(gca, 'FontSize', 16)
xlim([0.5, 6.5])
ylim([-4.5,8.5])
ylabel('Choice: Stage 3 - Stage 1(%)')
xtickangle(45);

set(gcf, 'Position', [89.6667   41.6667  574.0000  599.3333])
set(gca, 'FontSize', 20)

end