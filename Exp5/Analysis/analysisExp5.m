clc; close all; clear all;
rng('default')

%% Read data
dataNew = readtable('Exp5.csv');

%% Extract choices
nSubject = 200;
for iSubject= 1:nSubject
    dataSubject = dataNew(dataNew.subject == iSubject, :);
    for iChoice = 1:6
        choice(iSubject, iChoice) = sum(dataSubject.choice == iChoice);
    end
    condition(iSubject, 1) = dataSubject.groupAI(1);
end

groupAI = choice(condition == 1, :);
nSubjectAI = length(groupAI);
groupControl = choice(condition == 0, :);
nSubjectControl = length(groupControl);

%% Create Figure
figure
hold on

%% Parameters
color1 = [0, 70, 120]/256;
color2 = [173, 216, 230]/256;
color3 = [200, 53, 48]/256;
color4 = [244, 129, 97]/256;
color5 = [58, 58, 58]/256;
color6 = [190, 190, 190]/256;
width = 1.5;

%% Group means
groupMeans = [mean(groupAI(:, 1)),mean( groupControl(:, 1)),...
    mean(groupAI(:, 2)), mean(groupControl(:, 2)),...
    mean( groupAI(:, 3)),mean( groupControl(:, 3)),...
    mean(groupAI(:, 4)), mean(groupControl(:, 4)),...
    mean( groupAI(:, 5)), mean(groupControl(:, 5)),...
    mean( groupAI(:, 6)), mean(groupControl(:, 6))];

%% Means bar
xlabelValues = [1,1.5, 2.5,3, 4,4.5, 5.5,6, 7,7.5, 8.5,9];
b = bar(xlabelValues,groupMeans);
b.FaceColor = 'flat';
b.CData(1,:) = color1;
b.CData(2,:) = color1;
b.CData(3,:) = color2;
b.CData(4,:) = color2;
b.CData(5,:) = color3;
b.CData(6,:) = color3;
b.CData(7,:) = color4;
b.CData(8,:) = color4;
b.CData(9,:) = color5;
b.CData(10,:) = color5;
b.CData(11,:) = color6;
b.CData(12,:) = color6;
b.FaceAlpha = 0.675;

%% Error bars
errorbar(xlabelValues(1), mean(groupAI(:,1)), std(groupAI(:,1))/sqrt(nSubjectAI), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(2), mean(groupControl(:,1)), std(groupControl(:,1))/sqrt(nSubjectControl), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(3), mean(groupAI(:,2)), std(groupAI(:,2))/sqrt(nSubjectAI), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(4), mean(groupControl(:,2)), std(groupControl(:,2))/sqrt(nSubjectControl), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(5), mean(groupAI(:,3)), std(groupAI(:,3))/sqrt(nSubjectAI), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(6), mean(groupControl(:,3)), std(groupControl(:,3))/sqrt(nSubjectControl), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(7), mean(groupAI(:,4)), std(groupAI(:,4))/sqrt(nSubjectAI), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(8), mean(groupControl(:,4)), std(groupControl(:,5))/sqrt(nSubjectControl), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(9), mean(groupAI(:,5)), std(groupAI(:,5))/sqrt(nSubjectAI), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(10), mean(groupControl(:,5)), std(groupControl(:,5))/sqrt(nSubjectControl), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(11), mean(groupAI(:,6)), std(groupAI(:,6))/sqrt(nSubjectAI), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)
errorbar(xlabelValues(12), mean(groupControl(:,6)), std(groupControl(:,6))/sqrt(nSubjectControl), '.k', 'MarkerFaceColor', 'k', 'MarkerSize', 15,'LineWidth', width, 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'CapSize', 0)

%% Figure settings
set(gca, 'FontSize', 16)
xlim([0.5, 9.5])
ylim([0,42])
ylabel('Chosen as a Financial Manager (%)')
xtickangle(45);
set(gca, 'FontSize', 16)
b9 = bar(55,1);
b9.FaceColor = 'flat';
b9.CData(1,:) = [256,256,256]/256;
b10 = bar(55,1);
b10.FaceColor = 'flat';
b10.CData(1,:) = [256,256,256]/256;
set(gcf, 'Position', [89.6667   41.6667  586.6667  599.3333])
xlabelValues = [1,1.25,1.5, 2.5,2.75,3, 4,4.25,4.5, 5.5,5.75,6, 7,7.25,7.5, 8.5,8.75,9];
xticks(xlabelValues)
xticklabels({'','White Man', '',...
    '','White Woman', '',...
    '','Asian Man', '',...
    '','Asian Woman', '',...
    '','Black Man', '',...
    '','Black Woman', ''});


%% Add stripes
maxValue = mean(groupControl);
xValue = [1.5, 3,4.5, 6, 7.5,9];

for jGroup = 1:6
    for iBar = 1:maxValue(jGroup) + 3
        x1 = xValue(jGroup) - 0.2;
        x2 = xValue(jGroup) + 0.2;
        x3 = xValue(jGroup) - 0.2;
        x4 = xValue(jGroup) + 0.2;
        y1 =-1.5 + iBar - 1;
        y2 = 0 + iBar - 1;
        y3 = -1.4 + iBar - 1;
        y4 = 0.1 + iBar - 1;
        
        if y2 > maxValue(jGroup)
            x2 = xValue(jGroup) - 0.2 + (maxValue(jGroup) - y1)/(y2 - y1)*0.4;
            y2 = maxValue(jGroup);
            
        end
        
        if y4 > maxValue(jGroup)
            x4 = xValue(jGroup) - 0.2 + (maxValue(jGroup) - y3)/(y4 - y3)*0.4;
            
            y4 = maxValue(jGroup);
            
        end
        
        if y1 > maxValue(jGroup)
            x1 = 100; x2 = 100;
            y1 = 100; y2 = 100;
            
            
        end
        
        if y3 > maxValue(jGroup)
            x1 = 100; x2 = 100;
            y3 = 100; y4 = 100;
            
        end
        
        line([x1,x2], [y1,y2], 'Color', 'k', 'LineWidth', 1)
        line([x3,x4], [y3,y4], 'Color', 'k', 'LineWidth', 1)
        set(get(gca, 'XAxis'), 'TickLength', [0 0])
    end
end

legend([b9,b10],{'Stable Diffusion Interaction Group', 'Control Group'}, 'FontSize', 14)

%% t-tests
% White Men
[h1, p1] = ttest2(groupAI(:,1), groupControl(:,1));
% White Women
[h2, p2] = ttest2(groupAI(:,2), groupControl(:,2));
% Asian Men
[h3, p3] = ttest2(groupAI(:,3), groupControl(:,3));
% Asian Women
[h4, p4] = ttest2(groupAI(:,4), groupControl(:,4));
% Black Men
[h5, p5] = ttest2(groupAI(:,5), groupControl(:,5));
% Black Women
[h6, p6] = ttest2(groupAI(:,6), groupControl(:,6));