clc; clear all; close all;
rng('default')

%%%%%%%%%%%%%
%% Level 1 %%
%%%%%%%%%%%%%

%% Load data
allDataL1 = readtable('Exp1-Level1.csv');
nSubjectL1 = max(allDataL1.subject);

%% Extract data
for iSubject = 1:nSubjectL1
    
    %% Subject data
    data = allDataL1(allDataL1.subject == iSubject, :);
    nTrial = height(data);
    
    %% Sad response
    responseSadL1(iSubject, 1) = mean(data.responseSad);
    responseSadL1Block(iSubject, 1:4) = mean([data.responseSad(1:25),...
        data.responseSad(26:50),...
        data.responseSad(51:75),...
        data.responseSad(76:100)]);
    
end

%% Main analysis
[p1, d1, CI1] = permutationTest(responseSadL1, .5);

%% Baseline analysis
[p2, d2, CI2] = permutationTest(responseSadL1Block(:,1), mean(responseSadL1Block(:,2:4),2));
meanL1 = mean([responseSadL1Block(:,1), mean(responseSadL1Block(:,2:4),2)]);

%% Single Faces analysis
singleDataL1 = readtable('Exp1-Level1-SingleFace.csv');
modelSingleL1 = fitlme(singleDataL1, 'subjectFaceResponse ~ singleFace + (singleFace|subject)');
[~, ~, stats1] = fixedEffects(modelSingleL1,'DFMethod','satterthwaite','alpha',0.05)

%%%%%%%%%%%%%%%%
%% Level 2 AI %%
%%%%%%%%%%%%%%%%

%% Load data
allDataL3AI = readtable('Exp1-Level3-H-AI-H.csv');
allDataL2AI = allDataL3AI(151:450, [1:14, 18:19]);
allDataL2AI.responseSad = allDataL2AI.responseAISad;
allDataL2AI.responseHappy = 1 - allDataL2AI.responseSad;

[p3, d3, CI3] = permutationTest(responseSadL1, .65);

%%%%%%%%%%%%%%%%
%% Level 3 AI %%
%%%%%%%%%%%%%%%%

%% Load data
allDataL3AI = readtable('Exp1-Level3-H-AI-H.csv');
nSubject = max(allDataL3AI.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allDataL3AI(allDataL3AI.subject == iSubject, :);
    data.Block = [0*ones(150, 1); 1*ones(50, 1); 2*ones(50, 1); 3*ones(50, 1); 4*ones(50, 1); 5*ones(50, 1); 6*ones(50, 1)];
    nTrial = height(data);
    
    %% Sad response - Part 1
    responseBaselineL3AI(iSubject, 1) = mean(data.responseSad(1:150));
    responseBaselineBlockL3AI(iSubject, 1:5) = mean([data.responseSad(1:30),...
        data.responseSad(31:60),...
        data.responseSad(61:90),...
        data.responseSad(91:120),...
        data.responseSad(121:150)]);
    
    %% Sad response - Collaboration
    startTrial = 151;
    endTrial = 200;
    nBlock = 6;
    for iBlock = 1:nBlock
        responseCollaborationL3AI(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
        
    %% Percentage change mind: Agreement vs. Disagreement
    Change_Agreement_L3AI(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3AI(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);
    
    for iBlock = 1:nBlock
        Change_Disagreement_L3AIBlock(iSubject, iBlock) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1 & data.Block == iBlock)/sum(data.responseSad ~= data.responseAISad == 1 & data.Block == iBlock);
    end
    
end

%% Baseline bias
[p4, d4, CI4] = permutationTest(responseBaselineBlockL3AI(:,1), mean(responseBaselineBlockL3AI(:,2:5),2));

%% Trend analysis
responseAI = responseCollaborationL3AI - responseBaselineL3AI;
nSubject = 50;
nBlock = 6;
nItem = nSubject*nBlock;
responseAIL3AI = reshape(responseAI', nItem, 1);
blockL3AI = repmat([1:nBlock]', nSubject, 1);
subjectL3AI = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTableL3AI = table(subjectL3AI, blockL3AI, responseAIL3AI);
trendAnalysisL3AI = fitlme(trendTableL3AI,...
    'responseAIL3AI ~  1 + blockL3AI + (blockL3AI|subjectL3AI)')
[~, ~, stats2] = fixedEffects(trendAnalysisL3AI,'DFMethod','satterthwaite','alpha',0.05)

%% Change mind
[p5, d5, CI5] = permutationTest(Change_Disagreement_L3AI, Change_Agreement_L3AI);

%% L3 AI - effect
[p6, d6, CI6] = permutationTest(mean(responseCollaborationL3AI, 2) - responseBaselineL3AI, 0);

%%%%%%%%%%%%%%%%%%%
%% Level 2 Human %%
%%%%%%%%%%%%%%%%%%%

%% Load data
allDataL2 = readtable('Exp1-Level2-H-H-H.csv');
nSubject = max(allDataL2.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allDataL2(allDataL2.subject == iSubject, :);
    nTrial = height(data);
    
    %% Sad response
    responseSadL2(iSubject, 1) = nanmean(data.responseSad);
    
end

%% Human - Level 2 bias
[p7, d7, CI7] = permutationTest(responseSadL2, 0.5);
[p8, d8, CI8] = permutationTest(responseSadL2, responseSadL1);

%% Figure
figure
hold on
x = responseSadL1;
y = mean(allDataL2AI.responseSad);
z = responseSadL2;

% Colors
colorL1 = [0.4375, 0.6758, 0.2773];
colorL2AI = [0, 0.4375, 0.75];
colorL2Human = [0.9258, 0.4883, 0.1914];

% Error bars
errorbar(1, mean(x), std(x)/sqrt(size(x, 1)), '.','MarkerSize', 80,'LineWidth', 1, 'CapSize', 0, 'Color', colorL1, 'MarkerEdgeColor', colorL1, 'MarkerFaceColor', colorL1);
errorbar(2, mean(y), 0, '.','MarkerSize', 80,'LineWidth', 1, 'CapSize', 0, 'Color', colorL2AI, 'MarkerEdgeColor', colorL2AI, 'MarkerFaceColor', colorL2AI);
errorbar(2, mean(z), std(z)/sqrt(size(z, 1)), '.','MarkerSize', 80,'LineWidth', 1, 'CapSize', 0, 'Color', colorL2Human, 'MarkerEdgeColor', colorL2Human, 'MarkerFaceColor', colorL2Human);
e1 = scatter(1, mean(x), 180, colorL1, 'filled');
e2 = scatter(2, mean(y), 180, colorL2AI, 'filled');
e3 = scatter(2, mean(z), 180, colorL2Human, 'filled');

% Figure parameters
xlim([0.5, 2.5]);
ylim([0.5, 0.7])
ylabel('P(More sad)')
xticks([1:2])
xticklabels({'Layer 1'; 'Layer 2'})
text(1,mean(x) + .03, 'Human','HorizontalAlignment', 'center','FontSize', 13);
text(2,mean(y) + .03, 'AI','HorizontalAlignment', 'center','FontSize', 13);
text(2,mean(z) + .03, 'Human','HorizontalAlignment', 'center','FontSize', 13);
set(gca, 'FontSize', 16)
b = mean(y) - mean(x);
a = mean(y) - 2*b;
myarrow([1.1, 1.9], b.*[1.1, 1.9] + a)
b = mean(z) - mean(x);
a = mean(z) - 2*b;
myarrow([1.1, 1.9], b.*[1.1, 1.9] + a)

%%%%%%%%%%%%%%%%%%%
%% Level 3 Human %%
%%%%%%%%%%%%%%%%%%%

%% Load data
allDataL3Human = readtable('Exp1-Level3-H-H-H.csv');
nSubject = max(allDataL3Human.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allDataL3Human(allDataL3Human.subject == iSubject, :);
    data.Block = [0*ones(150, 1); 1*ones(50, 1); 2*ones(50, 1); 3*ones(50, 1); 4*ones(50, 1); 5*ones(50, 1); 6*ones(50, 1)];
    nTrial = height(data);
    
    %% Sad response - Part 1
    responseBaselineL3Human(iSubject, 1) = mean(data.responseSad(1:150));
    responseBaselineBlockL3Human(iSubject, 1:5) = mean([data.responseSad(1:30),...
        data.responseSad(31:60),...
        data.responseSad(61:90),...
        data.responseSad(91:120),...
        data.responseSad(121:150)]);
    
    %% Sad response - Collaboration
    startTrial = 151;
    endTrial = 200;
    nBlock = 6;
    for iBlock = 1:nBlock
        responseCollaborationL3Human(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    %% Percentage change mind: Agreement vs. Disagreement
    Change_Agreement_L3Human(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3Human(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);
    
    for iBlock = 1:nBlock
        Change_Disagreement_L3HumanBlock(iSubject, iBlock) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1 & data.Block == iBlock)/sum(data.responseSad ~= data.responseAISad == 1 & data.Block == iBlock);
    end
    
end

%% Figure - Part 1
[p9, d9, CI9] = permutationTest(mean(responseCollaborationL3Human, 2) - responseBaselineL3Human, 0);
mean(responseCollaborationL3Human - responseBaselineL3Human)
mean(mean(responseCollaborationL3Human - responseBaselineL3Human))

%% Baseline bias
[p10, d10, CI10] = permutationTest(responseBaselineBlockL3Human(:,1), mean(responseBaselineBlockL3Human(:,2:5),2));
mean([responseBaselineBlockL3Human(:,1), mean(responseBaselineBlockL3Human(:,2:5),2)])

%% Trend analysis
responseHuman = responseCollaborationL3Human - responseBaselineL3Human;
nSubject = 50;
nBlock = 6;
nItem = nSubject*nBlock;
responseHumanL3Human = reshape(responseHuman', nItem, 1);
blockL3Human = repmat([1:nBlock]', nSubject, 1);
subjectL3Human = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTableL3Human = table(subjectL3Human, blockL3Human, responseHumanL3Human);
trendAnalysisL3Human = fitlme(trendTableL3Human,...
    'responseHumanL3Human ~  1 + blockL3Human + (blockL3Human|subjectL3Human)')
[~,~,stats3] = fixedEffects(trendAnalysisL3Human,'DFMethod','satterthwaite','alpha',0.05)

%% Change mind
[p11, d11, CI11] = permutationTest(Change_Disagreement_L3Human, Change_Agreement_L3Human);

%%%%%%%%%%%%%%%%%%%%%
%% Level 3 HumanAI %%
%%%%%%%%%%%%%%%%%%%%%

%% Load data
allDataL3HumanAI = readtable('Exp1-Level3-H-AIH-H.csv');
nSubject = max(allDataL3HumanAI.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allDataL3HumanAI(allDataL3HumanAI.subject == iSubject, :);
    data.Block = [0*ones(150, 1); 1*ones(50, 1); 2*ones(50, 1); 3*ones(50, 1); 4*ones(50, 1); 5*ones(50, 1); 6*ones(50, 1)];
    nTrial = height(data);
    
    %% Sad response - Part 1
    responseBaselineL3HumanAI(iSubject, 1) = mean(data.responseSad(1:150));
    responseBaselineBlockL3HumanAI(iSubject, 1:5) = mean([data.responseSad(1:30),...
        data.responseSad(31:60),...
        data.responseSad(61:90),...
        data.responseSad(91:120),...
        data.responseSad(121:150)]);
    
    %% Sad response - Collaboration
    startTrial = 151;
    endTrial = 200;
    nBlock = 6;
    for iBlock = 1:nBlock
        responseCollaborationL3HumanAI(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    %% Percentage change mind: Agreement vs. Disagreement
    Change_Agreement_L3HumanAI(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3HumanAI(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);
    
    for iBlock = 1:nBlock
        Change_Disagreement_L3HumanAIBlock(iSubject, iBlock) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1 & data.Block == iBlock)/sum(data.responseSad ~= data.responseAISad == 1 & data.Block == iBlock);
    end
    
end

%% Baseline comparison
[p12, d12, CI1] = permutationTest(mean(responseCollaborationL3HumanAI, 2) - responseBaselineL3HumanAI, 0);
mean(responseCollaborationL3HumanAI - responseBaselineL3HumanAI)
mean(mean(responseCollaborationL3HumanAI - responseBaselineL3HumanAI))

%% Trend analysis
responseHumanAI = responseCollaborationL3HumanAI - responseBaselineL3HumanAI;
nSubject = max(allDataL3HumanAI.subject);
nBlock = 6;
nItem = nSubject*nBlock;
responseHumanAIL3HumanAI = reshape(responseHumanAI', nItem, 1);
blockL3HumanAI = repmat([1:nBlock]', nSubject, 1);
subjectL3HumanAI = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTableL3HumanAI = table(subjectL3HumanAI, blockL3HumanAI, responseHumanAIL3HumanAI);
trendAnalysisL3HumanAI = fitlme(trendTableL3HumanAI,...
    'responseHumanAIL3HumanAI ~  1 + blockL3HumanAI + (blockL3HumanAI|subjectL3HumanAI)')
[~,~,stats4] = fixedEffects(trendAnalysisL3HumanAI,'DFMethod','satterthwaite','alpha',0.05)

%% Change mind
[p13, d13, CI13] = permutationTest(Change_Disagreement_L3HumanAI, Change_Agreement_L3HumanAI);

%%%%%%%%%%%%%%%%%%%%%%
%% Level 3 HumanAI2 %%
%%%%%%%%%%%%%%%%%%%%%%

%% Load data
allDataL3HumanAI2 = readtable('Exp1-Level3-H-HAI-H.csv');
nSubject = max(allDataL3HumanAI2.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allDataL3HumanAI2(allDataL3HumanAI2.subject == iSubject, :);
    data.Block = [0*ones(150, 1); 1*ones(50, 1); 2*ones(50, 1); 3*ones(50, 1); 4*ones(50, 1); 5*ones(50, 1); 6*ones(50, 1)];
    nTrial = height(data);
    
    %% Sad response - Part 1
    responseBaselineL3HumanAI2(iSubject, 1) = mean(data.responseSad(1:150));
    responseBaselineBlockL3HumanAI2(iSubject, 1:5) = mean([data.responseSad(1:30),...
        data.responseSad(31:60),...
        data.responseSad(61:90),...
        data.responseSad(91:120),...
        data.responseSad(121:150)]);
    
    %% Sad response - Collaboration
    startTrial = 151;
    endTrial = 200;
    nBlock = 6;
    for iBlock = 1:6
        responseCollaborationL3HumanAI2(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    %% Percentage change mind: Agreement vs. Disagreement
    
    Change_Agreement_L3HumanAI2(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3HumanAI2(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);
    
    for iBlock = 1:nBlock
        Change_Disagreement_L3HumanAI2Block(iSubject, iBlock) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1 & data.Block == iBlock)/sum(data.responseSad ~= data.responseAISad == 1 & data.Block == iBlock);
    end
    
end

%% Figure - Part 1
[p14, d14, CI14] = permutationTest(mean(responseCollaborationL3HumanAI2, 2) - responseBaselineL3HumanAI2, 0);
mean(responseCollaborationL3HumanAI2 - responseBaselineL3HumanAI2)
mean(mean(responseCollaborationL3HumanAI2 - responseBaselineL3HumanAI2))

%% Baseline compairosn
[p15, d15, CI15] = permutationTest(mean(responseCollaborationL3HumanAI2, 2) - responseBaselineL3HumanAI2, 0);

%% Trend analysis
responseHumanAI2 = responseCollaborationL3HumanAI2 - responseBaselineL3HumanAI2;
nSubject = max(allDataL3HumanAI2.subject);
nBlock = 6;
nItem = nSubject*nBlock;
responseHumanAIL3HumanAI2 = reshape(responseHumanAI2', nItem, 1);
blockL3HumanAI2 = repmat([1:nBlock]', nSubject, 1);
subjectL3HumanAI2 = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTableL3HumanAI2 = table(subjectL3HumanAI2, blockL3HumanAI2, responseHumanAIL3HumanAI2);
trendAnalysisL3HumanAI2 = fitlme(trendTableL3HumanAI2,...
    'responseHumanAIL3HumanAI2 ~  1 + blockL3HumanAI2 + (blockL3HumanAI2|subjectL3HumanAI2)')
[~,~,stats5] = fixedEffects(trendAnalysisL3HumanAI2,'DFMethod','satterthwaite','alpha',0.05)

%% Change mind
[p16, d16, CI16] = permutationTest(Change_Disagreement_L3HumanAI2, Change_Agreement_L3HumanAI2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline - No interaction %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
allDataBaseline = readtable('Exp1-Level3-Baseline.csv');
nSubject = max(allDataBaseline.subject);

%% Extract data
for iSubject = 1:nSubject
    
    %% Subject data
    data = allDataBaseline(allDataBaseline.subject == iSubject, :);
    nTrial = height(data);
    
    %% Sad response - Part 1
    responseBaselineL3Baseline(iSubject, 1) = mean(data.responseSad(1:150));
    responseBaselineBlockL3Baseline(iSubject, 1:5) = mean([data.responseSad(1:30),...
        data.responseSad(31:60),...
        data.responseSad(61:90),...
        data.responseSad(91:120),...
        data.responseSad(121:150)]);
    
    %% Sad response - Collaboration
    startTrial = 151;
    endTrial = 200;
    nBlock = 6;
    for iBlock = 1:nBlock
        responseCollaborationL3Baseline(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        
        Change_Total_Baseline_Block(iSubject, iBlock) = sum(data.changeMind(startTrial:endTrial) == 1)/50;
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    % Percentage change mind: Agreement vs. Disagreement
    Change_Total_Baseline(iSubject, 1) = sum(data.changeMind == 1)/300;
    
end

%% Remove outliers
Change_Total_Baseline_Clean = Change_Total_Baseline(abs(zscore(Change_Total_Baseline)) < 3);

%% Comparisons to other conditions
[p17, d17, CI17] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Disagreement_L3AI,2));
[p18, d18, CI18] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Disagreement_L3Human,2));
[p19, d19, CI19] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Disagreement_L3HumanAI,2));
[p20, d20, CI20] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Disagreement_L3HumanAI2,2));

[p21, d21, CI21] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Agreement_L3AI,2));
[p22, d22, CI22] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Agreement_L3Human,2));
[p23, d23, CI23] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Agreement_L3HumanAI,2));
[p24, d24, CI24] = permutationTest2(Change_Total_Baseline_Clean, mean(Change_Agreement_L3HumanAI2,2));

%% Trend analysis
nSubject = 50;
nBlock = 6;
nItem = nSubject*nBlock;
responseChangeMind = reshape(Change_Total_Baseline_Block, nItem, 1);
blockL3Baseline = repmat([1:nBlock]', nSubject, 1);
subjectL3Baseline = reshape(repmat([1:nSubject], nBlock, 1), nItem, 1);
trendTableL3Baseline = table(subjectL3Baseline, blockL3Baseline, responseChangeMind);
trendAnalysisL3Baseline = fitlme(trendTableL3Baseline,...
    'responseChangeMind ~  1 + blockL3Baseline + (blockL3Baseline|subjectL3Baseline)')
[~,~,stats6] = fixedEffects(trendAnalysisL3Baseline,'DFMethod','satterthwaite','alpha',0.05)

%%%%%%%%%%%%%%%%%%
%% Human Vs. AI %%
%%%%%%%%%%%%%%%%%%

figure
hold on
x = 1:6;

% Colors
colorAI = [0, 0.4375, 0.75];
colorHuman = [0.9258, 0.4883, 0.1914];
colorHumanAI = [0.4033, 0.43995, 0.5586];
colorHumanAI2 = [0.80765, 0.441375, 0.56745];

% Responses
responseAI = responseCollaborationL3AI - responseBaselineL3AI;
responseHuman = responseCollaborationL3Human - responseBaselineL3Human;

% Human-AI2
y = mean(responseHumanAI2);
plot4 = plot(x, y, '-o', 'Color', colorHumanAI2, 'LineWidth', 2, 'MarkerSize',4 , 'MarkerFaceColor', colorHumanAI2)
xconf = [x x(end:-1:1)] ;
yconf = [y - std(responseHumanAI2)/sqrt(size(responseHumanAI2, 1)) y(end:-1:1) + fliplr(std(responseHumanAI2)/sqrt(size(responseHumanAI2, 1)))];
f33 = fill(xconf, yconf,'red');
f33.FaceColor = colorHumanAI2;
f33.EdgeColor = 'none';
f33.FaceAlpha = 0.2;
line([0.5, 6.5], [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')

% Human-AI
y = mean(responseHumanAI);
plot2 = plot(x, y, '-o', 'Color', colorHumanAI, 'LineWidth', 2, 'MarkerSize',4 , 'MarkerFaceColor', colorHumanAI)
xconf = [x x(end:-1:1)] ;
yconf = [y - std(responseHumanAI)/sqrt(size(responseHumanAI, 1)) y(end:-1:1) + fliplr(std(responseHumanAI)/sqrt(size(responseHumanAI, 1)))];
f33 = fill(xconf, yconf,'red');
f33.FaceColor = colorHumanAI;
f33.EdgeColor = 'none';
f33.FaceAlpha = 0.2;
line([0.5, 6.5], [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')

% Human
y = mean(responseHuman);
plot3 = plot(x, y, '-o', 'Color', colorHuman, 'LineWidth', 2, 'MarkerSize',4 , 'MarkerFaceColor', colorHuman)
xconf = [x x(end:-1:1)] ;
yconf = [y - std(responseHuman)/sqrt(size(responseHuman, 1)) y(end:-1:1) + fliplr(std(responseHuman)/sqrt(size(responseHuman, 1)))];
f22 = fill(xconf, yconf,'red');
f22.FaceColor = colorHuman;
f22.EdgeColor = 'none';
f22.FaceAlpha = 0.2;
line([0.5, 6.5], [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')
y = mean(responseHuman);
plot3 = plot(x, y, '-o', 'Color', colorHuman, 'LineWidth', 2, 'MarkerSize',4 , 'MarkerFaceColor', colorHuman)
xconf = [x x(end:-1:1)] ;
yconf = [y - std(responseHuman)/sqrt(size(responseHuman, 1)) y(end:-1:1) + fliplr(std(responseHuman)/sqrt(size(responseHuman, 1)))];
f22 = fill(xconf, yconf,'red');
f22.FaceColor = colorHuman;
f22.EdgeColor = 'none';
f22.FaceAlpha = 0.2;
line([0.5, 6.5], [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')

% AI
y = mean(responseAI);
xconf = [x x(end:-1:1)] ;
plot1 = plot(x, y,'-o' ,'Color', colorAI, 'LineWidth', 2, 'MarkerSize',4)
yconf = [y - std(responseAI)/sqrt(size(responseAI, 1)) y(end:-1:1) + fliplr(std(responseAI)/sqrt(size(responseAI, 1)))];
f11 = fill(xconf, yconf,'red');
f11.FaceColor = colorAI;
f11.EdgeColor = 'none';
f11.FaceAlpha = 0.2;
x = 1:6;
y = mean(responseAI);
xconf = [x x(end:-1:1)] ;
plot1 = plot(x, y,'-o' ,'Color', colorAI, 'LineWidth', 2, 'MarkerSize',4)
yconf = [y - std(responseAI)/sqrt(size(responseAI, 1)) y(end:-1:1) + fliplr(std(responseAI)/sqrt(size(responseAI, 1)))];
f11 = fill(xconf, yconf,'red');
f11.FaceColor = colorAI;
f11.EdgeColor = 'none';
f11.FaceAlpha = 0.2;

% Properties
axis([0.75 6.25 -0.03 0.15])
ylabel({'Induced bias'})
set(gca, 'FontSize', 16)
legend([plot1, plot2, plot4, plot3],{'Human-AI'; 'Human-AI perceived as human'; 'Human-Human perceived as AI'; 'Human-Human'}, 'FontSize', 14,  'EdgeColor', 'none','Color', 'none')
xlabel('\bf Level 3  Blocks \rm', 'FontSize', 20)

%%%%%%%%%%%%%%%%%%%
%% Summary Table %%
%%%%%%%%%%%%%%%%%%%

% Generate table
effectHumanAI = responseCollaborationL3HumanAI - responseBaselineL3HumanAI;
effectHumanAI2 = responseCollaborationL3HumanAI2 - responseBaselineL3HumanAI2;
effectHuman = responseCollaborationL3Human - responseBaselineL3Human;
effectAI = responseCollaborationL3AI - responseBaselineL3AI;
allEffect = [effectAI; effectHumanAI; effectHumanAI2; effectHuman];
codeGroup = [ones(50,1); 2*ones(50,1); 3*ones(50, 1); 4*ones(50, 1)];
actualAI = zeros(length(codeGroup), 1);
actualAI(codeGroup == 1 | codeGroup == 2) = 1;
purportedlyAI = zeros(length(codeGroup), 1);
purportedlyAI(codeGroup == 1 | codeGroup == 3) = 1;
changeDecision = [Change_Disagreement_L3AIBlock; Change_Disagreement_L3HumanAIBlock; Change_Disagreement_L3HumanAI2Block; Change_Disagreement_L3HumanBlock];
allEffect = [allEffect, changeDecision, codeGroup, actualAI, purportedlyAI];

% Save table
inducedBiasNames = arrayfun(@(x) sprintf('InducedBias%d', x), 1:6, 'UniformOutput', false);
changeDecisionNames = arrayfun(@(x) sprintf('ChangeDecision%d', x), 1:6, 'UniformOutput', false);
variableNames = [inducedBiasNames, changeDecisionNames, {'Condition', 'Input', 'Label'}];
dataTable = array2table(allEffect, 'VariableNames', variableNames);
writetable(dataTable, 'Summary_Table.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indifference point analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[IP1, CI_1, estimation1] = IPCI(allDataL1);
[IP2, CI_2, estimation2] = IPCI(allDataL2AI);
[IP3, CI_3, estimation3] = IPCI(allDataL3AI(~isnan(allDataL3AI.changeMind), :));
[IP4, CI_4, estimation4] = IPCI(allDataL3Human(~isnan(allDataL3Human.changeMind), :));
[IP5, CI_5, estimation5] = IPCI(allDataL3HumanAI(~isnan(allDataL3HumanAI.changeMind), :));
[IP6, CI_6, estimation6] = IPCI(allDataL3HumanAI2(~isnan(allDataL3HumanAI2.changeMind), :));

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

%% permutationTest 2
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
d = round((mean(x) - mean(y)) / sqrt(((sizeX - 1)*std(x)^2 + sizeY*std(y)^2) / (sizeX + sizeY - 2)), 2);

%% CI
meanSample = zeros(nPerm, 1);
for iPerm = 1:nPerm
    meanSample(iPerm, 1) = mean(x(randi(sizeX, sizeX, 1))) - mean(y(randi(sizeY, sizeY, 1)));
end

CI = round(quantile(meanSample, [.025, .975]), 2);

end

%% IPCI
function [IP, CI, estimation] = IPCI(dataTable)

tic
% Normalize
q = dataTable{:,2:13};
dataTable.meanEvidence = mean(q,2);
dataTable.meanEvidence = (dataTable.meanEvidence - 25.5)/24.5;
N = height(dataTable);
nIter = 1000;
estimation1 = zeros(nIter, 1);
estimation2 = zeros(nIter, 1);

for iIter = 1:nIter
    dataTableRandom = dataTable(randi(N, N, 1), :);
    glmeBaseline = fitglme(dataTableRandom,...
        'responseHappy ~  1 + meanEvidence + (meanEvidence|subject)',...
        'Distribution','Binomial','Link','logit');
    estimation1(iIter, 1) = glmeBaseline.Coefficients.Estimate(1);
    estimation2(iIter, 1) = glmeBaseline.Coefficients.Estimate(2);
end

estimation = -estimation1./estimation2;
estimation = sort(estimation);
CI(1) = quantile(estimation, 0.025);
CI(2) = quantile(estimation, 0.975);

glmeBaseline = fitglme(dataTable,...
    'responseHappy ~  1 + meanEvidence + (meanEvidence|subject)',...
    'Distribution','Binomial','Link','logit');
IP = -glmeBaseline.Coefficients.Estimate(1)/glmeBaseline.Coefficients.Estimate(2);
toc

end

%% myarrow
function myarrow(x,y)
ax = gca;
axpos = get(ax, 'Position');
X = get(gca,'XLim');
Y = get(gca,'YLim');
annotation('arrow',[axpos(1) + (x(1) - X(1))/(X(2) - X(1))*axpos(3),...
    axpos(1) + (x(2) - X(1))/(X(2) - X(1))*axpos(3)],...
    [axpos(2) + (y(1) - Y(1))/(Y(2) - Y(1))*axpos(4),...
    axpos(2) + (y(2) - Y(1))/(Y(2) - Y(1))*axpos(4)]);
end