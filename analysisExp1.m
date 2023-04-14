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

%%%%%%%%%%%%%
%% Level 2 %%
%%%%%%%%%%%%%

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

%% Figure 
[p3, d3, CI3] = permutationTest(responseSadL2, 0.5);

%%%%%%%%%%%%%%%%
%% Level 2 AI %%
%%%%%%%%%%%%%%%%

%% Load data
allDataL3AI = readtable('Exp1-Level3-H-AI-H.csv');
allDataL2AI = allDataL3AI(151:450, [1:14, 18:19]);
allDataL2AI.responseSad = allDataL2AI.responseAISad;
allDataL2AI.responseHappy = 1 - allDataL2AI.responseSad;

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
    for iBlock = 1:6
        responseCollaborationL3AI(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    %% Percentage type trial - 2 (Human sad/Human happy) X 2 (AI sad/AI happy)
    HS_AIS_L3AI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 1 & ~isnan(data.changeMind))/300;
    HS_AIH_L3AI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 0 & ~isnan(data.changeMind))/300;
    HH_AIS_L3AI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 1 & ~isnan(data.changeMind))/300;
    HH_AIH_L3AI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 0 & ~isnan(data.changeMind))/300;
    
    %% Percentage change mind - 2 (Human sad/Human happy) X 2 (AI sad/AI happy)
    Change_HS_AIS_L3AI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 1 & data.changeMind == 1)/sum(data.responseSad == 1 & data.responseAISad == 1 & ~isnan(data.changeMind));
    Change_HS_AIH_L3AI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 0 & data.changeMind == 1)/sum(data.responseSad == 1 & data.responseAISad == 0 & ~isnan(data.changeMind));
    Change_HH_AIS_L3AI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 1 & data.changeMind == 1)/sum(data.responseSad == 0 & data.responseAISad == 1 & ~isnan(data.changeMind));
    Change_HH_AIH_L3AI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 0 & data.changeMind == 1)/sum(data.responseSad == 0 & data.responseAISad == 0 & ~isnan(data.changeMind));
    
    %% Percentage change mind: Agreement vs. Disagreement 
    Change_Agreement_L3AI(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3AI(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);

    %% Total Sad (after change of mind)
    totalSadL3AI(iSubject, 1) = HS_AIS_L3AI(iSubject, 1) + HS_AIH_L3AI(iSubject, 1)...
        - Change_HS_AIS_L3AI(iSubject, 1)*HS_AIS_L3AI(iSubject, 1)...
        - Change_HS_AIH_L3AI(iSubject, 1)*HS_AIH_L3AI(iSubject, 1)...
        + Change_HH_AIS_L3AI(iSubject, 1)*HH_AIS_L3AI(iSubject, 1)...
        + Change_HH_AIH_L3AI(iSubject, 1)*HH_AIH_L3AI(iSubject, 1);
    
end

%% Figure - Part 1
[p4, d4, CI4] = permutationTest(mean(responseCollaborationL3AI, 2) - responseBaselineL3AI, 0);

%% Baseline bias
[p5, d5, CI5] = permutationTest(responseBaselineBlockL3AI(:,1), mean(responseBaselineBlockL3AI(:,2:5),2));

%% Colaboration
mean([HS_AIS_L3AI, HS_AIH_L3AI, HH_AIS_L3AI, HH_AIH_L3AI]);
mean([Change_HS_AIS_L3AI, Change_HS_AIH_L3AI, Change_HH_AIS_L3AI, Change_HH_AIH_L3AI]);

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

%% Change mind
[p6, d6, CI6] = permutationTest(Change_Disagreement_L3AI, Change_Agreement_L3AI);

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
    for iBlock = 1:6
        responseCollaborationL3Human(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    %% Percentage type trial - 2 (Human sad/Human happy) X 2 (Human sad/Human happy)
    HS_HumanS_L3Human(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 1 & ~isnan(data.changeMind))/300;
    HS_HumanH_L3Human(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 0 & ~isnan(data.changeMind))/300;
    HH_HumanS_L3Human(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 1 & ~isnan(data.changeMind))/300;
    HH_HumanH_L3Human(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 0 & ~isnan(data.changeMind))/300;
    
    %% Percentage change mind - 2 (Human sad/Human happy) X 2 (Human sad/Human happy)
    Change_HS_HumanS_L3Human(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 1 & data.changeMind == 1)/sum(data.responseSad == 1 & data.responseAISad == 1 & ~isnan(data.changeMind));
    Change_HS_HumanH_L3Human(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 0 & data.changeMind == 1)/sum(data.responseSad == 1 & data.responseAISad == 0 & ~isnan(data.changeMind));
    Change_HH_HumanS_L3Human(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 1 & data.changeMind == 1)/sum(data.responseSad == 0 & data.responseAISad == 1 & ~isnan(data.changeMind));
    Change_HH_HumanH_L3Human(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 0 & data.changeMind == 1)/sum(data.responseSad == 0 & data.responseAISad == 0 & ~isnan(data.changeMind));
    
    %% Percentage change mind: Agreement vs. Disagreement 
    Change_Agreement_L3Human(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3Human(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);

    %% Total Sad (after change of mind)
    totalSadL3Human(iSubject, 1) = HS_HumanS_L3Human(iSubject, 1) + HS_HumanH_L3Human(iSubject, 1)...
        - Change_HS_HumanS_L3Human(iSubject, 1)*HS_HumanS_L3Human(iSubject, 1)...
        - Change_HS_HumanH_L3Human(iSubject, 1)*HS_HumanH_L3Human(iSubject, 1)...
        + Change_HH_HumanS_L3Human(iSubject, 1)*HH_HumanS_L3Human(iSubject, 1)...
        + Change_HH_HumanH_L3Human(iSubject, 1)*HH_HumanH_L3Human(iSubject, 1);
    
end

%% Figure - Part 1
[p7, d7, CI7] = permutationTest(mean(responseCollaborationL3Human, 2) - responseBaselineL3Human, 0);
mean(responseCollaborationL3Human - responseBaselineL3Human)
mean(mean(responseCollaborationL3Human - responseBaselineL3Human))

%% Baseline bias
[p8, d8, CI8] = permutationTest(responseBaselineBlockL3Human(:,1), mean(responseBaselineBlockL3Human(:,2:5),2));
mean([responseBaselineBlockL3Human(:,1), mean(responseBaselineBlockL3Human(:,2:5),2)])

%% Colaboration
mean([HS_HumanS_L3Human, HS_HumanH_L3Human, HH_HumanS_L3Human, HH_HumanH_L3Human])
mean([Change_HS_HumanS_L3Human, Change_HS_HumanH_L3Human, Change_HH_HumanS_L3Human, Change_HH_HumanH_L3Human])

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

%% Change mind
[p9, d9, CI9] = permutationTest(Change_Disagreement_L3Human, Change_Agreement_L3Human);

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
    for iBlock = 1:6
        responseCollaborationL3HumanAI(iSubject, iBlock) = mean(data.responseSad(startTrial:endTrial));
        startTrial = startTrial + 50;
        endTrial = endTrial + 50;
    end
    
    %% Percentage type trial - 2 (HumanAI sad/HumanAI happy) X 2 (HumanAI sad/HumanAI happy)
    HS_HumanAIS_L3HumanAI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 1 & ~isnan(data.changeMind))/300;
    HS_HumanAIH_L3HumanAI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 0 & ~isnan(data.changeMind))/300;
    HH_HumanAIS_L3HumanAI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 1 & ~isnan(data.changeMind))/300;
    HH_HumanAIH_L3HumanAI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 0 & ~isnan(data.changeMind))/300;
    
    %% Percentage change mind - 2 (HumanAI sad/HumanAI happy) X 2 (HumanAI sad/HumanAI happy)
    Change_HS_HumanAIS_L3HumanAI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 1 & data.changeMind == 1)/sum(data.responseSad == 1 & data.responseAISad == 1 & ~isnan(data.changeMind));
    Change_HS_HumanAIH_L3HumanAI(iSubject, 1) = sum(data.responseSad == 1 & data.responseAISad == 0 & data.changeMind == 1)/sum(data.responseSad == 1 & data.responseAISad == 0 & ~isnan(data.changeMind));
    Change_HH_HumanAIS_L3HumanAI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 1 & data.changeMind == 1)/sum(data.responseSad == 0 & data.responseAISad == 1 & ~isnan(data.changeMind));
    Change_HH_HumanAIH_L3HumanAI(iSubject, 1) = sum(data.responseSad == 0 & data.responseAISad == 0 & data.changeMind == 1)/sum(data.responseSad == 0 & data.responseAISad == 0 & ~isnan(data.changeMind));
    
    %% Percentage change mind: Agreement vs. Disagreement 
    Change_Agreement_L3HumanAI(iSubject, 1) = sum(data.responseSad == data.responseAISad & data.changeMind == 1)/sum(data.responseSad == data.responseAISad == 1);
    Change_Disagreement_L3HumanAI(iSubject, 1) = sum(data.responseSad ~= data.responseAISad & data.changeMind == 1)/sum(data.responseSad ~= data.responseAISad == 1);

    %% Total Sad (after change of mind)
    totalSadL3HumanAI(iSubject, 1) = HS_HumanAIS_L3HumanAI(iSubject, 1) + HS_HumanAIH_L3HumanAI(iSubject, 1)...
        - Change_HS_HumanAIS_L3HumanAI(iSubject, 1)*HS_HumanAIS_L3HumanAI(iSubject, 1)...
        - Change_HS_HumanAIH_L3HumanAI(iSubject, 1)*HS_HumanAIH_L3HumanAI(iSubject, 1)...
        + Change_HH_HumanAIS_L3HumanAI(iSubject, 1)*HH_HumanAIS_L3HumanAI(iSubject, 1)...
        + Change_HH_HumanAIH_L3HumanAI(iSubject, 1)*HH_HumanAIH_L3HumanAI(iSubject, 1);
      
end

%% Figure - Part 1
[p10, d10, CI10] = permutationTest(mean(responseCollaborationL3HumanAI, 2) - responseBaselineL3HumanAI, 0);
mean(responseCollaborationL3HumanAI - responseBaselineL3HumanAI)
mean(mean(responseCollaborationL3HumanAI - responseBaselineL3HumanAI))

%% Baseline bias
[p11, d11, CI11] = permutationTest(responseBaselineBlockL3HumanAI(:,1), mean(responseBaselineBlockL3HumanAI(:,2:5),2));
mean([responseBaselineBlockL3HumanAI(:,1), mean(responseBaselineBlockL3HumanAI(:,2:5),2)])

%% Colaboration
mean([HS_HumanAIS_L3HumanAI, HS_HumanAIH_L3HumanAI, HH_HumanAIS_L3HumanAI, HH_HumanAIH_L3HumanAI])
mean([Change_HS_HumanAIS_L3HumanAI, Change_HS_HumanAIH_L3HumanAI, Change_HH_HumanAIS_L3HumanAI, Change_HH_HumanAIH_L3HumanAI])

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

%% Change mind
[p12, d12, CI12] = permutationTest(Change_Disagreement_L3HumanAI, Change_Agreement_L3HumanAI);

%%%%%%%%%%%%%%%%%%
%% Human Vs. AI %%
%%%%%%%%%%%%%%%%%%

%% Change mind
[p13, d13, CI13] = permutationTest(Change_Disagreement_L3AI, Change_Disagreement_L3Human);
[p14, d14, CI14] = permutationTest(Change_Disagreement_L3AI, Change_Disagreement_L3HumanAI);
[p15, d15, CI15] = permutationTest(Change_Disagreement_L3Human, Change_Disagreement_L3HumanAI);

mean(Change_Disagreement_L3AI)
mean(Change_Disagreement_L3Human)
mean(Change_Disagreement_L3HumanAI)

%% Figures
%% Bias
figure
hold on
colorAI = [0, 0.4375, 0.75];
colorHuman = [0.9258, 0.4883, 0.1914];
colorHumanAI = [0.5, 0, 0.5];
x = 1:3;
xticks(1:3)
xticklabels([])
y = [mean(Change_Disagreement_L3AI), mean(Change_Disagreement_L3Human), mean(Change_Disagreement_L3HumanAI)];
ystd = [std(Change_Disagreement_L3AI), std(Change_Disagreement_L3Human), std(Change_Disagreement_L3HumanAI)];
errorbar(1, y(1), ystd(1)/sqrt(50),'-o','MarkerSize', 12,'LineWidth', 1.5, 'Color', colorAI, 'MarkerEdgeColor', colorAI, 'MarkerFaceColor', colorAI, 'CapSize', 0)
errorbar(3, y(2), ystd(2)/sqrt(50),'-o','MarkerSize', 12,'LineWidth', 1.5, 'Color', colorHuman, 'MarkerEdgeColor', colorHuman, 'MarkerFaceColor', colorHuman, 'CapSize', 0)
errorbar(2, y(3), ystd(3)/sqrt(50),'-o','MarkerSize', 12,'LineWidth', 1.5, 'Color', colorHumanAI, 'MarkerEdgeColor', colorHumanAI, 'MarkerFaceColor', colorHumanAI, 'CapSize', 0)
ylabel({'P(change mind)'})
xlim([.5 3.5])
ylim([0 0.45])
set(gca, 'FontSize', 16)

%% Independent decision figure
figure
hold on
% Colors
colorAI = [0, 0.4375, 0.75];
colorHuman = [0.9258, 0.4883, 0.1914];
colorHumanAI = [0.5, 0, 0.5];

% Responses
responseAI = responseCollaborationL3AI - responseBaselineL3AI;
responseHuman = responseCollaborationL3Human - responseBaselineL3Human;
% AI
x = 1:6;                     
y = mean(responseAI);
xconf = [x x(end:-1:1)] ;
plot1 = plot(x, y,'-o' ,'Color', colorAI, 'LineWidth', 2, 'MarkerSize',4  ,'MarkerFaceColor', colorAI)
yconf = [y - std(responseAI)/sqrt(size(responseAI, 1)) y(end:-1:1) + fliplr(std(responseAI)/sqrt(size(responseAI, 1)))];
f11 = fill(xconf, yconf,'red');
f11.FaceColor = colorAI;      
f11.EdgeColor = 'none';  
f11.FaceAlpha = 0.2;

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

% Properties
axis([0.75 6.25 -0.03 0.15])
ylabel({'Induced bias'})
set(gca, 'FontSize', 16)
legend([plot1, plot2, plot3],{'Human-AI'; 'Human-AI perceived as human'; 'Human-human'}, 'FontSize', 12)
xlabel('\bf Layer 3 \rm Blocks', 'FontSize', 20)

effectHumanAI = responseCollaborationL3HumanAI - responseBaselineL3HumanAI;
effectHuman = responseCollaborationL3Human - responseBaselineL3Human;
effectAI = responseCollaborationL3AI - responseBaselineL3AI;
allEffect = [effectAI; effectHumanAI; effectHuman];
allEffect = [allEffect, [ones(50,1); 2*ones(50,1); 3*ones(50, 1)]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indifference point analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[IP1, CI11, estimation1] = IPCI(allDataL1);
[IP2, CI22, estimation2] = IPCI(allDataL2AI);
[IP3, CI33, estimation3] = IPCI(allDataL3AI(~isnan(allDataL3AI.changeMind), :));
[IP4, CI44, estimation4] = IPCI(allDataL3Human(~isnan(allDataL3Human.changeMind), :));
[IP5, CI55, estimation5] = IPCI(allDataL3HumanAI(~isnan(allDataL3HumanAI.changeMind), :));

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