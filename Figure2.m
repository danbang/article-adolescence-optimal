% Haller*, Bang*, Bahrami & Lau (2018) Group decision-making is optimal
% in adolescence [*equal contribution]
%
% Dan Bang danbang.db@gmail.com 2018

% fresh memory
clc;clear;close;

% add paths
addpath('helpers');

% load data
load('data.mat');

% vector unique subject IDs
sbj_v = unique(data.sbjID);

% initialise variable for logging statistics
i_log = 0;

% loop through subjects
for i_sbj = sbj_v
   
    % update data log
    i_log = i_log+1;
    
    % indices for current subject
    indx=find(data.sbjID==i_sbj);
    
    % load vector data
    stm_i_v         = data.stimInterval(indx);
    stm_k_v         = data.stimContrast(indx);
    stm_d_v         = data.stimDelta(indx);
    sbj_binary_v    = data.sbjChoice(indx);
    sbj_acc_v       = data.sbjAcc(indx);
    sbj_arbi        = data.sbjArbi(indx);
    dya_binary_v    = data.dyaChoice(indx);
    dya_disagree    = data.disagree(indx);
    
    % load scalar data
    sbj_age         = unique(data.sbjMONTHS(indx));
    sbj_group       = unique(data.groupID(indx));
    condition       = unique(data.condition(indx));
    
    % output measures
    % accuracy
    acc_mean        = mean(sbj_acc_v);
    % sensitivity
    slope           = quickSlope(stm_d_v',sbj_binary_v');
    % RT
    sbj_rt1 = data.sbjRT(indx);
    sbj_ses = data.session(indx);
    sbj_even = mod(2,2)==0; 
    rt1_mean        = nanmean(sbj_rt1(sbj_ses==sbj_even+1))./1000;
    % egocentric
    ego_mean        = sum(sbj_binary_v(sbj_arbi==1)==dya_binary_v(sbj_arbi==1))./sum(sbj_arbi);
    
    % store individual statistics (idata)
    idata.subject(i_log,1)   = i_sbj;
    idata.group(i_log,1)     = sbj_group;
    idata.condition(i_log,1) = condition;
    idata.age(i_log,1)       = sbj_age;
    idata.accuracy(i_log,1)  = acc_mean;
    idata.slope(i_log,1)     = slope;
    idata.egochoice(i_log,1) = ego_mean;
    idata.rt1(i_log,1)       = rt1_mean;
    idata.efficientacc(i_log,1)      = acc_mean/rt1_mean;
    idata.efficientslo(i_log,1)      = slope/rt1_mean;
    
end

% age group indices
cindx=idata.condition==1;
aindx=idata.condition==2;
gindx=idata.condition==3;

% PLOT RESULTS
figz=figure('color',[1 1 1]);
% ACCURACY
subplot(2,2,1);
group_mu = [mean(idata.accuracy(cindx)) mean(idata.accuracy(aindx))];
group_se = [std(idata.accuracy(cindx))/sqrt(sum((cindx))) std(idata.accuracy(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,idata.accuracy(cindx),'k.','MarkerSize',12);
plot(2.2+rand(1,sum(aindx))./10,idata.accuracy(aindx),'k.','MarkerSize',12);
box('off')
xlim([0 3]); ylim([.5 1]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[.5:.1:1]);
title('accuracy','FontWeight','normal')
ylabel('percent correct')
set(gca,'FontSize',16,'LineWidth',2)
% text(-1,1.05,'A','clipping','off','FontSize',30,'FontWeight','Bold')
% SENSITIVITY
subplot(2,2,2);
group_mu = [mean(idata.slope(cindx)) mean(idata.slope(aindx))];
group_se = [std(idata.slope(cindx))/sqrt(sum((cindx))) std(idata.slope(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,idata.slope(cindx),'k.','MarkerSize',12);
plot(2.2+rand(1,sum(aindx))./10,idata.slope(aindx),'k.','MarkerSize',12);
box('off')
xlim([0 3]); ylim([0 4]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[0:1:4]);
title('sensitivity','FontWeight','normal')
ylabel('slope')
set(gca,'FontSize',16,'LineWidth',2)
% text(-1,4.4,'B','clipping','off','FontSize',30,'FontWeight','Bold')
% EGOCENTRIC
subplot(2,2,3);
group_mu = [mean(idata.egochoice(cindx)) mean(idata.egochoice(aindx))];
group_se = [std(idata.egochoice(cindx))/sqrt(sum((cindx))) std(idata.egochoice(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
% plot([0 3],[.5 .5],'k-'); hold on
plot(.8+rand(1,sum(cindx))./10,idata.egochoice(cindx),'k.','MarkerSize',12);
plot(2.2+rand(1,sum(aindx))./10,idata.egochoice(aindx),'k.','MarkerSize',12);
box('off')
xlim([0 3]); %ylim([0 4]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[0:.2:1]);
title('egocentric bias','FontWeight','normal')
ylabel('follow own decision')
set(gca,'FontSize',16,'LineWidth',2)
% text(-1,1.05,'C','clipping','off','FontSize',30,'FontWeight','Bold')
% RT
subplot(2,2,4);
group_mu = [mean(idata.rt1(cindx)) mean(idata.rt1(aindx))];
group_se = [std(idata.rt1(cindx))/sqrt(sum((cindx))) std(idata.rt1(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,idata.rt1(cindx),'k.','MarkerSize',12);
plot(2.2+rand(1,sum(aindx))./10,idata.rt1(aindx),'k.','MarkerSize',12);
box('off')
xlim([0 3]); ylim([0 4]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[0:1:4]);
title('reaction time','FontWeight','normal')
ylabel('seconds')
set(gca,'FontSize',16,'LineWidth',2)
% text(-1,4.4,'D','clipping','off','FontSize',30,'FontWeight','Bold')
print('-djpeg','-r300',['matlab-Figure2']);