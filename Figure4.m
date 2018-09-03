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
    
    % loop through sessions
    for i_session = 1:2
        
    % indices for current subject and session
    indx=find(data.sbjID==i_sbj & data.session==i_session);
    
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
    % egocentric
    ego_mean        = sum(sbj_binary_v(sbj_arbi==1)==dya_binary_v(sbj_arbi==1))./sum(sbj_arbi);
    
    % store individual statistics (idata)
    idata.subject(i_log,i_session)   = i_sbj;
    idata.group(i_log,i_session)     = sbj_group;
    idata.condition(i_log,i_session) = condition;
    idata.age(i_log,i_session)       = sbj_age;
    idata.accuracy(i_log,i_session)  = acc_mean;
    idata.slope(i_log,i_session)     = slope;
    idata.egochoice(i_log,i_session) = ego_mean;
    
    end

end

% age group indices
cindx=idata.condition(:,1)==1;
aindx=idata.condition(:,1)==2;
gindx=idata.condition(:,1)==3;

% PLOT RESULTS
figz=figure('color',[1 1 1]);
% ACCURACY
subplot(2,2,1);
group_mu = [mean(idata.accuracy(cindx,:)); mean(idata.accuracy(aindx,:))];
group_se = [std(idata.accuracy(cindx,:))./sqrt(sum((cindx))); std(idata.accuracy(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([.65 .85]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[.65:.05:.85]);
title('accuracy','FontWeight','normal')
ylabel('percent correct')
set(gca,'FontSize',16,'LineWidth',2) 
% SENSITIVITY
subplot(2,2,2);
group_mu = [mean(idata.slope(cindx,:)); mean(idata.slope(aindx,:))];
group_se = [std(idata.slope(cindx,:))./sqrt(sum((cindx))); std(idata.slope(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([1 2.5]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[1:.5:2.5]);
title('sensitivity','FontWeight','normal')
ylabel('slope')
set(gca,'FontSize',16,'LineWidth',2) 
% EGOCENTRIC
subplot(2,2,3);
group_mu = [mean(idata.egochoice(cindx,:)); mean(idata.egochoice(aindx,:))];
group_se = [std(idata.egochoice(cindx,:))./sqrt(sum((cindx))); std(idata.egochoice(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([.4 .6]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[.4:.05:.6]);
title('egocentric bias','FontWeight','normal')
ylabel('follow own decision')
set(gca,'FontSize',16,'LineWidth',2) 
print('-djpeg','-r300',['matlab-Figure4']);