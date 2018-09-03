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

% vector unique group IDs
group_v = unique(data.groupID);

% initialise variable for logging statistics
i_log = 0;

% initialise output matrix
data_groups=[];

% loop through groups
for i_group = group_v
   
    % update data log
    i_log = i_log+1;
    
    % loop through sessions
    for i_session = 1:2
    
    % indices for current group and session
    dy_indx=find(data.groupID==i_group & data.sbjNUM==1 & data.session==i_session);
    s1_indx=find(data.groupID==i_group & data.sbjNUM==1 & data.session==i_session);
    s2_indx=find(data.groupID==i_group & data.sbjNUM==2 & data.session==i_session);
    
    % load vector data
    stm_i_v         = data.stimInterval(dy_indx);
    stm_k_v         = data.stimContrast(dy_indx);
    stm_d_v         = data.stimDelta(dy_indx); 
    s1_binary_v     = data.sbjChoice(s1_indx);
    s1_acc_v        = data.sbjAcc(s1_indx);
    s1_arbi         = data.sbjArbi(s1_indx);    
    s2_binary_v     = data.sbjChoice(s2_indx);
    s2_acc_v        = data.sbjAcc(s2_indx);
    s2_arbi         = data.sbjArbi(s2_indx);   
    dy_binary_v    = data.dyaChoice(dy_indx);
    dy_acc_v       = data.dyaAcc(dy_indx);
    dy_disagree    = data.disagree(dy_indx);
    
    % load scalar data
    s1_age          = unique(data.sbjMONTHS(s1_indx));
    s2_age          = unique(data.sbjMONTHS(s2_indx));
    dy_age          = mean([s1_age s2_age]);
    condition       = unique(data.condition(s1_indx));
    
    % output measures
    % accuracy
    s1_acc_mean     = mean(s1_acc_v);
    s2_acc_mean     = mean(s2_acc_v);
    dy_acc_mean     = mean(dy_acc_v);
    % sensitivity
    s1_slope        = quickSlope(stm_d_v',s1_binary_v');
    s2_slope        = quickSlope(stm_d_v',s2_binary_v');
    dy_slope        = quickSlope(stm_d_v',dy_binary_v');
    % similarity
    amin            = min([s1_acc_mean s2_acc_mean]);
    amax            = max([s1_acc_mean s2_acc_mean]);
    smin            = min([s1_slope s2_slope]);
    smax            = max([s1_slope s2_slope]);
    amean           = mean([s1_acc_mean s2_acc_mean]);
    smean           = mean([s1_slope s2_slope]);
    aminamax        = min([s1_acc_mean s2_acc_mean])/max([s1_acc_mean s2_acc_mean]);
    sminsmax        = min([s1_slope s2_slope])/max([s1_slope s2_slope]);
    % collective benefit: accuracy
    cb_acc_min      = dy_acc_mean/min([s1_acc_mean s2_acc_mean]);
    cb_acc_max      = dy_acc_mean/max([s1_acc_mean s2_acc_mean]);
    cb_acc_mean     = dy_acc_mean/mean([s1_acc_mean s2_acc_mean]);
    % collective benefit: sensitivity
    cb_slope_min    = dy_slope/min([s1_slope s2_slope]);
    cb_slope_max    = dy_slope/max([s1_slope s2_slope]);
    cb_slope_mean   = dy_slope/mean([s1_slope s2_slope]);
    % optimality
    optimality      = dy_slope/((s1_slope+s2_slope)/(2^.5));
    % deliberation time
    disagree = data.disagree(dy_indx);
    dyadcort = data.dyaRT(dy_indx);
    mean_dyad_rt = nanmean(dyadcort(disagree==1))./1000;
    
    % store group statistics (gdata)
    gdata.group(i_log,i_session)      = i_group;
    gdata.condition(i_log,i_session)  = condition;
    gdata.age(i_log,i_session)        = dy_age;
    gdata.amin(i_log,i_session)       = amin;
    gdata.amax(i_log,i_session)       = amax;
    gdata.smin(i_log,i_session)       = smin;
    gdata.smax(i_log,i_session)       = smax;
    gdata.dslo(i_log,i_session)       = dy_slope;
    gdata.sminsmax(i_log,i_session)   = sminsmax;
    gdata.dacc(i_log,i_session)       = dy_acc_mean;
    gdata.cbaccmax(i_log,i_session)   = cb_acc_max;
    gdata.cbslomax(i_log,i_session)   = cb_slope_max;
    gdata.optimality(i_log,i_session) = optimality;
    gdata.disagree(i_log,i_session)   = mean(data.disagree);
    gdata.deliberate(i_log,i_session) = mean_dyad_rt;
    
    end
    
    
end

% age group indices
cindx=gdata.condition(:,1)==1;
aindx=gdata.condition(:,1)==2;
gindx=gdata.condition(:,1)==3;

% PLOT RESULTS
figz=figure('color',[1 1 1]);
% SIMILARITY
subplot(2,2,1);
group_mu = [mean(gdata.sminsmax(cindx,:)); mean(gdata.sminsmax(aindx,:))];
group_se = [std(gdata.sminsmax(cindx,:))./sqrt(sum((cindx))); std(gdata.sminsmax(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([.50 .80]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[.50:.10:.80]);
title('similarity','FontWeight','normal')
ylabel('slope_m_i_n/slope_m_a_x')
set(gca,'FontSize',16,'LineWidth',2) 

% COLLECTIVE BENEFIT: SENSITIVITY
subplot(2,2,2);
group_mu = [mean(gdata.cbslomax(cindx,:)); mean(gdata.cbslomax(aindx,:))];
group_se = [std(gdata.cbslomax(cindx,:))./sqrt(sum((cindx))); std(gdata.cbslomax(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([.9 1.4]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[.9:.1:1.4]);
title('collective benefit','FontWeight','normal')
ylabel('slope_d_y_a_d/slope_m_a_x');
set(gca,'FontSize',16,'LineWidth',2) 
% OPTIMALITY
subplot(2,2,3);
group_mu = [mean(gdata.optimality(cindx,:)); mean(gdata.optimality(aindx,:))];
group_se = [std(gdata.optimality(cindx,:))./sqrt(sum((cindx))); std(gdata.cbslomax(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([.85 1.15]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[.85:.05:1.15]);
title('optimality','FontWeight','normal')
ylabel('slope_d_y_a_d/slope_W_C_S');
set(gca,'FontSize',16,'LineWidth',2) 
% DELIBERATION TIME
subplot(2,2,4);
group_mu = [mean(gdata.deliberate(cindx,:)); mean(gdata.deliberate(aindx,:))];
group_se = [std(gdata.deliberate(cindx,:))./sqrt(sum((cindx))); std(gdata.deliberate(aindx,:))./sqrt(sum((aindx)))];
plot([1 2],group_mu(1,:),'k-','LineWidth',2); hold on;
plot([1 2],group_mu(2,:),'k--','LineWidth',2); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(1,i)-group_se(1,i) group_mu(1,i)+group_se(1,i)],'k-','LineWidth',1); end;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(2,i)-group_se(2,i) group_mu(2,i)+group_se(2,i)],'k-','LineWidth',1); end;
plot([1 2],group_mu(1,:),'ko','MarkerFaceColor',[1 1 1],'LineWidth',1); hold on;
plot([1 2],group_mu(2,:),'ko','MarkerFaceColor',[0 0 0],'LineWidth',1); hold on;
box('off')
xlim([.5 2.5]); ylim([8 16]);
set(gca,'XTick',[1 2],'XTickLabel',{'session 1','session 2'})
set(gca,'YTick',[8:2:16]);
title('deliberation','FontWeight','normal')
ylabel('seconds');
set(gca,'FontSize',16,'LineWidth',2) 
print('-djpeg','-r300',['matlab-Figure5']);