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

% loop through groups
for i_group = group_v
   
    % update data log
    i_log = i_log+1;
    
    % indices for current group
    dy_indx=find(data.groupID==i_group & data.sbjNUM==1);
    s1_indx=find(data.groupID==i_group & data.sbjNUM==1);
    s2_indx=find(data.groupID==i_group & data.sbjNUM==2);
    
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
    gdata.group(i_log,1)     = i_group;
    gdata.condition(i_log,1) = condition;
    gdata.age(i_log,1)       = dy_age;
    gdata.amin(i_log,1)      = amin;
    gdata.amax(i_log,1)      = amax;
    gdata.smin(i_log,1)      = smin;
    gdata.smax(i_log,1)      = smax;
    gdata.dslo(i_log,1)      = dy_slope;
    gdata.dacc(i_log,1)      = dy_acc_mean;
    gdata.sminsmax(i_log,1)  = sminsmax;
    gdata.cbaccmax(i_log,1)  = cb_acc_max;
    gdata.cbslomax(i_log,1)  = cb_slope_max;
    gdata.optimality(i_log,1) = optimality;
    gdata.disagree(i_log,1)  = mean(data.disagree);
    gdata.deliberate(i_log,1) = mean_dyad_rt;
     
end

% age group indices
cindx=gdata.condition==1;
aindx=gdata.condition==2;
gindx=gdata.condition==3;

% PLOT RESULTS
figz=figure('color',[1 1 1]);
% SIMILARITY
subplot(2,2,1);
group_mu = [mean(gdata.sminsmax(cindx)) mean(gdata.sminsmax(aindx))];
group_se = [std(gdata.sminsmax(cindx))/sqrt(sum((cindx))) std(gdata.sminsmax(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,gdata.sminsmax(cindx),'k.');
plot(2.2+rand(1,sum(aindx))./10,gdata.sminsmax(aindx),'k.');
plot([0 3],[mean(gdata.sminsmax(gindx)) mean(gdata.sminsmax(gindx))],'r-');
box('off')
xlim([0 3]); ylim([0 1]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[0:.25:1]);
title('similarity','FontWeight','normal')
ylabel('slope_m_i_n/slope_m_a_x')
set(gca,'FontSize',14)
axis square
% COLLECTIVE BENEFIT: SENSITIVITY
subplot(2,2,2);
group_mu = [mean(gdata.cbslomax(cindx)) mean(gdata.cbslomax(aindx))];
group_se = [std(gdata.cbslomax(cindx))/sqrt(sum((cindx))) std(gdata.cbslomax(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,gdata.cbslomax(cindx),'k.');
plot(2.2+rand(1,sum(aindx))./10,gdata.cbslomax(aindx),'k.');
plot([0 3],[mean(gdata.cbslomax(gindx)) mean(gdata.cbslomax(gindx))],'r-');
box('off')
xlim([0 3]); ylim([.6 1.6]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[.6:.2:1.6]);
title('collective benefit','FontWeight','normal')
ylabel('slope_d_y_a_d/slope_m_a_x');
set(gca,'FontSize',14)
axis square
% OPTIMALITY
subplot(2,2,3);
group_mu = [mean(gdata.optimality(cindx)) mean(gdata.optimality(aindx))];
group_se = [std(gdata.optimality(cindx))/sqrt(sum((cindx))) std(gdata.optimality(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,gdata.optimality(cindx),'k.');
plot(2.2+rand(1,sum(aindx))./10,gdata.optimality(aindx),'k.');
plot([0 3],[mean(gdata.optimality(gindx)) mean(gdata.optimality(gindx))],'r-');
box('off')
xlim([0 3]); ylim([.7 1.2]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[.7:.1:1.2]);
title('optimality','FontWeight','normal')
ylabel('slope_d_y_a_d/slope_W_C_S');
set(gca,'FontSize',14)
axis square
% DELIBERATION TIME
subplot(2,2,4);
group_mu = [mean(gdata.deliberate(cindx)) mean(gdata.deliberate(aindx))];
group_se = [std(gdata.deliberate(cindx))/sqrt(sum((cindx))) std(gdata.deliberate(aindx))/sqrt(sum((aindx)))];
bar(group_mu,'FaceColor',[.7 .7 .7],'LineWidth',1); hold on;
i_v=[1 2];    for i=1:length(i_v); plot([i_v(i) i_v(i)],[group_mu(i)-group_se(i) group_mu(i)+group_se(i)],'k-','LineWidth',2); end;
plot(.8+rand(1,sum(cindx))./10,gdata.deliberate(cindx),'k.');
plot(2.2+rand(1,sum(aindx))./10,gdata.deliberate(aindx),'k.');
plot([0 3],[nanmean(gdata.deliberate(gindx)) nanmean(gdata.deliberate(gindx))],'r-');
box('off')
xlim([0 3]); ylim([5 25]);
set(gca,'XTick',[1 2],'XTickLabel',{'YA','OA'})
set(gca,'YTick',[5:5:25]);
title('deliberation','FontWeight','normal')
ylabel('seconds')
set(gca,'FontSize',14)
axis square