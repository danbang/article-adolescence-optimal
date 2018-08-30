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

sindx=(gdata.condition==1|gdata.condition==2);

% PLOT RESULTS
figz=figure('color',[1 1 1]);
% SIMILARITY
subplot(2,2,1);
x  = zscore(gdata.age(sindx));
y  = zscore(gdata.sminsmax(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot(x,y,'k.','MarkerSize',14); hold on;
plot([min(x):.1:max(x)],z1,'-','color','r','LineWidth',2); hold on; 
plot([min(x):.1:max(x)],z2,'-','color','b','LineWidth',2); hold on; 
if stats1.p(2)<.001; pval1='p<0.001'; 
else pval1=['p=',num2str(round(stats1.p(2)*1000)/1000)]; 
end;
if stats1.p(3)<.001; pval2='p<0.001'; 
else pval2=['p=',num2str(round(stats1.p(3)*1000)/1000)]; 
end;
legz=legend(pval1,pval2,'Location','NorthEast');
set(legz,'FontSize',10)
legend('boxoff');
xlim([-4 +4]); 
ylim([-4 +4]);
set(gca,'XTick',[-3:1.5:+3]);
set(gca,'YTick',[-3:1.5:+3]);
title('similarity','FontWeight','normal')
ylabel('slope_m_i_n/slope_m_a_x')
xlabel('age')
set(gca,'FontSize',14)
axis square
% COLLECTIVE BENEFIT: SENSITIVITY
subplot(2,2,2);
x  = zscore(gdata.age(sindx));
y  = zscore(gdata.cbslomax(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot(x,y,'k.','MarkerSize',14); hold on;
plot([min(x):.1:max(x)],z1,'-','color','r','LineWidth',2); hold on; 
plot([min(x):.1:max(x)],z2,'-','color','b','LineWidth',2); hold on; 
if stats1.p(2)<.001; pval1='p<0.001'; 
else pval1=['p=',num2str(round(stats1.p(2)*1000)/1000)]; 
end;
if stats1.p(3)<.001; pval2='p<0.001'; 
else pval2=['p=',num2str(round(stats1.p(3)*1000)/1000)]; 
end;
legz=legend(pval1,pval2,'Location','NorthEast');
set(legz,'FontSize',10)
legend('boxoff');
xlim([-4 +4]); 
ylim([-4 +4]);
set(gca,'XTick',[-3:1.5:+3]);
set(gca,'YTick',[-3:1.5:+3]);
title('collective benefit','FontWeight','normal')
ylabel('slope_d_y_a_d/slope_m_a_x');
xlabel('age')
set(gca,'FontSize',14)
axis square
% OPTIMALITY
subplot(2,2,3);
x  = zscore(gdata.age(sindx));
y  = zscore(gdata.optimality(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot(x,y,'k.','MarkerSize',14); hold on;
plot([min(x):.1:max(x)],z1,'-','color','r','LineWidth',2); hold on; 
plot([min(x):.1:max(x)],z2,'-','color','b','LineWidth',2); hold on; 
if stats1.p(2)<.001; pval1='p<0.001'; 
else pval1=['p=',num2str(round(stats1.p(2)*1000)/1000)]; 
end;
if stats1.p(3)<.001; pval2='p<0.001'; 
else pval2=['p=',num2str(round(stats1.p(3)*1000)/1000)]; 
end;
legz=legend(pval1,pval2,'Location','NorthEast');
set(legz,'FontSize',10)
legend('boxoff');
xlim([-4 +4]); 
ylim([-4 +4]);
set(gca,'XTick',[-3:1.5:+3]);
set(gca,'YTick',[-3:1.5:+3]);
title('optimality','FontWeight','normal')
ylabel('slope_d_y_a_d/slope_W_C_S');
xlabel('age')
set(gca,'FontSize',14)
axis square
% DELIBERATION TIME
subplot(2,2,4);
x  = zscore(gdata.age(sindx));
y  = zscore(gdata.deliberate(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot(x,y,'k.','MarkerSize',14); hold on;
plot([min(x):.1:max(x)],z1,'-','color','r','LineWidth',2); hold on; 
plot([min(x):.1:max(x)],z2,'-','color','b','LineWidth',2); hold on; 
if stats1.p(2)<.001; pval1='p<0.001'; 
else pval1=['p=',num2str(round(stats1.p(2)*1000)/1000)]; 
end;
if stats1.p(3)<.001; pval2='p<0.001'; 
else pval2=['p=',num2str(round(stats1.p(3)*1000)/1000)]; 
end;
legz=legend(pval1,pval2,'Location','NorthEast');
set(legz,'FontSize',10)
legend('boxoff');
xlim([-4 +4]); 
ylim([-4 +4]);
set(gca,'XTick',[-3:1.5:+3]);
set(gca,'YTick',[-3:1.5:+3]);
title('deliberation','FontWeight','normal')
ylabel('seconds')
xlabel('age')
set(gca,'FontSize',14)
axis square