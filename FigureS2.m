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
sindx=(idata.condition==1|idata.condition==2);

% PLOT RESULTS
figz=figure('color',[1 1 1]);
% ACCURACY
subplot(2,2,1);
x  = zscore(idata.age(sindx));
y  = zscore(idata.accuracy(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot([-10 10],[0 0],'k--'); hold on;
plot([0 0],[-10 10],'k--'); hold on;
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
title('accuracy','FontWeight','normal')
ylabel('percent correct')
xlabel('age')
set(gca,'FontSize',16,'LineWidth',2) 
% SENSITIVITY
subplot(2,2,2);
x  = zscore(idata.age(sindx));
y  = zscore(idata.slope(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot([-10 10],[0 0],'k--'); hold on;
plot([0 0],[-10 10],'k--'); hold on;
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
title('sensitivity','FontWeight','normal')
ylabel('slope')
xlabel('age')
set(gca,'FontSize',16,'LineWidth',2) 
% EGOCENTRIC
subplot(2,2,3);
x  = zscore(idata.age(sindx));
y  = zscore(idata.egochoice(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot([-10 10],[0 0],'k--'); hold on;
plot([0 0],[-10 10],'k--'); hold on;
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
title('egocentric bias','FontWeight','normal')
ylabel('follow own decision')
xlabel('age')
set(gca,'FontSize',16,'LineWidth',2) 
% REACTION TIME
subplot(2,2,4);
x  = zscore(idata.age(sindx));
y  = zscore(idata.rt1(sindx));
p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
z1 = polyval(p1,[min(x):.1:max(x)]);
z2 = polyval(p2,[min(x):.1:max(x)]);
[beta1,~,stats1]=glmfit([x x.^2],y);
plot([10 11],[10 11],'r-','LineWidth',2); hold on;
plot([10 11],[10 11],'b-','LineWidth',2); hold on;
plot([-10 10],[0 0],'k--'); hold on;
plot([0 0],[-10 10],'k--'); hold on;
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
title('reaction time','FontWeight','normal')
ylabel('seconds')
xlabel('age')
set(gca,'FontSize',16,'LineWidth',2);
print('-djpeg','-r300',['matlab-FigureS2']);