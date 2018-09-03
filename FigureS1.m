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
    
    % load scalar data
    gdata.group(i_log,1)     = i_group;
    gdata.s1_age(i_log,1)    = unique(data.sbjMONTHS(s1_indx));
    gdata.s2_age(i_log,1)    = unique(data.sbjMONTHS(s2_indx));
    gdata.mu_age(i_log,1)    = mean([unique(data.sbjMONTHS(s1_indx)) unique(data.sbjMONTHS(s2_indx))]);
    gdata.condition(i_log,1) = unique(data.condition(s1_indx));
    
end

% age group indices
cindx=gdata.condition==1;
aindx=gdata.condition==2;
gindx=gdata.condition==3;

% sort 
child_groupAge= gdata.mu_age(cindx);
child_sbj01Age= gdata.s1_age(cindx);
child_sbj02Age= gdata.s2_age(cindx);
[sortVal sortID]= sort(child_groupAge,'ascend');
for i= 1:length(child_sbj01Age)
    child_sbj01AgeS(i)= child_sbj01Age(sortID(i));
    child_sbj02AgeS(i)= child_sbj02Age(sortID(i));
end

% sort 
adols_groupAge= gdata.mu_age(aindx);
adols_sbj01Age= gdata.s1_age(aindx);
adols_sbj02Age= gdata.s2_age(aindx);
[sortVal sortID]= sort(adols_groupAge,'ascend');
for i= 1:length(adols_sbj01Age)
    adols_sbj01AgeS(i)= adols_sbj01Age(sortID(i));
    adols_sbj02AgeS(i)= adols_sbj02Age(sortID(i));
end

%% INDIVIDUAL PERFORMANCE
figz=figure('color',[1 1 1]);
hold on;
plot(1:length(child_sbj01AgeS),child_sbj01AgeS,'ro','MarkerSize',10);
plot(1:length(child_sbj02AgeS),child_sbj02AgeS,'r^','MarkerSize',10);
plot(length(child_sbj01AgeS)+2:length(child_sbj01AgeS)+length(adols_sbj01AgeS)+1,adols_sbj01AgeS,'bo','MarkerSize',10);
plot(length(child_sbj01AgeS)+2:length(child_sbj01AgeS)+length(adols_sbj01AgeS)+1,adols_sbj02AgeS,'b^','MarkerSize',10);
plot([19 19],[7 19],'k--','LineWidth',2);
xlabel('group');
ylabel('age');
set(gca,'YTick',8:2:18);
set(gca,'XTick',[]);
set(gca,'FontSize',30,'LineWidth',2)
ylim([7 19]);
xlim([0 36]);
legend('YA: subject 1', 'YA: subject 2', 'OA: subject 1', 'OA: subject 2','Location','NorthWest');
legend boxoff
print('-djpeg','-r300',['matlab-FigureS1']);