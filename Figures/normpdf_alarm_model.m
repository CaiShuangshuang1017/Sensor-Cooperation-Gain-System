%绘制正态分布报警模型
clear;
clc;

figure;

xlabel('H(m)','FontSize',12);hold on;
ylabel('PDF','FontSize',12);hold on;

%初始水压正态分布
x1=28:0.001:31;
y1=normpdf(x1,30,0.25);
plot(x1,y1,'color',[238/258,85/258,61/258],LineWidth=1.5);hold on;
grid on;
set(gca,'GridLineStyle','--')
%爆管后压力分布
delta_H=1;
x3=28:0.001:31;
y3=normpdf(x3,30-delta_H,0.25);
plot(x3,y3,'color',[86/258,121/258,186/258],LineWidth=1.5);hold on;
%报警阈值
n=norminv(0.01,30,0.25);
xline(n,'--','Alarm Threshold (H_i^{alarm})','FontName','Times New Roman','FontSize',14,'color',[016/258,139/258,150/258],LineWidth=1.5)
%误报率
x2=29:0.001:n-0.005;
y2=normpdf(x2,30,0.25);
area(x2,y2,'FaceColor',[238/258,85/258,61/258]);hold on;
%漏报率
x4=n+0.005:0.001:31-delta_H;
y4=normpdf(x4,30-delta_H,0.25);
area(x4,y4,'FaceColor',[86/258,121/258,186/258]);hold on;

%标注
str='\DeltaH';
text(29.45,0.65,str,'interpreter','tex','Color',[86/258,121/258,186/258],'FontSize',14);
xline(30,'--','color',[238/258,85/258,61/258],LineWidth=1.5);
xline(30-delta_H,'--','color',[86/258,121/258,186/258],LineWidth=1.5);
arrow([29,0.6],[30,0.6],'color',[86/258,121/258,186/258]);
arrow([30,0.6],[29,0.6],'color',[86/258,121/258,186/258]);
legend('Normal Condition','Burst Condition','Alarm Threshold','False Alarm Rate','Missed Alarm Rate');hold on;
ax=gca;
ax.YLim=[0,1.8];
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0)
