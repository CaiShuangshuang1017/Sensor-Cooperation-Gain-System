%%%绘制自变量为误报率，因变量为阈值和漏报率
colororder([238,85,61;86,121,186]./258)

x=linspace(0.001,1);
mean=30;
std=0.25;
delta_H=1;
y=mean+std*norminv(x/100);
yyaxis left
plot(x,y,'--',LineWidth=1.5);
xlabel('False Alarm Rate(%)')
ylabel('Alarm Threshold (m)');

z=(1-normcdf(norminv(x/100)+(1/std)*delta_H))*100;
yyaxis right
plot(x,z,LineWidth=1.5)
ylim([0,60]);
ylabel('Missed Alarm Rate(%)')
grid on;
set(gca,'GridLineStyle','--');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.0);