%用不同颜色标出物理传感器的覆盖范围、虚拟传感器的覆盖扩展范围

%单个传感器监测覆盖范围
clear;
clc;


start_toolkit;
wds = epanet('L-TOWN_extended_coverage_new.inp');

%%%爆管压力下降值
burst_pipe_H_decline=xlsread("C:\Users\蔡爽爽\Desktop\L_town_Case_Study\增益系统爆管压降.xlsx",'sheet1','H2:AIC14');

%%%物理传感器覆盖范围
single_physical_sensor_coverage_matrix=xlsread("C:\Users\蔡爽爽\Desktop\L_town_Case_Study\不同测量噪声模拟\单个监测点0.3Qmax爆管压力下降值.xlsx",'sheet1','B23:AHV30');

all_pipe_cover_matrix=zeros(13,905);%%%
for i=2:906
    for j=1:13%%%
        if burst_pipe_H_decline(j,i)<burst_pipe_H_decline(j,1)
            all_pipe_cover_matrix(j,i-1)=1;
        end
    end
end

virtual_sensor_coverage_matrix=[];
virtual_sensor_all_cover_num=0;
for i=1:905
    current_pipe_cover_num=0;
    for j=1:13%%%
        if all_pipe_cover_matrix(j,i)==1
            current_pipe_cover_num=current_pipe_cover_num+1;
        end
    end
    if current_pipe_cover_num>0
        virtual_sensor_coverage_matrix(1,i)=1;
        virtual_sensor_all_cover_num=virtual_sensor_all_cover_num+1;
    end
end

physical_sensor_coverage_matrix=[];
physical_sensor_all_cover_num=0;
for i=1:905
    current_pipe_cover_num=0;
    for j=1:8
        if single_physical_sensor_coverage_matrix(j,i)==1
            current_pipe_cover_num=current_pipe_cover_num+1;
        end
    end
    if current_pipe_cover_num>0
        physical_sensor_coverage_matrix(1,i)=1;
        physical_sensor_all_cover_num=physical_sensor_all_cover_num+1;
    end
end
extend_num=0;
lose_num=0;
all_num=0;
all_lose_num=0;
for i=1:905
    if physical_sensor_coverage_matrix(1,i)==0 && virtual_sensor_coverage_matrix(1,i)==1
        disp(i);
        wds.setLinkWallReactionCoeff(i,100);
        extend_num=extend_num+1;
    end
    if physical_sensor_coverage_matrix(1,i)==1 && virtual_sensor_coverage_matrix(1,i)==0
        disp(i);
        wds.setLinkWallReactionCoeff(i,0);
        lose_num=lose_num+1;
    end
    if physical_sensor_coverage_matrix(1,i)==1 && virtual_sensor_coverage_matrix(1,i)==1
        disp(i);
        wds.setLinkWallReactionCoeff(i,50);
        all_num=all_num+1;
    end
    if physical_sensor_coverage_matrix(1,i)==0 && virtual_sensor_coverage_matrix(1,i)==0
        disp(i);
        wds.setLinkWallReactionCoeff(i,0);
        all_lose_num=all_lose_num+1;
    end
end

wds.saveInputFile('L-TOWN_extended_coverage_new.inp');
