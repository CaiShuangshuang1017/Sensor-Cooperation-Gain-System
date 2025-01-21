%显示监测点位置
clear;
clc;

start_toolkit;
wds = epanet('L-TOWN_sensor_position.inp');

sensor_id=[47,95,139,344,8,75,257,251];
for i=1:length(sensor_id)
    wds.setNodeBaseDemands(sensor_id(i),100);
end

wds.saveInputFile('L-TOWN_sensor_position.inp');