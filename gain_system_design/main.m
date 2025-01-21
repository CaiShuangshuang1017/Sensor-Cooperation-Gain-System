clc;
clear;

tic;
pipe_virtual_best=[];
pipe_fitness_best=[];
physical_burst_H_decline=xlsread("C:\Users\蔡爽爽\Desktop\遗传算法设计增益系统\爆管后各监测点压力下降值.xlsx",'sheet1','A2:AHV9');
physical_normal_H_measurement=xlsread("C:\Users\蔡爽爽\Desktop\遗传算法设计增益系统\8个监测点报警阈值.xlsx",'sheet1','A2:ALS9');
for current_burst_pipe_index=1:905
    disp(current_burst_pipe_index);
    %% 基础参数
    virtual_num = 100;  %种群内个体数目
    physical_num = 8; %染色体节点数，也就是每个个体有多少条染色体，其实说白了就是看适应函数里有几个自变量。
    iter = 200; %迭代次数，也就是一共有多少代
    mut = 0.2;  %突变概率
    acr = 0.2; %交叉概率
    best = 1;
    virtual_matrix= zeros(virtual_num, physical_num);%存放染色体的矩阵
    fitness = zeros(virtual_num, 1);%存放染色体的适应度
    current_virtual_best = zeros(1, physical_num+1);%存放当前代的最优染色体与适应度
    fitness_ave = zeros(1, iter);%存放每一代的平均适应度
    fitness_best = zeros(1, iter);%存放每一代的最优适应度
    

    %% 初始化，这只是用于生成第一代个体，并计算其适应度函数
    virtual_matrix= Initialize(virtual_num, physical_num); %初始化染色体
    %求物理传感器关于管段的最佳灵敏度
    physical_sensitivity_matrix=-physical_burst_H_decline(:,current_burst_pipe_index+1)./physical_normal_H_measurement(:,3);
    [best_physical_sensitivity,best_physical_sensor_index]=max(physical_sensitivity_matrix);

    [fitness,virtual_matrix]=CalFitness(virtual_num, physical_num,virtual_matrix, current_burst_pipe_index,physical_burst_H_decline, physical_normal_H_measurement,best_physical_sensor_index); %计算适应度
    current_virtual_best =FindBest(virtual_matrix, fitness, physical_num); %寻找最优染色体
    fitness_best(1) = current_virtual_best(end); %将当前最优存入矩阵当中
    fitness_ave(1) = CalAveFitness(fitness); %将当前平均适应度存入矩阵当中

    %% 用于生成以下其余各代，一共迭代多少步就一共有多少代
    for t = 2:iter
        virtual_matrix = MutVirtualMatrix(virtual_matrix, mut, virtual_num, physical_num,best_physical_sensor_index,current_burst_pipe_index,physical_normal_H_measurement,physical_burst_H_decline); %变异
        virtual_matrix = AcrVirtualMatrix(virtual_matrix, acr, virtual_num, physical_num,best_physical_sensor_index,current_burst_pipe_index,physical_normal_H_measurement,physical_burst_H_decline); %交叉
        [fitness,virtual_matrix]= CalFitness(virtual_num, physical_num,virtual_matrix, current_burst_pipe_index,physical_burst_H_decline, physical_normal_H_measurement,best_physical_sensor_index); %计算适应度
        virtual_best_temp = FindBest(virtual_matrix, fitness, physical_num); %寻找最优染色体
        if virtual_best_temp(end)>current_virtual_best(end) %替换掉当前储存的最优
            current_virtual_best = virtual_best_temp;
        elseif virtual_best_temp(end)==current_virtual_best(end) && sum(virtual_best_temp(1:end-1),'all')<sum(current_virtual_best(1:end-1),'all')
            current_virtual_best = virtual_best_temp;
        end
        %%替换掉最劣
        [virtual_matrix, fitness] = ReplaceWorse(virtual_matrix, current_virtual_best, fitness);
        fitness_best(t) = current_virtual_best(end); %将当前最优存入矩阵当中
        fitness_ave(t) = CalAveFitness(fitness); %将当前平均适应度存入矩阵当中
    end

    % %% 作图
    % figure(1)
    % plot(1:iter, fitness_ave, 'r', 1:iter, fitness_best, 'b')
    % grid on
    % legend('平均适应度', '最优适应度')

    %% 输出结果
    disp(['最优染色体为', num2str(current_virtual_best(1:end-1))])
    disp(['最优适应度为', num2str(current_virtual_best(end))])
    pipe_virtual_best=[pipe_virtual_best;current_virtual_best(1:end-1)];
    pipe_fitness_best=[pipe_fitness_best;current_virtual_best(end)];
    toc;
end

pipe_virtual_list=[];
for pipe_index=1:size(pipe_virtual_best,1)
    current_pipe_virtual='';
    for sensor_index=1:8
        if pipe_virtual_best(pipe_index,sensor_index)==1
            current_pipe_virtual=strcat(current_pipe_virtual,num2str(sensor_index));
        end
    end
    current_pipe_virtual=str2num(current_pipe_virtual);
    pipe_virtual_list=[pipe_virtual_list;current_pipe_virtual];
end
results=tabulate(string(pipe_virtual_list));


%% 初始化种群
function virtual_matrix= Initialize(N, physical_num)
virtual_matrix= round(rand(N, physical_num));
end

%% 计算适应度
function [fitness,new_virtual_matrix] =CalFitness(virtual_num, physical_num,virtual_matrix, current_burst_pipe_index,physical_burst_H_decline, physical_normal_H_measurement,best_physical_sensor_index)

fitness = zeros(virtual_num,1);

%开始计算适应度
for current_virtual = 1:virtual_num
    void=0;
    %验证是否为空
    if sum(virtual_matrix(current_virtual,:),'all')==0
        void=1;
    else
        %%%%%条件1验证%%%%%%%%%%%%
        %计算虚拟传感器标准差
        current_virtual_normal_measurement=[];
        for virtual_normal_measurement_num=1:1000
            current_virtual_normal_measurement(virtual_normal_measurement_num)=sum(virtual_matrix(current_virtual,:)*physical_normal_H_measurement(:,7+virtual_normal_measurement_num),'all');
        end
        current_virtual_std=std(current_virtual_normal_measurement);
        %验证条件1
        sum_physical_std=sum(virtual_matrix(current_virtual,:)*physical_normal_H_measurement(:,3),'all');
        if current_virtual_std>sum_physical_std && abs(current_virtual_std-sum_physical_std)>0.0000001
            void=1;
        else
            %%%%%条件2验证%%%%%%%%%%%%
            %求虚拟传感器爆管压降值
            current_virtual_H_decline=sum(virtual_matrix(current_virtual,:)*physical_burst_H_decline(:,1+current_burst_pipe_index),'all');
            %验证条件2
            best_physical_std=physical_normal_H_measurement(best_physical_sensor_index,3);
            best_physical_H_decline=physical_burst_H_decline(best_physical_sensor_index,1+current_burst_pipe_index);
            H_decline_ratio=(-current_virtual_H_decline)/(-best_physical_H_decline);
            if H_decline_ratio<(current_virtual_std/best_physical_std) && abs(H_decline_ratio-(current_virtual_std/best_physical_std))>0.0000001  
                void=1;
            else
                if H_decline_ratio>(sum_physical_std/best_physical_std) && abs(H_decline_ratio-(sum_physical_std/best_physical_std))>0.0000001
                    void=1;
                end
            end
        end
    end
    
    
    %条件不通过重新生成新虚拟传感器
    while void==1
        %重新生成虚拟传感器
        new_virtual=round(rand(1, physical_num));

        if sum(new_virtual(1,:),'all')==0
            void=1;
        else
            %计算虚拟传感器标准差
            new_virtual_normal_measurement=[];
            for virtual_normal_measurement_num=1:1000
                new_virtual_normal_measurement(virtual_normal_measurement_num)=sum(new_virtual(1,:)*physical_normal_H_measurement(:,7+virtual_normal_measurement_num),'all');
            end
            new_virtual_std=std(new_virtual_normal_measurement);
            new_sum_physical_std=sum(new_virtual(1,:)*physical_normal_H_measurement(:,3),'all');
            if new_virtual_std>new_sum_physical_std  && abs(new_virtual_std-new_sum_physical_std)>0.0000001
                void=1;
            else
                %求虚拟传感器爆管压降值
                new_virtual_H_decline=sum(new_virtual(1,:)*physical_burst_H_decline(:,1+current_burst_pipe_index),'all');
                best_physical_std=physical_normal_H_measurement(best_physical_sensor_index,3);
                best_physical_H_decline=physical_burst_H_decline(best_physical_sensor_index,1+current_burst_pipe_index);
                new_H_decline_ratio=(-new_virtual_H_decline)/(-best_physical_H_decline);
                if new_H_decline_ratio<(new_virtual_std/best_physical_std) && abs(new_H_decline_ratio-(new_virtual_std/best_physical_std))>0.0000001
                    void=1;
                else
                    if new_H_decline_ratio>(new_sum_physical_std/best_physical_std) && abs(new_H_decline_ratio-(new_sum_physical_std/best_physical_std))>0.0000001
                        void=1;
                    else
                        void=0;
                        %替换原虚拟传感器
                        virtual_matrix(current_virtual,:)=new_virtual(1,:);
                        current_virtual_H_decline=new_virtual_H_decline;
                        current_virtual_std=new_virtual_std;
                    end
                end
            end
        end
    end
    %计算适应度
    fitness(current_virtual)=(-current_virtual_H_decline/current_virtual_std)-(-best_physical_H_decline/best_physical_std);
end
new_virtual_matrix=virtual_matrix;
end

%% 寻找当前最优虚拟传感器及对应适应度
function virtual_best = FindBest(virtual_matrix, fitness, physical_num)
virtual_best = zeros(1, physical_num+1);
[max_fitness, max_virtual_index] = max(fitness);%因为所有个体对应的适应度大小都被存放在fitness矩阵中
virtual_best(1:physical_num) =virtual_matrix(max_virtual_index, :);
virtual_best(end) = max_fitness;
end

%% 变异处理
function virtual_matrix_new = MutVirtualMatrix(virtual_matrix, mut, virtual_num, physical_num, best_physical_sensor_index,current_burst_pipe_index,physical_normal_H_measurement,physical_burst_H_decline)
for i = 1:virtual_num %%N是个体总数，也就是每一代有多少头袋鼠
    original_virtual=virtual_matrix(i,:);
    for j = 1:physical_num  %N_chrom是染色体节点数，就是有几条染色体
        mut_rand = rand; %随机生成一个数，代表自然里的基因突变，然后用改值来决定是否产生突变。
        if mut_rand <=mut%*(iter-t)/iter  %mut代表突变概率，即产生突变的阈值，如果小于0.2的基因突变概率阈值才进行基因突变处理，否者不进行突变处理
            virtual_matrix(i, j)=1-virtual_matrix(i, j);
        end
    end

    %检验是否越界
    void=0;
    if sum(virtual_matrix(i,:),'all')==0
        void=1;
    else
        %%%%%条件1验证%%%%%%%%%%%%
        %计算虚拟传感器标准差
        mut_virtual_normal_measurement=[];
        for virtual_normal_measurement_num=1:1000
            mut_virtual_normal_measurement(virtual_normal_measurement_num)=sum(virtual_matrix(i,:)*physical_normal_H_measurement(:,7+virtual_normal_measurement_num),'all');
        end
        mut_virtual_std=std(mut_virtual_normal_measurement);
        %验证条件1
        mut_sum_physical_std=sum(virtual_matrix(i,:)*physical_normal_H_measurement(:,3),'all');
        if mut_virtual_std>mut_sum_physical_std && abs(mut_virtual_std-mut_sum_physical_std)>0.0000001
            void=1;
        else
            %%%%%条件2验证%%%%%%%%%%%%
            %求虚拟传感器爆管压降值
            mut_virtual_H_decline=sum(virtual_matrix(i,:)*physical_burst_H_decline(:,1+current_burst_pipe_index),'all');
            %验证条件2
            best_physical_std=physical_normal_H_measurement(best_physical_sensor_index,3);
            best_physical_H_decline=physical_burst_H_decline(best_physical_sensor_index,1+current_burst_pipe_index);
            H_decline_ratio=mut_virtual_H_decline/best_physical_H_decline;
            if H_decline_ratio<(mut_virtual_std/best_physical_std) && abs(H_decline_ratio-(mut_virtual_std/best_physical_std))>0.0000001
                void=1;
            else
                if H_decline_ratio>(mut_sum_physical_std/best_physical_std) && abs(H_decline_ratio-(mut_sum_physical_std/best_physical_std))>0.0000001
                    void=1;
                end
            end
        end
    end
 
    if void==1
        virtual_matrix(i,:)=original_virtual;
    end
end
virtual_matrix_new = virtual_matrix;%%把变异处理完后的结果存在新矩阵里
end

%% 交叉处理
function virtual_matrix_new = AcrVirtualMatrix(virtual_matrix, acr, virtual_num, physical_num,best_physical_sensor_index,current_burst_pipe_index,physical_normal_H_measurement,physical_burst_H_decline)
for i= 1:virtual_num
    original_virtual=virtual_matrix(i,:);
    acr_rand = rand;%生成一个代表该个体是否产生交叉的概率大小，用于判别是否进行交叉处理
    if acr_rand<acr %如果该个体的交叉概率值大于产生交叉处理的阈值，则对该个体的染色体（两条，因为此案例中有两个自变量）进行交叉处理
        acr_chrom = floor((virtual_num-1)*rand+1); %要交叉的染色体
        acr_node = floor((physical_num-1)*rand+1); %要交叉的节点
        %交叉开始
        temp = virtual_matrix(i, acr_node);
        virtual_matrix(i, acr_node) = virtual_matrix(acr_chrom, acr_node); 
        virtual_matrix(acr_chrom, acr_node) = temp;
    end

    %检验是否越界
    void=0;
    if sum(virtual_matrix(i,:),'all')==0
        void=1;
    else
        %%%%%条件1验证%%%%%%%%%%%%
        %计算虚拟传感器标准差
        mut_virtual_normal_measurement=[];
        for virtual_normal_measurement_num=1:1000
            mut_virtual_normal_measurement(virtual_normal_measurement_num)=sum(virtual_matrix(i,:)*physical_normal_H_measurement(:,7+virtual_normal_measurement_num),'all');
        end
        mut_virtual_std=std(mut_virtual_normal_measurement);
        %验证条件1
        mut_sum_physical_std=sum(virtual_matrix(i,:)*physical_normal_H_measurement(:,3),'all');
        if mut_virtual_std>mut_sum_physical_std && abs(mut_virtual_std-mut_sum_physical_std)>0.0000001
            void=1;
        else
            %%%%%条件2验证%%%%%%%%%%%%
            %求虚拟传感器爆管压降值
            mut_virtual_H_decline=sum(virtual_matrix(i,:)*physical_burst_H_decline(:,1+current_burst_pipe_index),'all');
            %验证条件2
            best_physical_std=physical_normal_H_measurement(best_physical_sensor_index,3);
            best_physical_H_decline=physical_burst_H_decline(best_physical_sensor_index,1+current_burst_pipe_index);
            H_decline_ratio=mut_virtual_H_decline/best_physical_H_decline;
            if H_decline_ratio<(mut_virtual_std/best_physical_std) && abs(H_decline_ratio-(mut_virtual_std/best_physical_std))>0.0000001 
                void=1;
            else
                if H_decline_ratio>(mut_sum_physical_std/best_physical_std) && abs(H_decline_ratio-(mut_sum_physical_std/best_physical_std))>0.0000001
                    void=1;
                end
            end
        end
    end
    if void==1
        virtual_matrix(i,:)=original_virtual;
    end
end
virtual_matrix_new = virtual_matrix;
end

%% 最差个体统计
function [virtual_matrix_new, fitness_new] = ReplaceWorse(virtual_matrix, current_virtual_best, fitness)
max_num = max(fitness);
min_num = min(fitness);
limit = (max_num-min_num)*0.20+min_num;
replace_corr = fitness<limit;
replace_num = sum(replace_corr);
virtual_matrix(replace_corr, :) = ones(replace_num, 1)*current_virtual_best(1:end-1);
fitness(replace_corr) = ones(replace_num, 1)*current_virtual_best(end);
virtual_matrix_new = virtual_matrix;
fitness_new = fitness;
end

%% 计算平均适应度
function fitness_ave = CalAveFitness(fitness)
[N ,~] = size(fitness);
fitness_ave = sum(fitness)/N;
end

