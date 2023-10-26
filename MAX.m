clc;
clear all;
%%
%读取
% filename=('E:\FTW\IMU_Space\IMU_Space\CeShi\FT1-test\lz2-1.csv');
filename=('E:\FTW\record\sry1011-70hz-21-1.csv');
[num,txt] = xlsread(filename);
 a=txt(length(txt)-2);
tt=strsplit(a{:,1},':');%在“：”处把时间字符串拆分
frame=tt(2);;%一共采集了多少帧数据；
frame=str2num(char(frame));
%%
%提取时间值
start_row_Time=10;%开始的行
interval=67;%间隔
num_data_points=frame;%要提取的数据点数量
selected_rows=start_row_Time+(0:num_data_points-1)*interval;%选择行的索引
selected_data_Time=txt(selected_rows,:);%提取数据
%%提取压力值
% n=56;%每次提取的行数
% start_row_Pressure=4;%开始的行
% end_row_Pressure=59;
% interval=67;%间隔
% num_data_points=frame;%要提取的数据点数量
% selected_rows=start_row_Pressure+(0:num_data_points-1)*interval;%选择行的索引
% selected_data_Pressure=num(selected_rows:selected_rows+n-1,:);%提取数据

%% 提取每一帧的压力数据
num_rows=size(num,1);
n=56;%每次提取的行数
for i=4:interval:num_rows-n+1
    %从第4行，第4+interval行开始提取数据
   data_extract=num(i:i+n-1,:);
   L_pressure=data_extract(:,1:19);%左脚压力矩阵
   R_pressure=data_extract(:,21:end);%右脚压力矩阵
   L_pressure(54,14)=0;L_pressure(53,15)=0;L_pressure(52,16)=0;%去除异常点
   
   %% 接触面积
   points_L=nnz(L_pressure);%接触点数，计算矩阵中非零元素的数量
   L_points(i)=points_L;
   Area_L(i)=points_L*4.4*4.4;%一个点4.4mm*4.4mm

   points_R=nnz(R_pressure);%接触点数，计算矩阵中非零元素的数量
   R_points(i)=points_R;
   Area_R(i)=points_R*4.4*4.4;%一个点4.4mm*4.4mm
   %% 左脚分区
   L_T1=L_pressure(2:10,14:19);
   L_T2=L_pressure(2:10,1:13);
   L_Forefoot_Lateral=L_pressure(11:25,1:9);
   L_Forefoot_Medial=L_pressure(11:25,10:19);
   L_Midfoot_Lateral=L_pressure(26:38,1:9);
   L_Midfoot_Medial=L_pressure(26:38,10:19);
   L_Heel_Lateral=L_pressure(39:56,1:9);
   L_Heel_Medial=L_pressure(39:56,10:19);
   
   %% 右脚分区
   R_T1=R_pressure(2:10,1:6);
   R_T2=R_pressure(2:10,7:19);
   R_Forefoot_Lateral=R_pressure(11:25,10:19);
   R_Forefoot_Medial=R_pressure(11:25,1:9);
   R_Midfoot_Lateral=R_pressure(26:38,10:19);
   R_Midfoot_Medial=R_pressure(26:38,1:9);
   R_Heel_Lateral=R_pressure(39:56,10:19);
   R_Heel_Medial=R_pressure(39:56,1:9);
   
  %% 每个区的压强最大值
   L_T1_Max1=max(max(L_T1));
   L_T2_Max1=max(max(L_T2));
   L_Forefoot_Lateral_Max1=max(max(L_Forefoot_Lateral));
   L_Forefoot_Medial_Max1=max(max( L_Forefoot_Medial));
   L_Midfoot_Lateral_Max1=max(max(L_Midfoot_Lateral));
   L_Midfoot_Medial_Max1=max(max( L_Midfoot_Medial));
   L_Heel_Lateral_Max1=max(max( L_Heel_Lateral));
   L_Heel_Medial_Max1=max(max( L_Heel_Medial));

 L_T1_Max(i)=L_T1_Max1;
 L_T2_Max(i)=L_T2_Max1;
 L_Forefoot_Lateral_Max(i)=L_Forefoot_Lateral_Max1;
 L_Forefoot_Medial_Max(i)=L_Forefoot_Medial_Max1;
 L_Midfoot_Lateral_Max(i)=L_Midfoot_Lateral_Max1;
 L_Midfoot_Medial_Max(i)=L_Midfoot_Medial_Max1;
 L_Heel_Lateral_Max(i)=L_Heel_Lateral_Max1;
 L_Heel_Medial_Max(i)=L_Heel_Medial_Max1;
L_Max(i)=max(max(L_pressure));
 
   R_T1_Max1(i)=max(max(R_T1));
   R_T2_Max1(i)=max(max(R_T2));
   R_Forefoot_Lateral_Max1(i)=max(max(R_Forefoot_Lateral));
   R_Forefoot_Medial_Max1(i)=max(max( R_Forefoot_Medial));
   R_Midfoot_Lateral_Max1(i)=max(max(R_Midfoot_Lateral));
   R_Midfoot_Medial_Max1(i)=max(max( R_Midfoot_Medial));
   R_Heel_Lateral_Max1(i)=max(max( R_Heel_Lateral));
   R_Heel_Medial_Max1(i)=max(max( R_Heel_Medial));
   R_Max(i)=max(max(R_pressure));
  
end
%% 整合
 Pressure_L_T1_Max=L_T1_Max(4:interval:end)';%每一帧的每个区域的压强的最大值
 Pressure_L_T2_Max=L_T2_Max(4:interval:end)';
 Pressure_L_Forefoot_Lateral_Max=L_Forefoot_Lateral_Max(4:interval:end)';
 Pressure_L_Forefoot_Medial_Max=L_Forefoot_Medial_Max(4:interval:end)';
 Pressure_L_Midfoot_Lateral_Max=L_Midfoot_Lateral_Max(4:interval:end)';
 Pressure_L_Midfoot_Medial_Max=L_Midfoot_Medial_Max(4:interval:end)';
 Pressure_L_Heel_Lateral_Max=L_Heel_Lateral_Max(4:interval:end)';
 Pressure_L_Heel_Medial_Max=L_Heel_Medial_Max(4:interval:end)';
 Pressure_L_Max=L_Max(4:interval:end)';
 
 Pressure_R_T1_Max=R_T1_Max1(4:interval:end)';
 Pressure_R_T2_Max=R_T2_Max1(4:interval:end)';
 Pressure_R_Forefoot_Lateral_Max=R_Forefoot_Lateral_Max1(4:interval:end)';
 Pressure_R_Forefoot_Medial_Max=R_Forefoot_Medial_Max1(4:interval:end)';
 Pressure_R_Midfoot_Lateral_Max=R_Midfoot_Lateral_Max1(4:interval:end)';
 Pressure_R_Midfoot_Medial_Max=R_Midfoot_Medial_Max1(4:interval:end)';
 Pressure_R_Heel_Lateral_Max=R_Heel_Lateral_Max1(4:interval:end)';
 Pressure_R_Heel_Medial_Max=R_Heel_Medial_Max1(4:interval:end)';
 Pressure_R_Max=R_Max(4:interval:end)';
%% 压力值
Stress_L_T1_Max=Pressure_L_T1_Max*4.4*4.4;
Stress_L_T2_Max=Pressure_L_T2_Max*4.4*4.4;
Stress_L_Forefoot_Lateral_Max=Pressure_L_Forefoot_Lateral_Max*4.4*4.4;
Stress_L_Forefoot_Medial_Max=Pressure_L_Forefoot_Medial_Max*4.4*4.4;
Stress_L_Midfoot_Lateral_Max=Pressure_L_Midfoot_Lateral_Max*4.4*4.4;
Stress_L_Midfoot_Medial_Max=Pressure_L_Midfoot_Medial_Max*4.4*4.4;
Stress_L_Heel_Lateral_Max=Pressure_L_Heel_Lateral_Max*4.4*4.4;
Stress_L_Heel_Medial_Max=Pressure_L_Heel_Medial_Max*4.4*4.4;
Stress_L_Max=Pressure_L_Max*4.4*4.4;

Stress_R_T1_Max=Pressure_R_T1_Max*4.4*4.4;
Stress_R_T2_Max=Pressure_R_T2_Max*4.4*4.4;
Stress_R_Forefoot_Lateral_Max=Pressure_R_Forefoot_Lateral_Max*4.4*4.4;
Stress_R_Forefoot_Medial_Max=Pressure_R_Forefoot_Medial_Max*4.4*4.4;
Stress_R_Midfoot_Lateral_Max=Pressure_R_Midfoot_Lateral_Max*4.4*4.4;
Stress_R_Midfoot_Medial_Max=Pressure_R_Midfoot_Medial_Max*4.4*4.4;
Stress_R_Heel_Lateral_Max=Pressure_R_Heel_Lateral_Max*4.4*4.4;
Stress_R_Heel_Medial_Max=Pressure_R_Heel_Medial_Max*4.4*4.4;
Stress_R_Max=Pressure_R_Max*4.4*4.4;

   L_points=L_points(4:interval:end)';%左脚接触点数变化
   Area_L=Area_L(4:interval:end)';%左脚接触面积
   Area_R=Area_R(4:interval:end)';%右脚接触面积
%% 支撑相中最大值的判定
%通过时间判定哪些帧是支撑相，寻找最大值队列里面的相应帧，并对其进行对比分析。



%% 画图


plot(Pressure_L_T1_Max,'LineWidth', 1);
title('压强值变化图','FontSize',16);
xlabel('帧数','FontSize',16);
ylabel('压强/kPa','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细


plot(Area_R,'LineWidth', 1);
title('接触面积变化','FontSize',16);
xlabel('帧数','FontSize',16);
ylabel('面积/mm^2','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细