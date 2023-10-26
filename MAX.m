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
   L_pressure=data_extract(:,1:20);%左脚压力矩阵
   R_pressure=data_extract(:,21:end);%右脚压力矩阵
  L_pressure(54,14)=0;L_pressure(53,15)=0;L_pressure(52,16)=0;
  R_pressure(17,1)=0;R_pressure(52,5)=0;R_pressure(53,6)=0;R_pressure(54,7)=0;%去除异常点
   
   %% 接触面积
   points_L=nnz(L_pressure);%接触点数，计算矩阵中非零元素的数量
   L_points(i)=points_L;
   Area_L(i)=points_L*4.4*4.4;%一个点4.4mm*4.4mm

   points_R=nnz(R_pressure);%接触点数，计算矩阵中非零元素的数量
   R_points(i)=points_R;
   Area_R(i)=points_R*4.4*4.4;%一个点4.4mm*4.4mm
   %% 左脚分区,单位kPa
   L_T1=L_pressure(2:10,14:19);
   L_T2=L_pressure(2:10,1:13);
   L_M1=L_pressure(11:25,17:20);
   L_M2=L_pressure(11:25,13:16);
   L_M3=L_pressure(11:25,9:12);
   L_M4=L_pressure(11:25,5:8);
   L_M5=L_pressure(11:25,1:4);
   L_Midfoot_Lateral=L_pressure(26:38,1:9);
   L_Midfoot_Medial=L_pressure(26:38,10:19);
   L_Heel_Lateral=L_pressure(39:56,1:9);
   L_Heel_Medial=L_pressure(39:56,10:19);
   
   %% 右脚分区
   R_T1=R_pressure(2:10,1:6);
   R_T2=R_pressure(2:10,7:19);
   R_M1=L_pressure(11:25,1:4);
   R_M2=L_pressure(11:25,5:8);
   R_M3=L_pressure(11:25,9:12);
   R_M4=L_pressure(11:25,13:16);
   R_M5=L_pressure(11:25,17:20);
   R_Midfoot_Lateral=R_pressure(26:38,10:19);
   R_Midfoot_Medial=R_pressure(26:38,1:9);
   R_Heel_Lateral=R_pressure(39:56,10:19);
   R_Heel_Medial=R_pressure(39:56,1:9);
   
  %% 每个区的压强最大值
   L_T1_Max(i)=max(max(L_T1));
   L_T2_Max(i)=max(max(L_T2));
   L_M1_Max(i)=max(max(L_M1));
   L_M2_Max(i)=max(max(L_M2));
   L_M3_Max(i)=max(max(L_M3));
   L_M4_Max(i)=max(max(L_M4));
   L_M5_Max(i)=max(max(L_M5));
   L_Midfoot_Lateral_Max(i)=max(max(L_Midfoot_Lateral));
   L_Midfoot_Medial_Max(i)=max(max( L_Midfoot_Medial));
   L_Heel_Lateral_Max(i)=max(max( L_Heel_Lateral));
   L_Heel_Medial_Max(i)=max(max( L_Heel_Medial));
   L_Max(i)=max(max(L_pressure));
 
   R_T1_Max(i)=max(max(R_T1));
   R_T2_Max(i)=max(max(R_T2));
   R_M1_Max(i)=max(max(R_M1));
   R_M2_Max(i)=max(max(R_M2));
   R_M3_Max(i)=max(max(R_M3));
   R_M4_Max(i)=max(max(R_M4));
   R_M5_Max(i)=max(max(R_M5));
   R_Midfoot_Lateral_Max(i)=max(max(R_Midfoot_Lateral));
   R_Midfoot_Medial_Max(i)=max(max( R_Midfoot_Medial));
   R_Heel_Lateral_Max(i)=max(max( R_Heel_Lateral));
   R_Heel_Medial_Max(i)=max(max( R_Heel_Medial));
   R_Max(i)=max(max(R_pressure));
   
   %% 区域接触力,单位kg
Stress_L_T1(i)=sum(sum(L_T1.*19.36*0.0001));
Stress_L_T2(i)=sum(sum(L_T2.*19.36*0.0001));
Stress_L_M1(i)=sum(sum(L_M1.*19.36*0.001));
Stress_L_M2(i)=sum(sum(L_M2.*19.36*0.0001));
Stress_L_M3(i)=sum(sum(L_M3.*19.36*0.0001));
Stress_L_M4(i)=sum(sum(L_M4.*19.36*0.0001));
Stress_L_M5(i)=sum(sum(L_M5.*19.36*0.0001));
Stress_L_Midfoot_Lateral(i)=sum(sum(L_Midfoot_Lateral.*19.36*0.0001));
Stress_L_Midfoot_Medial(i)=sum(sum(L_Midfoot_Medial.*19.36*0.0001));
Stress_L_Heel_Lateral(i)=sum(sum(L_Heel_Lateral.*19.36*0.0001));
Stress_L_Heel_Medial(i)=sum(sum(L_Heel_Medial.*19.36*0.0001));
Stress_L(i)=sum(sum(L_pressure.*19.36*0.0001));

Stress_R_T1(i)=sum(sum(R_T1.*19.36*0.0001));
Stress_R_T2(i)=sum(sum(R_T2.*19.36*0.0001));
Stress_R_M1(i)=sum(sum(R_M1.*19.36*0.0001));
Stress_R_M2(i)=sum(sum(R_M2.*19.36*0.0001));
Stress_R_M3(i)=sum(sum(R_M3.*19.36*0.0001));
Stress_R_M4(i)=sum(sum(R_M4.*19.36*0.0001));
Stress_R_M5(i)=sum(sum(R_M5.*19.36*0.0001));
Stress_R_Midfoot_Lateral(i)=sum(sum(R_Midfoot_Lateral.*19.36*0.0001));
Stress_R_Midfoot_Medial(i)=sum(sum(R_Midfoot_Medial.*19.36*0.0001));
Stress_R_Heel_Lateral(i)=sum(sum(R_Heel_Lateral.*19.36*0.0001));
Stress_R_Heel_Medial(i)=sum(sum(R_Heel_Medial.*19.36*0.0001));
Stress_R(i)=sum(sum(R_pressure.*19.36*0.0001));

  
end
%% 整合
 Pressure_L_T1_Max=L_T1_Max(4:interval:end)';%每一帧的每个区域的压强的最大值（一个区域内一个点的压强）
 Pressure_L_T2_Max=L_T2_Max(4:interval:end)';
 Pressure_L_M1_Max=L_M1_Max(4:interval:end)';
 Pressure_L_M2_Max=L_M2_Max(4:interval:end)';
 Pressure_L_M3_Max=L_M3_Max(4:interval:end)';
 Pressure_L_M4_Max=L_M4_Max(4:interval:end)';
 Pressure_L_M5_Max=L_M5_Max(4:interval:end)';
 Pressure_L_Midfoot_Lateral_Max=L_Midfoot_Lateral_Max(4:interval:end)';
 Pressure_L_Midfoot_Medial_Max=L_Midfoot_Medial_Max(4:interval:end)';
 Pressure_L_Heel_Lateral_Max=L_Heel_Lateral_Max(4:interval:end)';
 Pressure_L_Heel_Medial_Max=L_Heel_Medial_Max(4:interval:end)';
 Pressure_L_Max=L_Max(4:interval:end)';
 
 Pressure_R_T1_Max=R_T1_Max(4:interval:end)';
 Pressure_R_T2_Max=R_T2_Max(4:interval:end)';
 Pressure_R_M1_Max=R_M1_Max(4:interval:end)';
 Pressure_R_M2_Max=R_M2_Max(4:interval:end)';
 Pressure_R_M3_Max=R_M3_Max(4:interval:end)';
 Pressure_R_M4_Max=R_M4_Max(4:interval:end)';
 Pressure_R_M5_Max=R_M5_Max(4:interval:end)';
 Pressure_R_Midfoot_Lateral_Max=R_Midfoot_Lateral_Max(4:interval:end)';
 Pressure_R_Midfoot_Medial_Max=R_Midfoot_Medial_Max(4:interval:end)';
 Pressure_R_Heel_Lateral_Max=R_Heel_Lateral_Max(4:interval:end)';
 Pressure_R_Heel_Medial_Max=R_Heel_Medial_Max(4:interval:end)';
 Pressure_R_Max=R_Max(4:interval:end)';
 
 
%% 接触力值
Stress_L_T1=Stress_L_T1(4:interval:end)';%每一帧的每个区域的接触力（一整个区域所受的总力）
Stress_L_T2=Stress_L_T2(4:interval:end)';
Stress_L_M1=Stress_L_M1(4:interval:end)';
Stress_L_M2=Stress_L_M2(4:interval:end)';
Stress_L_M3=Stress_L_M3(4:interval:end)';
Stress_L_M4=Stress_L_M4(4:interval:end)';
Stress_L_M5=Stress_L_M5(4:interval:end)';
Stress_L_Midfoot_Lateral=Stress_L_Midfoot_Lateral(4:interval:end)';
Stress_L_Midfoot_Medial=Stress_L_Midfoot_Medial(4:interval:end)';
Stress_L_Heel_Lateral=Stress_L_Heel_Lateral(4:interval:end)';
Stress_L_Heel_Medial=Stress_L_Heel_Medial(4:interval:end)';
Stress_L=Stress_L(4:interval:end)';

Stress_R_T1=Stress_R_T1(4:interval:end)';
Stress_R_T2=Stress_R_T2(4:interval:end)';
Stress_R_M1=Stress_R_M1(4:interval:end)';
Stress_R_M2=Stress_R_M2(4:interval:end)';
Stress_R_M3=Stress_R_M3(4:interval:end)';
Stress_R_M4=Stress_R_M4(4:interval:end)';
Stress_R_M5=Stress_R_M5(4:interval:end)';
Stress_R_Midfoot_Lateral=Stress_R_Midfoot_Lateral(4:interval:end)';
Stress_R_Midfoot_Medial=Stress_R_Midfoot_Medial(4:interval:end)';
Stress_R_Heel_Lateral=Stress_R_Heel_Lateral(4:interval:end)';
Stress_R_Heel_Medial=Stress_R_Heel_Medial(4:interval:end)';
Stress_R=Stress_R(4:interval:end)';

%% 接触面积,单位cm^2
   L_points=L_points(4:interval:end)';%左脚接触点数变化
   Area_L=Area_L(4:interval:end)'.*0.01;%左脚接触面积
   Area_R=Area_R(4:interval:end)'.*0.01;%右脚接触面积
%% 支撑相中最大值的判定
%通过时间判定哪些帧是支撑相，寻找最大值队列里面的相应帧，并对其进行对比分析。



%% 画图

% 
% plot(Pressure_L_T1_Max,'LineWidth', 1);
% title('压强值变化图','FontSize',16);
% xlabel('帧数','FontSize',16);
% ylabel('压强/kPa','FontSize',16);
% ax = gca;%获取当前坐标轴对象
% set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
% set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细
% 
% 
% plot(Area_R,'LineWidth', 1);
% title('接触面积变化','FontSize',16);
% xlabel('帧数','FontSize',16);
% ylabel('面积/cm^2','FontSize',16);
% ax = gca;%获取当前坐标轴对象
% set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
% set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细

plot(Stress_R);hold on
plot(Stress_R_Heel_Lateral);hold on
plot(Stress_R_Heel_Medial);hold on
plot(Stress_R_M1);hold on
plot(Stress_R_M2);hold on
plot(Stress_R_M3);hold on
plot(Stress_R_M4);hold on
plot(Stress_R_M5);hold on
plot(Stress_R_Midfoot_Lateral);hold on
plot(Stress_R_Midfoot_Medial);hold on
plot(Stress_R_T1);hold on
plot(Stress_R_T2);hold on
legend('总压力', '足后跟外侧','足后跟内侧','第一跖骨',...
    '第二跖骨','第三跖骨','第四跖骨','第五跖骨','足弓外侧',...
    '足弓内侧', '大拇指', '其余四指');
title('接触力变化','FontSize',16);
xlabel('帧数','FontSize',16);
ylabel('接触力/kg','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细