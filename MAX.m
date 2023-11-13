clc;
clear all;
%%
%读取
% filename=('E:\FTW\IMU_Space\IMU_Space\CeShi\FT1-test\lz2-1.csv');
filename=('E:\FTW\record\DOE_pressure\fama\3-6-1.csv');
% filename=('E:\FTW\record\Area\dct_R-4-1.csv');
% filename=('E:\FTW\record\DOE_pressure\11qu\zhf-2-1.csv');
% filename=('E:\FTW\record\DOE_pressure\fenqu\sry0926-70hz-12-1.csv');
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
%% 提取每一帧的压力数据
num_rows=size(num,1);
n=56;%每次提取的行数
%% 提取时间
start_row_Time=10;%开始的行
interval=67;%间隔
num_data_points=frame;%要提取的数据点数量
selected_rows=start_row_Time+(0:num_data_points-1)*interval;%选择行的索引
selected_data_Time=txt(selected_rows,:);%提取数据
t=selected_data_Time(:,1);%WIFI版的IMU时间数据
numt=Time(t,frame);%单位s
numt=numt*1000;%单位ms
%% 提取有数据的行数和列数（开始的行列、结束的行列），以已有的数据进行比例划分区域
for i=4:interval:num_rows-n+1
   %从第4行，第4+interval行开始提取数据
   data_extract=num(i:i+n-1,:);
   L_pressure=data_extract(:,1:20);%左脚压力矩阵
   R_pressure=data_extract(:,21:end);%右脚压力矩阵
%% 行列
%    [row_L,col_L]=find(L_pressure);
%    [row_R,col_R]=find(R_pressure);
% 
%    if isempty(row_L)
%     % 如果find没有找到满足条件的元素，给数组一个特定的值
%     row_L = [row_L,5]; % 为数组A添加一个5
%    end
%    if isempty(col_L)
%     col_L = [col_L,5]; 
%    end
%    if isempty(row_R)
%     row_R = [row_R,5]; 
%    end
%       if isempty(col_R)
%     col_R = [col_R,5]; 
%      end
%    
%    Start_Row_L=min(row_L); %开始的行列是有数据的行列数的最小值
%    Start_Col_L=min(col_L);
%    Start_Row_R=min(row_R);
%    Start_Col_R=min(col_R);
%    
%    Start_Row_L1(i)=Start_Row_L; %开始的行列是有数据的行列数的最小值
%    Start_Col_L1(i)=Start_Col_L; 
%    Start_Row_R1(i)=Start_Row_R; %开始的行列是有数据的行列数的最小值
%    Start_Col_R1(i)=Start_Col_R; 
%   
%    End_Row_L=max(row_L); %结束的行列是有数据的行列数的最大值
%    End_Col_L=max(col_L);
%    End_Row_R=max(row_R);
%    End_Col_R=max(col_R);
%    
%    End_Row_L1(i)=End_Row_L; %开始的行列是有数据的行列数的最小值
%    End_Col_L1(i)=End_Col_L; 
%    End_Row_R1(i)=End_Row_R; %开始的行列是有数据的行列数的最小值
%    End_Col_R1(i)=End_Col_R; 
% fama
  L_pressure=data_extract(:,1:20);%左脚压力矩阵
 Stress_L(i)=sum(sum(L_pressure.*19.36*0.0001*10));
   %% 接触面积
   points_L=nnz(L_pressure);%接触点数，计算矩阵中非零元素的数量
   L_points(i)=points_L;
   Area_L(i)=points_L*4.4*4.4;%一个点4.4mm*4.4mm

   points_R=nnz(R_pressure);%接触点数，计算矩阵中非零元素的数量
   R_points(i)=points_R;
   Area_R(i)=points_R*4.4*4.4;%一个点4.4mm*4.4mm
end
% fama
Stress_L=Stress_L(4:interval:end)';

%% 接触面积,单位cm^2
   L_points=L_points(4:interval:end)';%左脚接触点数变化
   Area_L=Area_L(4:interval:end)'.*0.01;%左脚接触面积
   Area_R=Area_R(4:interval:end)'.*0.01;%右脚接触面积
 
   
   %% 整合行列
   Start_Row_L1=Start_Row_L1(4:interval:end)';
   Start_Col_L1=Start_Col_L1(4:interval:end)';
   Start_Row_R1=Start_Row_R1(4:interval:end)'; %开始的行列是有数据的行列数的最小值
   Start_Col_R1=Start_Col_R1(4:interval:end)';
   End_Row_L1=End_Row_L1(4:interval:end)';
   End_Col_L1=End_Col_L1(4:interval:end)';
   End_Row_R1=End_Row_R1(4:interval:end)'; %开始的行列是有数据的行列数的最小值
   End_Col_R1=End_Col_R1(4:interval:end)';   
   %% 获取足部接触面积开始的行列
   startrow_L=min(Start_Row_L1);
   startcol_L=min(Start_Col_L1);
   startrow_R=min(Start_Row_R1);
   startcol_R=min(Start_Col_R1);
   
   endrow_L=max(End_Row_L1);
   endcol_L=max(End_Col_L1);
   endrow_R=max(End_Row_R1);
   endcol_R=max(End_Col_R1);   
   
   row_L=endrow_L-startrow_L+1;%统计有效足压矩阵一共多少行，多少列
   row_R=endrow_R-startrow_R+1;
   col_L=endcol_L-startcol_L+1;
   col_R=endcol_R-startcol_R+1;
   
  %% 得到接触的面积行列数之后进行分区和接触力的计算 
   for i=4:interval:num_rows-n+1
    %从第4行，第4+interval行开始提取数据
   data_extract=num(i:i+n-1,:);
   L_pressure=data_extract(:,1:20);%左脚压力矩阵
   R_pressure=data_extract(:,21:end);%右脚压力矩阵 
   L_pressure(54,14)=0;L_pressure(53,15)=0;L_pressure(52,16)=0;
   R_pressure(17,1)=0;R_pressure(52,5)=0;R_pressure(53,6)=0;R_pressure(54,7)=0;%去除异常点
   R_pressure(20,2)=0;   R_pressure(27,18)=0;
   %% 左脚分区,单位（kPa）
   L_T1=L_pressure(startrow_L:round(1/5* row_L),round(2/3*col_L)+1:endcol_L);
   L_T2=L_pressure(startrow_L:round(1/5* row_L),startcol_L:round(2/3*col_L));
   L_M5=L_pressure(round(1/5*row_L)+1:round(2/5*row_L),startcol_L:round(12/50*col_L));
   L_M4=L_pressure(round(1/5*row_L)+1:round(2/5*row_L),round(12/50*col_L)+1:round(19/50*col_L));
   L_M3=L_pressure(round(1/5*row_L)+1:round(2/5*row_L),round(19/50*col_L)+1:round(26/50*col_L));
   L_M2=L_pressure(round(1/5*row_L)+1:round(2/5*row_L),round(26/50*col_L)+1:round(37/50*col_L));
   L_M1=L_pressure(round(1/5*row_L)+1:round(2/5*row_L),round(37/50*col_L)+1:endcol_L);
   L_Midfoot_Lateral=L_pressure(round(2/5*row_L)+1:round(7/10*row_L),startcol_L:round(1/2*col_L));
   L_Midfoot_Medial=L_pressure(round(2/5*row_L)+1:round(7/10*row_L),round(1/2*col_L)+1:endcol_L);
   L_Heel_Lateral=L_pressure(round(7/10*row_L)+1:endrow_L,startcol_L:round(1/2*col_L));
   L_Heel_Medial=L_pressure(round(7/10*row_L)+1:endrow_L,round(1/2*col_L)+1:endcol_L);
   
   %% 右脚分区
   R_T1=R_pressure(startrow_R:round(1/5* row_R),startcol_R:round(2/3*col_R));
   R_T2=R_pressure(startrow_R:round(1/5* row_R),round(2/3*col_R)+1:endcol_R);
   R_M1=R_pressure(round(1/5*row_R)+1:round(2/5*row_R),startcol_R:round(17/50*col_R));
   R_M2=R_pressure(round(1/5*row_R)+1:round(2/5*row_R),round(17/50*col_R)+1:round(24/50*col_R));
   R_M3=R_pressure(round(1/5*row_R)+1:round(2/5*row_R),round(24/50*col_R)+1:round(31/50*col_R));
   R_M4=R_pressure(round(1/5*row_R)+1:round(2/5*row_R),round(31/50*col_R)+1:round(38/50*col_R));
   R_M5=R_pressure(round(1/5*row_R)+1:round(2/5*row_R),round(38/50*col_R)+1:endcol_R);
   R_Midfoot_Medial=R_pressure(round(2/5*row_R)+1:round(7/10*row_R),startcol_R:round(1/2*col_R));
   R_Midfoot_Lateral=R_pressure(round(2/5*row_R)+1:round(7/10*row_R),round(1/2*col_R)+1:endcol_R);
   R_Heel_Medial=R_pressure(round(7/10*row_R)+1:endrow_R,startcol_R:round(1/2*col_R));   
   R_Heel_Lateral=R_pressure(round(7/10*row_R)+1:endrow_R,round(1/2*col_R)+1:endcol_R);
 
   
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
Stress_L_M1(i)=sum(sum(L_M1.*19.36*0.0001));
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
 
 
%% 接触力值（kg）
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


%% 支撑相中最大值的判定
%通过时间判定哪些帧是支撑相，寻找最大值队列里面的相应帧，并对其进行对比分析。



%% 画图

% 
% plot(Stress_L*10,'LineWidth', 1);
% title('接触力','FontSize',16);
% xlabel('帧数','FontSize',16);
% ylabel('接触力/N','FontSize',16);
% ax = gca;%获取当前坐标轴对象
% set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
% set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细
% ylim([0, 50]);
% 
% plot(Area_R,'LineWidth', 1);
% title('接触面积变化','FontSize',16);
% xlabel('帧数','FontSize',16);
% ylabel('面积/cm^2','FontSize',16);
% ax = gca;%获取当前坐标轴对象
% set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
% set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细

%% 右
plot(numt,Stress_R_Heel_Lateral,'LineWidth', 1);hold on
plot(numt,Stress_R_Heel_Medial,'LineWidth', 1);hold on
plot(numt,Stress_R_M1,'LineWidth', 1);hold on
plot(numt,Stress_R_M2,'LineWidth', 1);hold on
plot(numt,Stress_R_M3,'LineWidth', 1);hold on
plot(numt,Stress_R_M4,'LineWidth', 1);hold on
plot(numt,Stress_R_M5,'LineWidth', 1);hold on
plot(numt,Stress_R_Midfoot_Lateral,'LineWidth', 1);hold on
plot(numt,Stress_R_Midfoot_Medial,'LineWidth', 1);hold on
plot(numt,Stress_R_T1,'LineWidth', 1);hold on
plot(numt,Stress_R_T2,'LineWidth', 1);hold on
legend( '足后跟外侧','足后跟内侧','第一跖骨',...
    '第二跖骨','第三跖骨','第四跖骨','第五跖骨','足弓外侧',...
    '足弓内侧', '大拇指', '其余四指');
title('右脚区域接触力变化','FontSize',16);
xlabel('毫秒（ms）','FontSize',16);
ylabel('接触力/kg','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细
ylim([0, 35]);

%% 左
figure;%再建立一个窗口

plot(numt,Stress_L_Heel_Lateral,'LineWidth', 1);hold on
plot(numt,Stress_L_Heel_Medial,'LineWidth', 1);hold on
plot(numt,Stress_L_M1,'LineWidth', 1);hold on
plot(numt,Stress_L_M2,'LineWidth', 1);hold on
plot(numt,Stress_L_M3,'LineWidth', 1);hold on
plot(numt,Stress_L_M4,'LineWidth', 1);hold on
plot(numt,Stress_L_M5,'LineWidth', 1);hold on
plot(numt,Stress_L_Midfoot_Lateral,'LineWidth', 1);hold on
plot(numt,Stress_L_Midfoot_Medial,'LineWidth', 1);hold on
plot(numt,Stress_L_T1,'LineWidth', 1);hold on
plot(numt,Stress_L_T2,'LineWidth', 1);hold on
legend( '足后跟外侧','足后跟内侧','第一跖骨',...
    '第二跖骨','第三跖骨','第四跖骨','第五跖骨','足弓外侧',...
    '足弓内侧', '大拇指', '其余四指');
title('左脚区域接触力变化','FontSize',16);
xlabel('毫秒（ms）','FontSize',16);
ylabel('接触力/kg','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细
ylim([0, 35]);

%% 截取曲线的一部分（右R）
figure
x_min =1520;
x_max =2520;
% 找出要截取的数据点索引
indices_to_keep = (numt>= x_min) & (numt<= x_max);
color_T1=[162/255,32/255,236/255];
color_T2=[255/255,0/255,255/255];
color_M1=[254/255,190/255,251/255];
color_M2=[254/255,164/255,0/255];
color_M3=[0/255,253/255,2/255];
color_M4=[2/255,100/255,1/255];
color_M5=[0/255,2/255,254/255];
color_CM=[135/255,136/255,138/255];
color_CL=[135/255,206/255,251/255];
color_HM=[3/255,138/255,142/255];
color_HL=[136/255,71/255,17/255];
% 仅绘制被截取的数据点
plot(numt(indices_to_keep), Stress_R_Heel_Lateral(indices_to_keep),'LineWidth', 1,'Color',color_HL);hold on
plot(numt(indices_to_keep), Stress_R_Heel_Medial(indices_to_keep),'LineWidth', 1,'Color',color_HM);hold on
plot(numt(indices_to_keep), Stress_R_M1(indices_to_keep),'LineWidth', 1,'Color',color_M1);hold on
plot(numt(indices_to_keep), Stress_R_M2(indices_to_keep),'LineWidth', 1,'Color',color_M2);hold on
plot(numt(indices_to_keep), Stress_R_M3(indices_to_keep),'LineWidth', 1,'Color',color_M3);hold on
plot(numt(indices_to_keep), Stress_R_M4(indices_to_keep),'LineWidth', 1,'Color',color_M4);hold on
plot(numt(indices_to_keep), Stress_R_M5(indices_to_keep),'LineWidth', 1,'Color',color_M5);hold on
plot(numt(indices_to_keep), Stress_R_Midfoot_Lateral(indices_to_keep),'LineWidth', 1,'Color',color_CL);hold on
plot(numt(indices_to_keep), Stress_R_Midfoot_Medial(indices_to_keep),'LineWidth', 1,'Color',color_CM);hold on
plot(numt(indices_to_keep), Stress_R_T1(indices_to_keep),'LineWidth', 1,'Color',color_T2);hold on
plot(numt(indices_to_keep), Stress_R_T2(indices_to_keep),'LineWidth', 1,'Color',color_T1);hold on

% 自定义 x 轴刻度值和标签
custom_x_values = [x_min:250: x_max];
custom_x_labels = {'0', '250', '500','750','1000'};
% 设置 x 轴刻度值和标签
xticks(custom_x_values);
xticklabels(custom_x_labels);
legend( '足后跟外侧','足后跟内侧','第一跖骨',...
    '第二跖骨','第三跖骨','第四跖骨','第五跖骨','足弓外侧',...
    '足弓内侧', '大拇指', '其余四指');
title('右脚区域接触力变化','FontSize',16);
xlabel('毫秒（ms）','FontSize',16);
ylabel('接触力/kg','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细
ylim([0, 30]);

%% 截取曲线的一部分（左L）
figure
x_min=3880;
x_max=4670;
%找出要截取的数据点索引
indices_to_keep = (numt>= x_min) & (numt<= x_max);
color_T1=[162/255,32/255,236/255];
color_T2=[255/255,0/255,255/255];
color_M1=[254/255,190/255,251/255];
color_M2=[254/255,164/255,0/255];
color_M3=[0/255,253/255,2/255];
color_M4=[2/255,100/255,1/255];
color_M5=[0/255,2/255,254/255];
color_CM=[135/255,136/255,138/255];
color_CL=[135/255,206/255,251/255];
color_HM=[3/255,138/255,142/255];
color_HL=[136/255,71/255,17/255];
% 仅绘制被截取的数据点
plot(numt(indices_to_keep), Stress_L_Heel_Lateral(indices_to_keep),'LineWidth', 1,'Color',color_HL);hold on
plot(numt(indices_to_keep), Stress_L_Heel_Medial(indices_to_keep),'LineWidth', 1,'Color',color_HM);hold on
plot(numt(indices_to_keep), Stress_L_M1(indices_to_keep),'LineWidth', 1,'Color',color_M1);hold on
plot(numt(indices_to_keep), Stress_L_M2(indices_to_keep),'LineWidth', 1,'Color',color_M2);hold on
plot(numt(indices_to_keep), Stress_L_M3(indices_to_keep),'LineWidth', 1,'Color',color_M3);hold on
plot(numt(indices_to_keep), Stress_L_M4(indices_to_keep),'LineWidth', 1,'Color',color_M4);hold on
plot(numt(indices_to_keep), Stress_L_M5(indices_to_keep),'LineWidth', 1,'Color',color_M5);hold on
plot(numt(indices_to_keep), Stress_L_Midfoot_Lateral(indices_to_keep),'LineWidth', 1,'Color',color_CL);hold on
plot(numt(indices_to_keep), Stress_L_Midfoot_Medial(indices_to_keep),'LineWidth', 1,'Color',color_CM);hold on
plot(numt(indices_to_keep), Stress_L_T1(indices_to_keep),'LineWidth', 1,'Color',color_T2);hold on
plot(numt(indices_to_keep), Stress_L_T2(indices_to_keep),'LineWidth', 1,'Color',color_T1);hold on

% 自定义 x 轴刻度值和标签
custom_x_values = [x_min:197.5: x_max];
custom_x_labels = {'0', '197.5', '395','592.5','790'};
% 设置 x 轴刻度值和标签
xticks(custom_x_values);
xticklabels(custom_x_labels);
legend( '足后跟外侧','足后跟内侧','第一跖骨',...
    '第二跖骨','第三跖骨','第四跖骨','第五跖骨','足弓外侧',...
    '足弓内侧', '大拇指', '其余四指');
title('左脚区域接触力变化','FontSize',16);
xlabel('毫秒（ms）','FontSize',16);
ylabel('接触力/kg','FontSize',16);
ax = gca;%获取当前坐标轴对象
set(ax.XAxis, 'LineWidth', 1);  % 设置X轴的粗细
set(ax.YAxis, 'LineWidth', 1);  % 设置Y轴的粗细
ylim([0,22]);