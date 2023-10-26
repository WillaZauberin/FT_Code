clc;
clear all;
%%
%��ȡ
% filename=('E:\FTW\IMU_Space\IMU_Space\CeShi\FT1-test\lz2-1.csv');
filename=('E:\FTW\record\sry1011-70hz-21-1.csv');
[num,txt] = xlsread(filename);
 a=txt(length(txt)-2);
tt=strsplit(a{:,1},':');%�ڡ���������ʱ���ַ������
frame=tt(2);;%һ���ɼ��˶���֡���ݣ�
frame=str2num(char(frame));
%%
%��ȡʱ��ֵ
start_row_Time=10;%��ʼ����
interval=67;%���
num_data_points=frame;%Ҫ��ȡ�����ݵ�����
selected_rows=start_row_Time+(0:num_data_points-1)*interval;%ѡ���е�����
selected_data_Time=txt(selected_rows,:);%��ȡ����
%%��ȡѹ��ֵ
% n=56;%ÿ����ȡ������
% start_row_Pressure=4;%��ʼ����
% end_row_Pressure=59;
% interval=67;%���
% num_data_points=frame;%Ҫ��ȡ�����ݵ�����
% selected_rows=start_row_Pressure+(0:num_data_points-1)*interval;%ѡ���е�����
% selected_data_Pressure=num(selected_rows:selected_rows+n-1,:);%��ȡ����

%% ��ȡÿһ֡��ѹ������
num_rows=size(num,1);
n=56;%ÿ����ȡ������
for i=4:interval:num_rows-n+1
    %�ӵ�4�У���4+interval�п�ʼ��ȡ����
   data_extract=num(i:i+n-1,:);
   L_pressure=data_extract(:,1:19);%���ѹ������
   R_pressure=data_extract(:,21:end);%�ҽ�ѹ������
   L_pressure(54,14)=0;L_pressure(53,15)=0;L_pressure(52,16)=0;%ȥ���쳣��
   
   %% �Ӵ����
   points_L=nnz(L_pressure);%�Ӵ���������������з���Ԫ�ص�����
   L_points(i)=points_L;
   Area_L(i)=points_L*4.4*4.4;%һ����4.4mm*4.4mm

   points_R=nnz(R_pressure);%�Ӵ���������������з���Ԫ�ص�����
   R_points(i)=points_R;
   Area_R(i)=points_R*4.4*4.4;%һ����4.4mm*4.4mm
   %% ��ŷ���
   L_T1=L_pressure(2:10,14:19);
   L_T2=L_pressure(2:10,1:13);
   L_Forefoot_Lateral=L_pressure(11:25,1:9);
   L_Forefoot_Medial=L_pressure(11:25,10:19);
   L_Midfoot_Lateral=L_pressure(26:38,1:9);
   L_Midfoot_Medial=L_pressure(26:38,10:19);
   L_Heel_Lateral=L_pressure(39:56,1:9);
   L_Heel_Medial=L_pressure(39:56,10:19);
   
   %% �ҽŷ���
   R_T1=R_pressure(2:10,1:6);
   R_T2=R_pressure(2:10,7:19);
   R_Forefoot_Lateral=R_pressure(11:25,10:19);
   R_Forefoot_Medial=R_pressure(11:25,1:9);
   R_Midfoot_Lateral=R_pressure(26:38,10:19);
   R_Midfoot_Medial=R_pressure(26:38,1:9);
   R_Heel_Lateral=R_pressure(39:56,10:19);
   R_Heel_Medial=R_pressure(39:56,1:9);
   
  %% ÿ������ѹǿ���ֵ
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
%% ����
 Pressure_L_T1_Max=L_T1_Max(4:interval:end)';%ÿһ֡��ÿ�������ѹǿ�����ֵ
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
%% ѹ��ֵ
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

   L_points=L_points(4:interval:end)';%��ŽӴ������仯
   Area_L=Area_L(4:interval:end)';%��ŽӴ����
   Area_R=Area_R(4:interval:end)';%�ҽŽӴ����
%% ֧���������ֵ���ж�
%ͨ��ʱ���ж���Щ֡��֧���࣬Ѱ�����ֵ�����������Ӧ֡����������жԱȷ�����



%% ��ͼ


plot(Pressure_L_T1_Max,'LineWidth', 1);
title('ѹǿֵ�仯ͼ','FontSize',16);
xlabel('֡��','FontSize',16);
ylabel('ѹǿ/kPa','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ


plot(Area_R,'LineWidth', 1);
title('�Ӵ�����仯','FontSize',16);
xlabel('֡��','FontSize',16);
ylabel('���/mm^2','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ