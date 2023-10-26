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
   L_pressure=data_extract(:,1:20);%���ѹ������
   R_pressure=data_extract(:,21:end);%�ҽ�ѹ������
  L_pressure(54,14)=0;L_pressure(53,15)=0;L_pressure(52,16)=0;
  R_pressure(17,1)=0;R_pressure(52,5)=0;R_pressure(53,6)=0;R_pressure(54,7)=0;%ȥ���쳣��
   
   %% �Ӵ����
   points_L=nnz(L_pressure);%�Ӵ���������������з���Ԫ�ص�����
   L_points(i)=points_L;
   Area_L(i)=points_L*4.4*4.4;%һ����4.4mm*4.4mm

   points_R=nnz(R_pressure);%�Ӵ���������������з���Ԫ�ص�����
   R_points(i)=points_R;
   Area_R(i)=points_R*4.4*4.4;%һ����4.4mm*4.4mm
   %% ��ŷ���,��λkPa
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
   
   %% �ҽŷ���
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
   
  %% ÿ������ѹǿ���ֵ
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
   
   %% ����Ӵ���,��λkg
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
%% ����
 Pressure_L_T1_Max=L_T1_Max(4:interval:end)';%ÿһ֡��ÿ�������ѹǿ�����ֵ��һ��������һ�����ѹǿ��
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
 
 
%% �Ӵ���ֵ
Stress_L_T1=Stress_L_T1(4:interval:end)';%ÿһ֡��ÿ������ĽӴ�����һ�����������ܵ�������
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

%% �Ӵ����,��λcm^2
   L_points=L_points(4:interval:end)';%��ŽӴ������仯
   Area_L=Area_L(4:interval:end)'.*0.01;%��ŽӴ����
   Area_R=Area_R(4:interval:end)'.*0.01;%�ҽŽӴ����
%% ֧���������ֵ���ж�
%ͨ��ʱ���ж���Щ֡��֧���࣬Ѱ�����ֵ�����������Ӧ֡����������жԱȷ�����



%% ��ͼ

% 
% plot(Pressure_L_T1_Max,'LineWidth', 1);
% title('ѹǿֵ�仯ͼ','FontSize',16);
% xlabel('֡��','FontSize',16);
% ylabel('ѹǿ/kPa','FontSize',16);
% ax = gca;%��ȡ��ǰ���������
% set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
% set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ
% 
% 
% plot(Area_R,'LineWidth', 1);
% title('�Ӵ�����仯','FontSize',16);
% xlabel('֡��','FontSize',16);
% ylabel('���/cm^2','FontSize',16);
% ax = gca;%��ȡ��ǰ���������
% set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
% set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ

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
legend('��ѹ��', '�������','�����ڲ�','��һ�Ź�',...
    '�ڶ��Ź�','�����Ź�','�����Ź�','�����Ź�','�㹭���',...
    '�㹭�ڲ�', '��Ĵָ', '������ָ');
title('�Ӵ����仯','FontSize',16);
xlabel('֡��','FontSize',16);
ylabel('�Ӵ���/kg','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ