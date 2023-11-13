clc;
clear all;
%%
%��ȡ
% filename=('E:\FTW\IMU_Space\IMU_Space\CeShi\FT1-test\lz2-1.csv');
filename=('E:\FTW\record\DOE_pressure\fama\3-6-1.csv');
% filename=('E:\FTW\record\Area\dct_R-4-1.csv');
% filename=('E:\FTW\record\DOE_pressure\11qu\zhf-2-1.csv');
% filename=('E:\FTW\record\DOE_pressure\fenqu\sry0926-70hz-12-1.csv');
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
%% ��ȡÿһ֡��ѹ������
num_rows=size(num,1);
n=56;%ÿ����ȡ������
%% ��ȡʱ��
start_row_Time=10;%��ʼ����
interval=67;%���
num_data_points=frame;%Ҫ��ȡ�����ݵ�����
selected_rows=start_row_Time+(0:num_data_points-1)*interval;%ѡ���е�����
selected_data_Time=txt(selected_rows,:);%��ȡ����
t=selected_data_Time(:,1);%WIFI���IMUʱ������
numt=Time(t,frame);%��λs
numt=numt*1000;%��λms
%% ��ȡ�����ݵ���������������ʼ�����С����������У��������е����ݽ��б�����������
for i=4:interval:num_rows-n+1
   %�ӵ�4�У���4+interval�п�ʼ��ȡ����
   data_extract=num(i:i+n-1,:);
   L_pressure=data_extract(:,1:20);%���ѹ������
   R_pressure=data_extract(:,21:end);%�ҽ�ѹ������
%% ����
%    [row_L,col_L]=find(L_pressure);
%    [row_R,col_R]=find(R_pressure);
% 
%    if isempty(row_L)
%     % ���findû���ҵ�����������Ԫ�أ�������һ���ض���ֵ
%     row_L = [row_L,5]; % Ϊ����A���һ��5
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
%    Start_Row_L=min(row_L); %��ʼ�������������ݵ�����������Сֵ
%    Start_Col_L=min(col_L);
%    Start_Row_R=min(row_R);
%    Start_Col_R=min(col_R);
%    
%    Start_Row_L1(i)=Start_Row_L; %��ʼ�������������ݵ�����������Сֵ
%    Start_Col_L1(i)=Start_Col_L; 
%    Start_Row_R1(i)=Start_Row_R; %��ʼ�������������ݵ�����������Сֵ
%    Start_Col_R1(i)=Start_Col_R; 
%   
%    End_Row_L=max(row_L); %�����������������ݵ������������ֵ
%    End_Col_L=max(col_L);
%    End_Row_R=max(row_R);
%    End_Col_R=max(col_R);
%    
%    End_Row_L1(i)=End_Row_L; %��ʼ�������������ݵ�����������Сֵ
%    End_Col_L1(i)=End_Col_L; 
%    End_Row_R1(i)=End_Row_R; %��ʼ�������������ݵ�����������Сֵ
%    End_Col_R1(i)=End_Col_R; 
% fama
  L_pressure=data_extract(:,1:20);%���ѹ������
 Stress_L(i)=sum(sum(L_pressure.*19.36*0.0001*10));
   %% �Ӵ����
   points_L=nnz(L_pressure);%�Ӵ���������������з���Ԫ�ص�����
   L_points(i)=points_L;
   Area_L(i)=points_L*4.4*4.4;%һ����4.4mm*4.4mm

   points_R=nnz(R_pressure);%�Ӵ���������������з���Ԫ�ص�����
   R_points(i)=points_R;
   Area_R(i)=points_R*4.4*4.4;%һ����4.4mm*4.4mm
end
% fama
Stress_L=Stress_L(4:interval:end)';

%% �Ӵ����,��λcm^2
   L_points=L_points(4:interval:end)';%��ŽӴ������仯
   Area_L=Area_L(4:interval:end)'.*0.01;%��ŽӴ����
   Area_R=Area_R(4:interval:end)'.*0.01;%�ҽŽӴ����
 
   
   %% ��������
   Start_Row_L1=Start_Row_L1(4:interval:end)';
   Start_Col_L1=Start_Col_L1(4:interval:end)';
   Start_Row_R1=Start_Row_R1(4:interval:end)'; %��ʼ�������������ݵ�����������Сֵ
   Start_Col_R1=Start_Col_R1(4:interval:end)';
   End_Row_L1=End_Row_L1(4:interval:end)';
   End_Col_L1=End_Col_L1(4:interval:end)';
   End_Row_R1=End_Row_R1(4:interval:end)'; %��ʼ�������������ݵ�����������Сֵ
   End_Col_R1=End_Col_R1(4:interval:end)';   
   %% ��ȡ�㲿�Ӵ������ʼ������
   startrow_L=min(Start_Row_L1);
   startcol_L=min(Start_Col_L1);
   startrow_R=min(Start_Row_R1);
   startcol_R=min(Start_Col_R1);
   
   endrow_L=max(End_Row_L1);
   endcol_L=max(End_Col_L1);
   endrow_R=max(End_Row_R1);
   endcol_R=max(End_Col_R1);   
   
   row_L=endrow_L-startrow_L+1;%ͳ����Ч��ѹ����һ�������У�������
   row_R=endrow_R-startrow_R+1;
   col_L=endcol_L-startcol_L+1;
   col_R=endcol_R-startcol_R+1;
   
  %% �õ��Ӵ������������֮����з����ͽӴ����ļ��� 
   for i=4:interval:num_rows-n+1
    %�ӵ�4�У���4+interval�п�ʼ��ȡ����
   data_extract=num(i:i+n-1,:);
   L_pressure=data_extract(:,1:20);%���ѹ������
   R_pressure=data_extract(:,21:end);%�ҽ�ѹ������ 
   L_pressure(54,14)=0;L_pressure(53,15)=0;L_pressure(52,16)=0;
   R_pressure(17,1)=0;R_pressure(52,5)=0;R_pressure(53,6)=0;R_pressure(54,7)=0;%ȥ���쳣��
   R_pressure(20,2)=0;   R_pressure(27,18)=0;
   %% ��ŷ���,��λ��kPa��
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
   
   %% �ҽŷ���
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
 
 
%% �Ӵ���ֵ��kg��
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


%% ֧���������ֵ���ж�
%ͨ��ʱ���ж���Щ֡��֧���࣬Ѱ�����ֵ�����������Ӧ֡����������жԱȷ�����



%% ��ͼ

% 
% plot(Stress_L*10,'LineWidth', 1);
% title('�Ӵ���','FontSize',16);
% xlabel('֡��','FontSize',16);
% ylabel('�Ӵ���/N','FontSize',16);
% ax = gca;%��ȡ��ǰ���������
% set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
% set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ
% ylim([0, 50]);
% 
% plot(Area_R,'LineWidth', 1);
% title('�Ӵ�����仯','FontSize',16);
% xlabel('֡��','FontSize',16);
% ylabel('���/cm^2','FontSize',16);
% ax = gca;%��ȡ��ǰ���������
% set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
% set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ

%% ��
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
legend( '�������','�����ڲ�','��һ�Ź�',...
    '�ڶ��Ź�','�����Ź�','�����Ź�','�����Ź�','�㹭���',...
    '�㹭�ڲ�', '��Ĵָ', '������ָ');
title('�ҽ�����Ӵ����仯','FontSize',16);
xlabel('���루ms��','FontSize',16);
ylabel('�Ӵ���/kg','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ
ylim([0, 35]);

%% ��
figure;%�ٽ���һ������

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
legend( '�������','�����ڲ�','��һ�Ź�',...
    '�ڶ��Ź�','�����Ź�','�����Ź�','�����Ź�','�㹭���',...
    '�㹭�ڲ�', '��Ĵָ', '������ָ');
title('�������Ӵ����仯','FontSize',16);
xlabel('���루ms��','FontSize',16);
ylabel('�Ӵ���/kg','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ
ylim([0, 35]);

%% ��ȡ���ߵ�һ���֣���R��
figure
x_min =1520;
x_max =2520;
% �ҳ�Ҫ��ȡ�����ݵ�����
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
% �����Ʊ���ȡ�����ݵ�
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

% �Զ��� x ��̶�ֵ�ͱ�ǩ
custom_x_values = [x_min:250: x_max];
custom_x_labels = {'0', '250', '500','750','1000'};
% ���� x ��̶�ֵ�ͱ�ǩ
xticks(custom_x_values);
xticklabels(custom_x_labels);
legend( '�������','�����ڲ�','��һ�Ź�',...
    '�ڶ��Ź�','�����Ź�','�����Ź�','�����Ź�','�㹭���',...
    '�㹭�ڲ�', '��Ĵָ', '������ָ');
title('�ҽ�����Ӵ����仯','FontSize',16);
xlabel('���루ms��','FontSize',16);
ylabel('�Ӵ���/kg','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ
ylim([0, 30]);

%% ��ȡ���ߵ�һ���֣���L��
figure
x_min=3880;
x_max=4670;
%�ҳ�Ҫ��ȡ�����ݵ�����
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
% �����Ʊ���ȡ�����ݵ�
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

% �Զ��� x ��̶�ֵ�ͱ�ǩ
custom_x_values = [x_min:197.5: x_max];
custom_x_labels = {'0', '197.5', '395','592.5','790'};
% ���� x ��̶�ֵ�ͱ�ǩ
xticks(custom_x_values);
xticklabels(custom_x_labels);
legend( '�������','�����ڲ�','��һ�Ź�',...
    '�ڶ��Ź�','�����Ź�','�����Ź�','�����Ź�','�㹭���',...
    '�㹭�ڲ�', '��Ĵָ', '������ָ');
title('�������Ӵ����仯','FontSize',16);
xlabel('���루ms��','FontSize',16);
ylabel('�Ӵ���/kg','FontSize',16);
ax = gca;%��ȡ��ǰ���������
set(ax.XAxis, 'LineWidth', 1);  % ����X��Ĵ�ϸ
set(ax.YAxis, 'LineWidth', 1);  % ����Y��Ĵ�ϸ
ylim([0,22]);