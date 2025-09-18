%求解评价指标的，10.12
clc
clear
close all
load qkpl.mat;
load kpl.mat;
load nsga.mat;
load vns.mat;
load gsa.mat;

QKPL=qkpl;
KPL=kpl;
NSGA=nsga;
VNS=vns;
GSA=gsa;
indices = fit(:, 1) < 300;
values = fit(indices, :);
disp(values);
format long
aa=[];bb=[];cc=[];dd=[];ee=[];
c=0;
for i=1:20
    aa=[aa;QKPL{1,i}];
    bb=[bb;KPL{1,i}];
    cc=[cc;NSGA{1,i}];
    dd=[dd;VNS{1,i}];
    ee=[ee;GSA{1,i}];
end

% aa(:,3) = aa(:,3)*-1;
% bb(:,3) = bb(:,3)*-1;
% cc(:,3) = cc(:,3)*-1;
% dd(:,3) = dd(:,3)*-1;
% ee(:,3) = ee(:,3)*-1;

format short
% digits(6);

aaaa1=[min(aa(:,1:3));mean(aa)]
aaaa2=[min(bb(:,1:3));mean(bb)]
aaaa3=[min(cc(:,1:3));mean(cc)]
aaaa4=[min(dd(:,1:3));mean(dd)]
aaaa5=[min(ee(:,1:3));mean(ee)]

%下面的代码是为了确保参数指标也是满足要求的即IGD，RPF，HV均符合要求
% 获取最大值和最小值
% 1.把数组转成矩阵
QKPL_P1_m=trans_to_m(QKPL);
KPL_P1_m=trans_to_m(KPL);
NSGA_P1_m=trans_to_m(NSGA);
VNS_P1_m=trans_to_m(VNS);
GSA_P1_m=trans_to_m(GSA);
aa=[1,2,3,4
    0,1,5,6];

% QKPL_P1_m(:,4) = QKPL_P1_m(:,4)*-1;
% KPL_P1_m(:,4) = KPL_P1_m(:,4)*-1;
% NSGA_P1_m(:,4) = NSGA_P1_m(:,4)*-1;
% VNS_P1_m(:,4) = VNS_P1_m(:,4)*-1;
% GSA_P1_m(:,4) = GSA_P1_m(:,4)*-1;

a=[];
a=[a;[ones(size(QKPL_P1_m,1),1),QKPL_P1_m]];
a=[a;[2*ones(size(KPL_P1_m,1),1),KPL_P1_m]];
a=[a;[3*ones(size(NSGA_P1_m,1),1),NSGA_P1_m]];
a=[a;[4*ones(size(VNS_P1_m,1),1),VNS_P1_m]];
a=[a;[5*ones(size(GSA_P1_m,1),1),GSA_P1_m]];
fit = a(:,3:5);
outpopulationpareto=a;
frontvalue = Newranking(fit);
outpopulationpareto(:,7)=frontvalue;
selected_rows = outpopulationpareto(outpopulationpareto(:, 7) == 1, :);
% 使用unique函数去除重复行
pareto = unique(selected_rows, 'rows', 'stable');
%归一化
[aa,~]=min(a(:,3:5));
[bb,~]=max(a(:,3:5));
outpopulation=a;
outpopulation(:,3)=(outpopulation(:,3)-aa(1))/(bb(1)-aa(1));
outpopulation(:,4)=(outpopulation(:,4)-aa(2))/(bb(2)-aa(2));
outpopulation(:,5)=(outpopulation(:,5)-aa(3))/(bb(3)-aa(3));


outpopulation_pareto_a=feizhipei(a);

%1.求非支配解
outpopulation_pareto=feizhipei(outpopulation);
a1=cell(1,20);
a2=cell(1,20);
a3=cell(1,20);
a4=cell(1,20);
a5=cell(1,20);

for i=1:20
    b1=find(outpopulation(:,1)==1);
    b2=find(outpopulation(:,2)==i);
    b3=intersect(b1,b2);
    a1{1,i}=outpopulation(b3,3:5);
end
for i=1:20
    b1=find(outpopulation(:,1)==2);
    b2=find(outpopulation(:,2)==i);
    b3=intersect(b1,b2);
    a2{1,i}=outpopulation(b3,3:5);
end
for i=1:20
    b1=find(outpopulation(:,1)==3);
    b2=find(outpopulation(:,2)==i);
    b3=intersect(b1,b2);
    a3{1,i}=outpopulation(b3,3:5);
end
for i=1:20
    b1=find(outpopulation(:,1)==4);
    b2=find(outpopulation(:,2)==i);
    b3=intersect(b1,b2);
    a4{1,i}=outpopulation(b3,3:5);
end
for i=1:20
    b1=find(outpopulation(:,1)==5);
    b2=find(outpopulation(:,2)==i);
    b3=intersect(b1,b2);
    a5{1,i}=outpopulation(b3,3:5);
end

%计算评价指标
zhibiao_last=zeros(20,10);
for f=1:20
    a_1=a1{1,f}(:,1:3);
    b_1=a2{1,f}(:,1:3);
    c_1=a3{1,f}(:,1:3);
    d_1=a4{1,f}(:,1:3);
    e_1=a5{1,f}(:,1:3);
    pop=[a_1;b_1;c_1;d_1;e_1];
    zhibiao1=zb(a_1,outpopulation_pareto(:,4:6));
    zhibiao2=zb(b_1,outpopulation_pareto(:,4:6));
    zhibiao3=zb(c_1,outpopulation_pareto(:,4:6));
    zhibiao4=zb(d_1,outpopulation_pareto(:,4:6));
    zhibiao5=zb(e_1,outpopulation_pareto(:,4:6));
    zhzhibiao=[zhibiao1,zhibiao2,zhibiao3,zhibiao4,zhibiao5];
    zhibiao_last(f,:)=zhzhibiao;
end

% for f=1:20
%     a_1=a1{1,f}(:,1:4);
%     b_1=a2{1,f}(:,1:4);
%     c_1=a3{1,f}(:,1:4);
%     d_1=a4{1,f}(:,1:4);
%     e_1=a5{1,f}(:,1:4);
%     pop=[a_1;b_1;c_1;d_1;e_1];
%     zhibiao1=calculate_igd(a_1,outpopulation_pareto(:,3:6));
%     zhibiao2=calculate_igd(b_1,outpopulation_pareto(:,3:6));
%     zhibiao3=calculate_igd(c_1,outpopulation_pareto(:,3:6));
%     zhibiao4=calculate_igd(d_1,outpopulation_pareto(:,3:6));
%     zhibiao5=calculate_igd(e_1,outpopulation_pareto(:,3:6));
%     zhzhibiao=[zhibiao1,zhibiao2,zhibiao3,zhibiao4,zhibiao5];
%     zhibiao_last(f,:)=zhzhibiao;
% end


% end
zhibiao_last1=zhibiao_last;

IGD_m=zhibiao_last1(:,[1,3,5,7,9]);
HV_m=zhibiao_last1(:,[2,4,6,8,10]);

Fddd=min(IGD_m);     %IGD最小值
Fddd1=mean(IGD_m);   %IGD平均值
Feee=max(HV_m);      %HV最大值
Feee1=mean(HV_m);    %HV平均值

zhibiao111=[Fddd;Fddd1;Feee;Feee1];

HV_3=HV_m;

[~,a]=min(min(IGD_m));
[~,a1]=min(mean(IGD_m));

[~,c]=max(max(HV_m));
[~,c1]=max(mean(HV_m));

IGD_m   
HV_m
%% 处理IGD
means_igd = mean(IGD_m )
mins_igd = min(IGD_m )


%% 处理HV
means_hv = mean(HV_m)
max_hv = max(HV_m)


[a,a1,c,c1]    %排名

[~,ww1]=sort(min(IGD_m));
[~,ww2]=sort(mean(IGD_m));
[~,ww3]=sort(max(HV_m));
[~,ww4]=sort(mean(HV_m));

format long   %kstest2
[h1,p1,k]=ranksum(HV_m(:,1),HV_m(:,2),0.05);
[h2,p2,k]=ranksum(HV_m(:,1),HV_m(:,3),0.05);
[h3,p3,k]=ranksum(HV_m(:,1),HV_m(:,4),0.05);
[h4,p4,k]=ranksum(HV_m(:,1),HV_m(:,5),0.05);
w9=[p1,p2,p3,p4];
m9=[h1,h2,h3,h4];
kHV=double([m9])

[h1,p1,k]=ranksum(IGD_m(:,1),IGD_m(:,2),0.05);
[h2,p2,k]=ranksum(IGD_m(:,1),IGD_m(:,3),0.05);
[h3,p3,k]=ranksum(IGD_m(:,1),IGD_m(:,4),0.05);
[h4,p4,k]=ranksum(IGD_m(:,1),IGD_m(:,5),0.05);
w9=[p1,p2,p3,p4];
z9=[h1,h2,h3,h4];
kIGD=double([z9])


figure 
zhibiao=zhibiao_last1;
IGD=zeros(20,5);
IGD(:,1)=zhibiao(:,1);
IGD(:,2)=zhibiao(:,3);
IGD(:,3)=zhibiao(:,5);
IGD(:,4)=zhibiao(:,7);
IGD(:,5)=zhibiao(:,9);

HV=zeros(20,5);
HV(:,1)=zhibiao(:,2);
HV(:,2)=zhibiao(:,4);
HV(:,3)=zhibiao(:,6);
HV(:,4)=zhibiao(:,8);
HV(:,5)=zhibiao(:,10);

subplot(1,2,1)
h=boxplot(IGD);
g = gca;
line=g.Children.Children;
line(1).Color ='b';
line(2).Color ='b';
line(3).Color ='b';
line(4).Color ='b'; %最下面的线
line(5).Color ='b';
line(6).Color ='r';
line(7).Color ='r';
line(8).Color ='r';
line(9).Color ='r';
line(10).Color ='r';
line(11).Color ='b'; %最下面的线
line(12).Color ='b';
line(13).Color ='b';
line(14).Color ='b';
line(15).Color ='b';
line(16).Color ='k';
line(17).Color ='k';
line(18).Color ='k'; %最下面的线
line(19).Color ='k';
line(20).Color ='k';
line(21).Color ='k';
set(h,'Linewidth',1.5');
set(gca,'XTickLabel',{'QDS-KOA','KOA','NSGAII','VNS-NSGAII','GSA'})
ylabel('IGD')
subplot(1,2,2)
n=boxplot(HV);
g2 = gca;
line1=g2.Children.Children;
line(1).Color ='b';
line(2).Color ='b';
line(3).Color ='b';
line(4).Color ='b'; %最下面的线
line(5).Color ='b';
line(6).Color ='r';
line(7).Color ='r';
line(8).Color ='r';
line(9).Color ='r';
line(10).Color ='r';
line(11).Color ='b'; %最下面的线
line(12).Color ='b';
line(13).Color ='b';
line(14).Color ='b';
line(15).Color ='b';
line(16).Color ='k';
line(17).Color ='k';
line(18).Color ='k'; %最下面的线
line(19).Color ='k';
line(20).Color ='k';
line(21).Color ='k';
set(n,'Linewidth',1.5');
set(gca,'XTickLabel',{'QDS-KOA','KOA','NSGAII','VNS-NSGAII','GSA'})
% set(gca,'YTickLabel','S');
ylabel('HV')




%% 计算指标
function  zhibiao=zb(nowpop,obj)
pop=obj;
%求反世代距离
dist=zeros(1,size(pop,1));
for i=1:size(pop,1)
    dist1=zeros(1,size(nowpop,1));
    [a1,~]=min(pop(:,1:3));
    [b1,~]=max(pop(:,1:3));
    for j=1:size(nowpop,1)
        f1=(pop(i,1)-nowpop(j,1));
        f2=(pop(i,2)-nowpop(j,2));
        f3=(pop(i,3)-nowpop(j,3));
        f4=f1*f1+f2*f2+f3*f3;
        dist1(j)=sqrt(f4);
    end
    [a3,~]=min(dist1);
    dist(i)=a3;
end
IGD=sum(dist)/size(pop,1);

%计算HIV，参考点(0,0,0)
%就是坐标点围城的体积-------------参考点（1，1，1, 1）
%不是简单的体积相加，是围成的立方体的面积
HV_sum=0;
%先对x坐标进行升序
nowpop=[nowpop;[1,1,1]];
[a,b]=sort(nowpop(:,1));
%a是排序后的顺序，b是对应的位置
for i=1:(size(nowpop,1)-1)
    %c=(nowpop(b(i),1)-nowpop(b(i+1),1))*(1-nowpop(b(i),2));
    % 计算相邻点在第一个维度上的距离乘积，这里假设其他维度都是固定的，例如第二维度是 y，第三维度是 z，第四维度是 w
    c = (nowpop(b(i),1)-nowpop(b(i+1),1))*(1-nowpop(b(i),2))*(1-nowpop(b(i),3));      
    HV_sum=HV_sum+abs(c);
end
HV=HV_sum;

zhibiao=[IGD,HV];
end

%% 筛选得到非劣解
function outpopulation=feizhipei(population)
aaa=1:size(population,1);
lastpopulation=[aaa',population];
popuvalue=population(:,3:4);
% population=population(:,3:end);
% fzpweizhi=[];   %%用来存放非支配解位置的
fzpweizhi_p=zeros(1,size(population,1));
for i=1:size(population,1)
    fzplabel=0;
    popuvalue1=popuvalue;
    for j=1:size(population,1)
        if i==j
            continue
        end
        chazhipopulation=popuvalue1(i,:)-popuvalue1(j,:);
        if chazhipopulation(1)<0.000000001&&chazhipopulation(1)>-0.000000001
            chazhipopulation(1)=0;
        end
        if chazhipopulation(2)<0.000000001&&chazhipopulation(2)>-0.000000001
            chazhipopulation(2)=0;
        end
        fzplabel2=length(find(chazhipopulation>=0));%%看是否该个体被支配
        fzplabel3=length(find(chazhipopulation>0));
        fzplabel4=length(find(chazhipopulation==0));
        if fzplabel2==2&&(fzplabel3>=1)
            fzplabel=1;
            break
        end
        if fzplabel4==2&&(i<j)  %%%%此条件是为了删除相同的值
            fzplabel=1;
            break
        end
    end
    if fzplabel==0
        %            fzpweizhi=[fzpweizhi;i];
        fzpweizhi_p(1,i)=1;
    end
    %     fzplabel=0;
end
fzpweizhi_last_p=find(fzpweizhi_p==1)';
outpopulation=lastpopulation(fzpweizhi_last_p,:);
end

%% 转成矩阵
function NSGA2VNS_data_last_m=trans_to_m(NSGA2VNS_data_last)
data_m=[];
for i=1:20
    a=NSGA2VNS_data_last{1,i};
    b=NSGA2VNS_data_last{1,i}(:,3);
    c=find(b>0);
    d=a(c,:);
    NSGA2VNS_data_last{1,i}=d;
    a=size(NSGA2VNS_data_last{1,i}(:,1:3),1);
    data_m=[data_m;[i*ones(a,1),NSGA2VNS_data_last{1,i}(:,1:3)]];
end
NSGA2VNS_data_last_m=data_m;
end

function igd_value = calculate_igd( pareto_set,pareto_frontier)

    % 计算真实 Pareto 前沿集合的大小
    pareto_frontier_size = size(pareto_frontier, 1);
    
    % 计算算法生成的 Pareto 解集的大小
    pareto_set_size = size(pareto_set, 1);
    
    % 初始化 IGD 值为 0
    igd_value = 0.0;
    
    % 遍历算法生成的 Pareto 解集中的每个解
    for i = 1:pareto_set_size
        min_distance = Inf; % 初始化最小距离为正无穷大
        % 计算当前解与真实 Pareto 前沿集合中每个解之间的距离
        for j = 1:pareto_frontier_size
            % 计算 Euclidean 距离
            distance = norm(pareto_set(i, :) - pareto_frontier(j, :));
            % 更新最小距离
            if distance < min_distance
                min_distance = distance;
            end
        end
        % 将最小距离加入 IGD 值
        igd_value = igd_value + min_distance;
    end
    
    % 计算平均最小距离
    igd_value = igd_value / pareto_set_size;
end

%% frontvalue = Newranking(ObjV)                         %非支配排序
function frontvalue = Newranking(ObjV)                         %非支配排序

fnum=0;                                                        %当前分配的前沿面编号
cz=false(1,size(ObjV,1));                              %记录个体是否已被分配编号
frontvalue=zeros(size(cz));                            %每个个体的前沿面编号
[functionvalue_sorted,newsite]=sortrows(ObjV);         %对种群按第一维目标值大小进行升序排序
while ~all(cz)                                         %开始迭代判断每个个体的前沿面,采用改进的deductive sort
    fnum=fnum+1;
    d=cz;
    for i=1:size(ObjV,1)
        if ~d(i)
            for j=i+1:size(ObjV,1) %判断i对应的所有集合里面的支配和非支配的解，被i支配则为1，不被i支配则为0
                if ~d(j)
                    k=1;
                    for m=2:size(ObjV,2) %判断是否支配，找到个体p不支配的个体，标记为k=0
                        if functionvalue_sorted(i,m)>functionvalue_sorted(j,m)%判断i,j是否支配，如果成立i,j互不支配
                            k=0;%i,j互不支配
                            break
                        end
                    end
                    if k
                        d(j)=true;%那么p所支配的个体k=1并记录在d里，则后面该个体已被支配就不能在这一层里进行判断
                    end
                end
            end
            frontvalue(newsite(i))=fnum;%实际位置的非支配层赋值
            cz(i)=true;
        end
    end
end
end