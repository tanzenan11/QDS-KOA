clc;
clear all;
% 三目标：平滑程度，利润，能耗
P = [0,0,0,0,0,0,0,0,0,0;
    -1,0,0,0,0,0,0,-1,-1,-1;
    -1,0,0,0,0,0,0,-1,-1,-1;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,1,1,0,0,0,0;
    0,0,0,1,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;];
%需求属性
%需求属性
DA = [0,1,0,0,0,0,1,0,1,0,];
%复杂属性
CA = [0,0,0,0,0,0,0,0,1,0];
%危险属性
HA = [0,0,0,0,0,0,1,0,0,0];
% 拆卸时间
Td = [14;10;12;18;23;16;20;36;14;10];
%预计节拍时间
C=50;
%任务利润
[n,a] = size(P);
Pr = round(20 + (100 - 20) * rand(1, a),2);
Pr = transpose(Pr);
%开启一个工作站的成本
Cw = 0.5;
%购买机器人的固定成本
Cr = 2;
%使用工人的单位时间固定成本
Cpt = 6/3600;
%使用机器人分解的单位时间成本
Crt = 3/3600;
%处理危险任务的额外成本
Ch = 0.01;
%使用工人的固定的单位时间任务能耗
Ept = 0.02;
%使用机器人分解的单位时间任务能耗
Ert = 0.03;
%处理危险任务额外能耗
Eh = 0.01;
% 规定参数
tic %开始计时
np = 50;  %行星数
Tmax = 50; %最大迭代次数
dim = 3; %目标数
[n,A] = size(P); % n代表零件数
favg = 0;   %存储平均适应度值
Q = zeros(20,10);
pareto_frontier = [0 0 0];  %计算IGD值真实点
ub = 1;
lb = 0;
IGD_one_sum = 0;    %第一代个体的IGD值求和
IGD_one_max = 0;    %第一代个体的IGD值最大值
IGD_one_S2 = 0;
Gamma = 0.6;%Q-learning远见系数
Alpha = 0.4;%Q-learning学习系数
%Q-learning 领域搜索比例
Proportion=[
    0.3, 0.1, 0.4, 0.2
    0.2, 0.3, 0.2, 0.3
    0.1, 0.4, 0.3, 0.2
    0.3, 0.2, 0.2, 0.3
    0.2, 0.2, 0.3, 0.3
    0.1, 0.3, 0.3, 0.3
    0.4, 0.1, 0.3, 0.2
    0.2, 0.4, 0.1, 0.3
    0.3, 0.3, 0.1, 0.3
    0.1, 0.2, 0.4, 0.3];
%% 规定参数
tic %开始计时
np = 70;  %种群数
pc = 0.8;  %交叉率
pm = 0.25;  %变异率
[n,A] = size(P);


[Chrom,A1] = Shengcheng(np,n,P,DA,HA);
a_time=toc;   %记录时间
time = 10;
gen = 1;
MaxGen = 10;
while gen < MaxGen
    % %计算目标值
    [Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值
    [F,Level] = non_domination_sort(np,fit);  %非支配排序
    Chrom(:,2*n+1) = Level;
    Chrom(:,2*n+2) = fit(:,1);
    Chrom(:,2*n+3) = fit(:,2);
    Chrom(:,2*n+4) = fit(:,3);
    Chrom = crowding_distance_sort(F,Chrom,n,fit);  %根据适应度值拥挤距离进行选择
    Chrom_son_jc = Jiaocha(Chrom,n,pc);      %交叉操作
    Chrom_son_jc = aberranceJm(Chrom,np,n,pm,P);
    Chrom_son_by = aberranceJm(Chrom_son_jc,np,n,pm,P);   %变异操作
    %将两代个体合并
    combine_Chrom = zeros(2*np,2*n);
    combine_Chrom(1:np,:) = Chrom(1:np,1:2*n);
    combine_Chrom(np+1:2*np,:) = Chrom_son_by(1:np,1:2*n);
    [Chrom,fit,PL_Fit] = Cal(combine_Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值;
    [F,Level] = non_domination_sort(2*np,fit);
    combine_Chrom(:,2*n+1) = Level;
    combine_Chrom(:,2*n+2) = fit(:,1);
    combine_Chrom(:,2*n+3) = fit(:,2);
    combine_Chrom(:,2*n+4) = fit(:,3);
    combine_Chrom = crowding_distance_sort(F,combine_Chrom,n,fit);  %根据适应度值拥挤距离进行选择
    %精英保留产生下一代种群
    Chrom = elitism(np,combine_Chrom,n,fit);
    a_time=toc;   %记录时间
    %计算这一代的IGD和以及IGD最大值
    IGD_t_sum = 0;
    IGD_t_max = 0;
    %把第t代个体的IGD值求和
    for j = 1:np
        IGD_t_sum = IGD_t_sum + calculate_igd(pareto_frontier,fit(j,:));
        if calculate_igd(pareto_frontier,fit(j,:)) > IGD_t_max
            IGD_t_max = calculate_igd(pareto_frontier,fit(j,:));   %求出第t代个体的IGD最大值
        end
    end
    s1 = IGD_t_sum/IGD_one_sum;
    % 将种群多样性进行归一化
    IGD_t_S2 = 0;
    for j = 1:np
        IGD_t_S2 = IGD_t_S2 + abs(calculate_igd(pareto_frontier,fit(j,:)) - IGD_t_sum/np);
    end
    s2 = IGD_t_S2/IGD_one_S2;     %种群多样性进行归一化结果
    %将最佳适应度归一化
    s3 = IGD_t_max / IGD_one_max;
    s = 0.35*s1 + 0.35*s2 + 0.3*s3;
    [pc,pm,Q,randomIndex,status] = Q_learning(MaxGen,s,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max,gen,Q,Proportion);
    gen = gen + 1;
end
%%%帕累托解
L=1;
para=[];  % 帕累托解
for i=1:np
    if  Chrom(i,2*n+1)==1
        para(L,:)=Chrom(i,:);
        L=L+1;
    end
end
% 去除重复个体
[unique_para, idx, ~] = unique(para, 'rows', 'stable');
disp(unique_para);
para = unique_para(:,n*2+2:n*2+4);
Foutput = para;

%% 生成种群
function [Chrom,A1] = Shengcheng(np,n,P,DA,HA)
Chrom = zeros(np,n*2);
% 生成操作顺序
for i = 1:np
    P0 = P;
    bused = [];   %存放已选择的任务
    for j = 1:n
        A1 = [];
        A2 = [];
        % 生成or集合
        a = 1;
        for k = 1:n   %遍历行
            if ismember(-1,P0(k,:)) == 1 %存在or约束,记录行数（前驱零件）
                A1(1,a) = k;
                a = a + 1;
            end
        end
        for k = 1:n   %遍历列
            if ismember(-1,P0(:,k)) == 1 %存在or约束,记录列数（前驱零件）
                A1(1,a) = k;
                a = a + 1;
            end
        end
        % 遍历列生成优先任务集合
        c = 1;
        for k = 1:n
            if any(P0(:,k)) == 0     %k零件没有任何前驱约束
                A2(1,c) = k;
                c = c + 1;
            end
        end
        %保证选出来的任务没有被安排过
        A2 = setdiff(A2,bused);
        [~,n1] = size(A2);
        val = A2(1,randi(n1)); %在可以加工的任务中随机选取选取一个任务
        p1 = [];
        if ismember(val,A1) == 1 %若所选任务存在or关系
            c = find(P0(val,:) == -1); % 找到val or约束的位置
            [~,n2] = size(c);         % or约束对象个数
            %消除or关系影响
            for k = 1:n2
                p1 = P0(:,c(1,k));
                p1(p1 == -1) = 0;
                P0(:,c(1,k)) = p1;
            end
            %消除and关系影响
            p1 = P0(val,:);  %将val对其后继零件影响消除
            p1 = 0;
            P0(val,:) = p1;
        else
            %消除and关系影响
            p1 = P0(val,:);  %将val对其后继零件影响消除
            p1 = 0;
            P0(val,:) = p1;
        end
        Chrom(i,j) = val;
        % 判断选出的工件属性是

        bused(1,j) = val;   %将已选择的任务进行保存，保证下次不会选中
    end
end
P0 = P;
indices_ca = find(DA == 1);  % 需求操作
indices_ha = find(HA == 1);  % 危险操作
indices = [indices_ca,indices_ha];    %所有有需求和危险的操作
len = length(indices);
A = indices;
dis = [];    %拆卸集合
while ~isempty(A2)
    len = length(indices);
    for i = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
        temp = indices(i);   %操作
        A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
        A1 = A1';
        A = [A,A1];
    end
    A = unique(A);   %必须拆卸的步骤
    dis = [dis,A];
    %获得前置步骤的前置步骤
    A2 = A(~ismember(A, indices));
    indices = A2;
    A = indices;
end
dis = unique(dis, 'rows', 'stable');   %必须拆卸的步骤
for j = 1:np
    notdis = setdiff(1:n, dis);               %其余步骤
    len = length(notdis);
    if len>0
        random_number = randi([1, len]); % 生成 0 到列表长度的随机数
        indices = randperm(len, random_number); % 生成 random_number 个不重复的随机索引
        random_dis = notdis(indices); % 随机取出的拆卸操作
        B = random_dis;               %重新存储一次随机取出的拆卸操作
        random_dis_fm = random_dis;
        A2 = [0];    %保证进入循环
        while ~isempty(A2)
            len = length(random_dis);
            for i = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                temp = random_dis(i);   %操作
                A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
                A1 = A1';
                A = [A,A1];
            end
            A = unique(A);   %必须拆卸的步骤
            random_dis_fm = [random_dis_fm,A];
            %获得前置步骤的前置步骤
            A2 = A(~ismember(A, random_dis));
            random_dis = A2;
            A = random_dis;
        end
        random_dis_fm = unique(random_dis_fm, 'rows', 'stable');   %必须拆卸的步骤
        ChromOne = Chrom(j,:);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
        [~, idx2] = ismember(random_dis_fm, ChromOne);  %随机的拆卸步骤
        Chrom(j,n+idx2) = 1;
    else
        ChromOne = Chrom(j,:);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
    end
end
end


%% 计算目标值
function [Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh)
[np,~] = size(Chrom);
fit = zeros(np,3);          %目标值

for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    Temp = Chromone(indices);   %需要拆卸的步骤
    len = length(Temp);
    Profit =  sum(Pr(Temp));   %拆卸得到的利润
    %对需要拆卸的步骤进行解码
    g = 1;           %工厂号
    Mt = zeros(100,1);  %工厂时间
    Mr = 0;   %人工站数目
    Mm = 0;   %机器站数目
    cost = 0; %成本
    profit = 0;   %利润
    ec = 0;       %能耗
    MM=zeros(15,15);
    MT=zeros(15,15);
    MP=zeros(15,1);
    for j = 1:len
        a = Temp(j);   %拆卸步骤
        da = DA(a);    %需求属性
        ca = CA(a);    %复杂属性
        ha = HA(a);    %危害属性
        t = Td(a);      %拆卸时间
        if ca~= 1 && ha ~= 1     %无属性拆卸给机器工厂
            if j ==1
                p = 1;
                h=1;
                Mm = Mm + 1;%机器工厂数加1
            end
            newP = 1;
            if newP ~= p                 %和前一个工厂属性不同
                g = g + 1;  %新开一个工厂
                h=1;
                Mm = Mm + 1;  %机器工厂数加1
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t;   %计算成本
                ec = ec + Ert*t;       %计算能耗
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t;
                ec = ec + Ert*t;       %计算能耗
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Crt*t;
                    ec = ec - Ert*t;
                    g = g + 1;
                    h=1;
                    Mm = Mm + 1;  %机器工厂数加1
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Crt*t;
                    ec = ec + Ert*t;       %计算能耗
                end
            end
        elseif ca == 1 && ha ~= 1    %复杂属性工人工厂
            if j == 1
                p = -1;
                h=1;
                Mr = Mr + 1;  %工人工厂数加1
            end
            newP = -1;       %工厂属性
            if newP ~= p                %和前一个工厂属性不同
                g = g + 1;  %新开一个工厂
                h=1;
                Mr = Mr + 1; %工人工厂数加1
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Cpt*t;
                ec = ec + Ept*t;       %计算能耗
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Cpt*t;
                ec = ec + Ept*t;       %计算能耗
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Cpt*t;
                    ec = ec - Ept*t;       %计算能耗
                    g = g + 1;
                    h=1;
                    Mr = Mr + 1;  %工人厂数加1
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Cpt*t;
                    ec = ec + Ept*t;       %计算能耗
                end
            end
        else             %危险属性机器人工厂，注意如果既复杂又危险得任务同样属于危险任务
            if j == 1
                p = 1;
                h=1;
                Mm = Mm + 1;  %机器工厂数加1
            end
            newP = 1;       %工厂属性
            if newP ~= p                   %和前一个工厂不同
                g = g + 1;  %新开一个工厂
                h=1;
                Mm = Mm + 1;  %机器工厂数加1
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Crt*t - Ch*t;  %注意加上处理危险任务的额外成本
                    ec = ec - Ept*t - Eh*t;       %计算加上处理危险任务的额外能耗
                    g = g + 1;
                    h=1;
                    Mm = Mm + 1;  %机器工厂数加1
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                    ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
                end
            end
        end
        p = newP;          %存下本次工厂属性
        %Chrom(i,2*n+j) = p;
        MM(g,h)=Temp(j);
        MT(g,h)=Td(Temp(j));
        MP(g,1)=p;
        h=h+1;
    end
    Mt = Mt(1:g,:);
    % 计算负载均衡
    F = 0;
    Ct = max(Mt);
    EC=0;
    for j = 1:g
        F = F + (Ct - Mt(j,1))^2;
        % if MP(j,1) == 1
        %     EC= EC+ (C - Mt(j,1))*0.8;
        % else
        %     EC= EC+ (C - Mt(j,1))*0.2;
        % end
    end
    %计算利润
    profit = Profit - cost - g * Cw - Mm * Cr;
    fit(i,1) = sqrt(F);
    fit(i,2) = -profit;
    fit(i,3) = ec+EC;
    PL_Fit(i,1) = sum(fit(i,:))/3;   %加权平均值
end
end

%% 非支配排序
function [F,Level] = non_domination_sort(np,fit)
%non_domination_sort 初始种群的非支配排序
%初始化pareto等级为1
pareto_rank = 1;
F(pareto_rank).ss = [];%pareto等级为pareto_rank的集合,F存储的是各等级解的索引
p = [];%每个个体p的集合
[~,a] = size(fit);  % a为目标数
for i = 1:np
    %%%计算出种群中每个个体p的被支配个数n和该个体支配的解的集合s
    p(i).n=0;%被支配个体数目n
    p(i).s=[];%支配解的集合s
    for j = 1:np
        less=0; %目前个体的目标函数值小于所有个体的目标函数值数目
        equal=0; %目前个体的目标函数值等于所有个体的目标函数值数目
        greater=0; %目前个体的目标函数值大于所有个体的目标函数值数目
        for k = 1:a   %k==1是为f1,k==2为f2
            if(fit(i,k) < fit(j,k)) %代表支配别人
                less = less + 1;
            elseif(fit(i,k) == fit(j,k))
                equal = equal +1;
            else %代表被别人支配
                greater = greater + 1;
            end
        end
        if(less == 0 && equal ~= a) %equal ~= 2 两个值都相同
            p(i).n = p(i).n + 1; %被支配个体数目+1
        elseif(greater==0 && equal ~= a)
            p(i).s = [p(i).s j]; %支配解
        end
    end
    %%%将种群中p(i).n = 1（非支配解）的个体放入集合F(1)中
    if(p(i).n == 0)
        Level(i,1) = 1;%储存个体的等级信息
        F(pareto_rank).ss = [F(pareto_rank).ss i];
    end
end
%%%求pareto等级为pareto_rank+1的个体
while ~isempty(F(pareto_rank).ss)
    temp = [];
    for i = 1:length(F(pareto_rank).ss)
        if ~isempty(p(F(pareto_rank).ss(i)).s)%若该pareto_rank的当前个体支配的解集不为空
            for j = 1:length(p(F(pareto_rank).ss(i)).s)
                p(p(F(pareto_rank).ss(i)).s(j)).n = p(p(F(pareto_rank).ss(i)).s(j)).n - 1;%nl=nl-1
                if p(p(F(pareto_rank).ss(i)).s(j)).n == 0
                    Level(p(F(pareto_rank).ss(i)).s(j),1) = pareto_rank+1;%储存个体的等级信息
                    temp = [temp p(F(pareto_rank).ss(i)).s(j)];
                end
            end
        end
    end
    pareto_rank = pareto_rank+1;
    F(pareto_rank).ss = temp;
end
end

%% 根据适应度值拥挤距离进行选择
function Chrom = crowding_distance_sort(F,Chrom,n0,fit)
%计算拥挤度
%%%按照pareto等级对种群中的个体进行排序
[np,~] = size(Chrom);
[~,index] = sort(Chrom(:,2*n0+1));
[~,mm1] = size(Chrom);
temp = zeros(np,mm1);
[~,a] = size(fit);     %a为目标数
for i = 1:length(index)%=np
    temp(i,:) = Chrom(index(i),:);%按照pareto等级排序后种群
end
%%%对于每个等级的个体开始计算拥挤度
current_index = 0;
for pareto_rank = 1:(length(F)-1)%计算F的循环时多了一次空，所以减掉
    %%拥挤度初始化为0
    nd = [];
    nd(:,1) = zeros(length(F(pareto_rank).ss),1);
    %y=[];%储存当前处理的等级的个体;u储存视图组合
    [~,mm2] = size(temp);
    y = zeros(length(F(pareto_rank).ss),mm2);%储存当前处理的等级的个体
    previous_index = current_index + 1;
    for i = 1:length(F(pareto_rank).ss)
        y(i,:) = temp(current_index + i,:);
    end
    current_index = current_index + i;
    %%对于每一个目标函数fm
    for i = 1:a
        %%根据该目标函数值对该等级的个体进行排序
        [~,index_objective] = sort(y(:,2*n0+1+i));
        %objective_sort=[];%通过目标函数排序后的个体
        [~,mm3] = size(y);
        objective_sort = zeros(length(index_objective),mm3);%通过目标函数排序后的个体
        objective_x = {};
        for j = 1:length(index_objective)
            objective_sort(j,:) = y(index_objective(j),:);
        end
        %%记fmax为最大值，fmin为最小值
        fmin = objective_sort(1,2 * n0 + 1 + i);
        fmax = objective_sort(length(index_objective),2 * n0 + 1 + i);
        %%对排序后的两个边界拥挤度设为1d和nd设为无穷
        y(index_objective(1),2*n0+a+1+i) = Inf;
        y(index_objective(length(index_objective)),2*n0+a+1+i) = Inf;
        %%计算nd = nd+(fm(i+1)-fm(i-1))/(fmax-fmin)
        for j = 2:(length(index_objective)-1)
            pre_f = objective_sort(j-1,2*n0 + 1 + i); %目标函数值
            next_f = objective_sort(j+1,2*n0 + 1 + i); %目标函数值
            if (fmax - fmin == 0)
                y(index_objective(j),2*n0 + a + 1 + i) = Inf;
            else
                y(index_objective(j),2*n0 + a + 1 + i) = (next_f - pre_f)/(fmax - fmin);
            end
        end
    end
    %多个目标函数拥挤度求和
    for i = 1:a
        nd(:,1) = nd(:,1)+y(:,2*n0+a+1+i);
    end
    %第2列保存拥挤度，其他的覆盖掉
    y(:,2*n0+1+a+1) = nd;
    y = y(:,1:(2*n0+1+a+1));
    temp_two(previous_index:current_index,:) = y;
end
Chrom = temp_two;
end

%% 交叉
function Chrom_son_jc = Jiaocha(Chrom,n,pc)
%两点映射交叉
% 随机选出两个父代，利用二元锦标赛策略选出最优父代
% 随机生成两段交叉片段，1号个体选择出的基因片段根据2号个体这几个基因中出现顺序进行改变
%举例 个体为1 3 2 5 6 8 7 4
%           1 2 5 3 6 8 7 4
% 选择第一个片段为2-3 一号个体基因为3 2 ，而这两个基因在二号个体中出现顺序为先2后3
% 所以最后一号个体为 1 2 3 6 8 7 4 5
%  注意点：因为二号个体中存在先2后3的顺序，说明2可以在3之间进行操作，所以不需考虑合理性
[np,index] = size(Chrom);
Chrom_son_jc = Chrom(:,1:n);
%找到每个个体需要拆卸的步骤保存起来
Temp = zeros(np,n);
for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    len = length(indices);
    Temp(i,1:len) = Chromone(indices);   %需要拆卸的步骤
end
for i = 1:np
    if rand < pc
        a = unidrnd(np);
        while a == i  %取出的另一个体不能是当前个体
            a = unidrnd(np);
        end
        p1(1,:) = Chrom_son_jc(i,:);
        p1(2,:) = Chrom_son_jc(a,:);
        %随机生成四个交叉点，从而生成我们需要的两段交叉片段
        p = randperm(n);
        p = p(1:2);
        p = sort(p);
        b=p(1,1);c=p(1,2);
        n1(1,1) = c - b + 1;   %基因片段基因数
        p01 = p1(1,b:c);
        %交叉片段一在对面顺序提取
        for k = 1:n1(1,1)
            w = 1;
            g = 1;
            while g == 1
                if p1(2,w)==p01(1,k)
                    p01(2,k) = w;  %在第二段中的位置
                    g = 0;
                else
                    w = w + 1;
                end
            end
        end
        p001 = [];
        [~,index] = sort(p01(2,:));
        for k = 1:length(index)
            p001(1,k) = p01(1,index(k)); %按照另一个体中的基因顺序排列
        end
        p1(1,b:c) = p001;
        Chrom_son_jc(i,1:n) = p1(1,:);
    end
end
for i = 1:np
    Chrom_son_jc(i,n+1:n*2) = 0;
end
%根据前面存储的Temp，重新置0，1
for i = 1:np
    ChromOne = Chrom_son_jc(i,:);
    % 获取第一行不等于0的元素的逻辑索引
    nonZeroIndices = Temp(i, :) ~= 0;
    % 使用逻辑索引获取需要拆卸的步骤
    dis = Temp(i, nonZeroIndices);
    [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
    Chrom_son_jc(i,n+idx1) = 1;
end
end

%% 变异
function Chrom_son_by = aberranceJm(Chrom,np,n,pm,P)
% 基因突变
% 随机取出变异点基因，得到其前置任务集和后置任务集,得到此基因可移动位置，再从可移动位置中随机选择一点进行插入操作
% 举例：1 2 5 7 9 10 8 3 6 4
% 随机突变基因为7，紧接在前的任务集是{3，4，5，6}，其紧接在后的任务是任务8。所以7，9，10这一段就是它可以插入的位置
% 若选择插在9之后得到的个体为1 2 5 9 7 10 8 3 6 4
Chrom_use = Chrom;
P0 = P;
%找到每个个体需要拆卸的步骤保存起来
Temp = zeros(np,n);
for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    len = length(indices);
    Temp(i,1:len) = Chromone(indices);   %需要拆卸的步骤
end

for i = 1:np
    if rand < pm
        Chrom_son = Chrom_use(i,1:n);
        len = length(Chrom_son);
        index = unidrnd(len);
        temp = Chrom_son(index);   % 取出变异点基因
        A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
        A2 = find(P0(temp, :) == -1 | P0(temp, :) == 1);  %后置与或集合
        % 找到A1中数字在Chrom_son中的索引位置
        indices = ismember(Chrom_son, A1);
        indexes_A1 = find(indices);
        % 找到A2中数字在Chrom_son中的索引位置
        indices = ismember(Chrom_son, A2);
        indexes_A2 = find(indices);
        if isempty(A1)
            if ~isempty(A2)                  %无前置有后置
                later_indices = indexes_A2(indexes_A2 > index);   %找到后置索引
                later_indices = min(later_indices);
                % 1到later_indices-1为插入区间
                random_number = randi([1, later_indices - 1]);  %插入位置
                if random_number < index
                    Chrom_son(random_number+1:index) = Chrom_son(random_number:index-1);
                    Chrom_son(random_number) = temp;
                else
                    Chrom_son(index:random_number-1) = Chrom_son(index+1:random_number);
                    Chrom_son(random_number) = temp;
                end
            end
        elseif isempty(A2)
            if ~isempty(A1)                  %无后置有前置
                former_indices = indexes_A1(indexes_A1 < index); %找到前置索引
                former_indices = max(former_indices);
                % former_indices+1到len是插入区间
                random_number = randi([former_indices + 1, len]); %插入位置
                if random_number < index
                    Chrom_son(random_number+1:index) = Chrom_son(random_number:index-1);
                    Chrom_son(random_number) = temp;
                else
                    Chrom_son(index:random_number-1) = Chrom_son(index+1:random_number);
                    Chrom_son(random_number) = temp;
                end
            end
        else
            former_indices = indexes_A1(indexes_A1 < index);
            later_indices = indexes_A2(indexes_A2 > index);
            former_indices = max(former_indices);
            later_indices = min(later_indices);
            % former_indices+1到later_indices-1是插入区间
            if any(former_indices ~= 0) && any(later_indices ~= 0)
                random_number = randi([former_indices + 1, later_indices - 1]); %插入位置
                if random_number < index
                    Chrom_son(random_number+1:index) = Chrom_son(random_number:index-1);
                    Chrom_son(random_number) = temp;
                else
                    Chrom_son(index:random_number-1) = Chrom_son(index+1:random_number);
                    Chrom_son(random_number) = temp;
                end
            end
        end
        Chrom_use(i,1:n) = Chrom_son;
    end
end
%根据前面存储的Temp，重新置0，1
for i = 1:np
    ChromOne = Chrom_use(i,:);
    % 获取第一行不等于0的元素的逻辑索引
    nonZeroIndices = Temp(i, :) ~= 0;
    % 使用逻辑索引获取需要拆卸的步骤
    dis = Temp(i, nonZeroIndices);
    [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
    Chrom_use(i,n+idx1) = 1;
end
Chrom_son_by = Chrom_use;
end

%% 精英保留策略
function Chrom=elitism(np,combine_Chrom,n,fit)
[~,a] = size(fit); %a为目标数
Chrom = zeros(np,(2 * n + a + 1 + 1));
% Chrom_rank = zeros(pop2,temp);
% %根据pareto等级从高到低进行排序
% [a,index_rank] = sort(combine_Chrom(:,2*n+1));
% for i = 1:pop2
%     Chrom_rank(i,:) = combine_Chrom(index_rank(i),:);
% end
Chrom_rank = combine_Chrom;
%求出最高的Pareto等级
max_rank = max(Chrom_rank(:,2*n+1));
%根据排序后的顺序，将等级相同的种群整个放入父代种群中，直到某一层不能全部放下为止
prev_index=0;   %记录已存放种群数
for i = 1:max_rank
    %寻找当前等级i个体里的最大索引
    current_index = max(find(Chrom_rank(:,(2*n+1)) == i));
    %不能放下为止
    if(current_index > np)  %放不下了
        %剩余种群数
        remain_pop = np - prev_index;
        %记录下等级为i的个体
        temp = Chrom_rank((prev_index + 1:current_index),:);
        %根据拥挤度从大到小排序
        [~,index_crowd] = sort(temp(:,n*2+4),'descend');
        %填满父代
        for j = 1:remain_pop
            Chrom(prev_index + j,:) = temp(index_crowd(j),:);
        end
        return;
    elseif(current_index < np)   % 还能放下
        Chrom((prev_index + 1:current_index),:) = Chrom_rank((prev_index + 1:current_index),:);
    else    %正好放满退出循环
        Chrom((prev_index + 1:current_index),:) = Chrom_rank((prev_index + 1:current_index),:);
        return;
    end
    prev_index = current_index;

end
end

%% Qlearning过程
function [Pc,Pm,Q,randomIndex,i] = Q_learning(MaxGen,s,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max,iter,Q,Proportion)
Epsilon = 0.9 - 0.85*iter/(MaxGen*10);
[state,action] = size(Q);
if s == 0
    i = 1; %表示此时状态
    if Epsilon < rand
        % 找到第一行的最大值及其索引
        [maxValue, maxIndex] = max(Q(1, :));
        % 找到所有最大值的索引
        maxIndices = find(Q(1, :) == maxValue);
        % 如果有多个最大值，随机选择一个索引
        randomIndex = maxIndices(randi(length(maxIndices)));   %动作指引
        selectedRow = Proportion(randomIndex,:);
        % 分别赋值给 Pc 和 Pm
        Pc = selectedRow(1);
        Pm = selectedRow(2);
    else
        randomIndex = unidrnd(action);    %动作指引
        selectedRow = Proportion(randomIndex,:);
        % 分别赋值给 Pc 和 Pm
        Pc = selectedRow(1);
        Pm = selectedRow(2);
    end

elseif s > 1
    i = 20;
    if Epsilon < rand
        % 找到第一行的最大值及其索引
        [maxValue, maxIndex] = max(Q(1, :));
        % 找到所有最大值的索引
        maxIndices = find(Q(1, :) == maxValue);
        % 如果有多个最大值，随机选择一个索引
        randomIndex = maxIndices(randi(length(maxIndices)));   %动作指引
        selectedRow = Proportion(randomIndex,:);
        % 分别赋值给 Pc 和 Pm
        Pc = selectedRow(1);
        Pm = selectedRow(2);
    else
        randomIndex = unidrnd(action);    %动作指引
        selectedRow = Proportion(randomIndex,:);
        % 分别赋值给 Pc 和 Pm
        Pc = selectedRow(1);
        Pm = selectedRow(2);
    end

else
    for i = 1:20
        if s > 1.5/(20)*(i-1) && s <= 1.5/(20)*(i)
            if Epsilon < rand
                % 找到第一行的最大值及其索引
                [maxValue, maxIndex] = max(Q(i, :));
                % 找到所有最大值的索引
                maxIndices = find(Q(i, :) == maxValue);
                % 如果有多个最大值，随机选择一个索引
                randomIndex = maxIndices(randi(length(maxIndices)));   %动作指引
                selectedRow = Proportion(randomIndex,:);
                % 分别赋值给 Pc 和 Pm
                Pc = selectedRow(1);
                Pm = selectedRow(2);
            else
                randomIndex = unidrnd(action);    %动作指引
                selectedRow = Proportion(randomIndex,:);
                % 分别赋值给 Pc 和 Pm
                Pc = selectedRow(1);
                Pm = selectedRow(2);
            end

        end
    end
end
end