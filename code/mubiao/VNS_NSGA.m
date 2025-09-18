function Foutput=VNS_NSGA(P,DA,CA,HA,Td,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,time,C)
%% 规定参数
tic %开始计时
np = 70;  %种群数
pc = 0.8;  %交叉率
pm = 0.1;  %变异率
[n,A] = size(P);


[Chrom,A1] = Shengcheng(np,n,P,DA,HA);
[Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh); 
a_time=toc;   %记录时间
indices = fit(:, 1) < 80;
values = fit(indices, :);
disp(values);
while a_time < time
    Chrom_son_jc = aberranceJm(Chrom,np,n,pm,P);      %交叉操作
    Chrom_son_by = aberranceJm(Chrom_son_jc,np,n,pm,P);   %变异操作
    %将两代个体合并
    combine_Chrom = zeros(2*np,2*n);
    combine_Chrom(1:np,:) = Chrom(1:np,1:2*n);
    combine_Chrom(np+1:2*np,:) = Chrom_son_by(1:np,1:2*n);
    [Chrom,fit,PL_Fit] = Cal(combine_Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值;
    lastpopulation1=[fit,Chrom];
    frontlabel=quickrank(lastpopulation1,fit);  %得到非支配解，包括重复值----------------------
    lastpopulation2=frontlabel(:,2:end);
    [combine_Chrom] = linyusearch(frontlabel,n,P,DA,HA,Td,C,CA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
    [F,Level] = non_domination_sort_1(np,combine_Chrom(:,2:4));  %非支配排序
    frontlabel=quickrank(lastpopulation2,lastpopulation2(:,1:3));
    if size(frontlabel,1)>2*np
        frontlabel=frontlabel(1:2*np,:);
    else
        while ~(size(frontlabel,1)==2*np)  %费时3秒
            if size(frontlabel,1)<2*np
                [buchongpop,A1]=Shengcheng(2*np-size(frontlabel,1),n,P,DA,HA);
                [fitvalue,~]=Cal(buchongpop,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
                buchongpop=[fitvalue,buchongpop];
                lastpopution=[frontlabel(:,2:end);buchongpop];  %合成种群
                frontvalue=lastpopution(:,1:3);  %%%三个目标值
                frontlabel=quickrank(lastpopution,frontvalue);
            end
        end
    end
    on_inferior=feizhipei(frontlabel(:,2:end),frontlabel(:,2:4));
    min(on_inferior(:,1));
    if (length(find(frontlabel(:,1)==1)))<np  %----------------单目标极值并没有被保存下来
        fnum=0;                                                             %当前前沿面
        while numel(frontlabel(:,1),frontlabel(:,1)<=fnum+1)<=np         %判断前多少个面的个体能完全放入下一代种群
            fnum=fnum+1;
        end   %至少需要前fnum层的个体
        newnum=numel(frontlabel(:,1),frontlabel(:,1)<=fnum);             %前fnum个面的个体数
        Chrom=zeros(np,n*2);
        Chrom(1:newnum,:)=frontlabel(1:newnum,5:end);                   %将前fnum个面的个体复制入下一代 -------
        a=find(frontlabel(:,1)==(fnum+1));
        a_population=frontlabel(a,:);
        a_population1=crowded(a_population(:,2:end),np-newnum); %------对第fnum+1的个体进行排序，选择个体
        Chrom(newnum+1:np,:)=a_population1(:,1:end); %%已经筛选出了popsize个个体 。
    else
        Chrom=crowded(frontlabel(:,2:end),np);
    end

    a_time=toc;   %记录时间
end
[Chrom,fit,PL_Fit1] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
frontvalue = Newranking(fit);
Chrom(:,2*n+1) = frontvalue;
Chrom(:,2*n+2) = fit(:,1);
Chrom(:,2*n+3) = fit(:,2);
Chrom(:,2*n+4) = fit(:,3);
Chrom=select(np,frontvalue,Chrom,fit);
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
% disp('para解集');
% disp(unique_para);
para = unique_para(:,n*2+2:n*2+4);
Foutput = para;
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
    matrix_Chrom = reshape(Chromone, 1, []);
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    Temp = Chromone(indices);   %需要拆卸的步骤
    len = length(Temp);
    Profit =  sum(Pr(Temp));   %拆卸得到的利润
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
        if ca~= 1 && ha ~= 1    
            if j ==1
                p = 1;
                h=1;
                Mm = Mm + 1;
            end
            newP = 1;
            if newP ~= p                 %和前一个工厂属性不同
                g = g + 1; 
                h=1;
                Mm = Mm + 1; 
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
                    Mm = Mm + 1; 
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Crt*t;
                    ec = ec + Ert*t;       %计算能耗
                end
            end
        elseif ca == 1 && ha ~= 1   
            if j == 1
                p = -1;
                h=1;
                Mr = Mr + 1; 
            end
            newP = -1;       
            if newP ~= p                %和前一个工厂属性不同
                g = g + 1;  
                h=1;
                Mr = Mr + 1; 
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
                    Mr = Mr + 1; 
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Cpt*t;
                    ec = ec + Ept*t;       %计算能耗
                end
            end
        else           
            if j == 1
                p = 1;
                h=1;
                Mm = Mm + 1;
            end
            newP = 1;       %工厂属性
            if newP ~= p                   %和前一个工厂不同
                g = g + 1;
                h=1;
                Mm = Mm + 1; 
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t + Ch*t;  
                ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Crt*t - Ch*t;  %注意加上处理危险任务的额外成本
                    ec = ec - Ept*t - Eh*t;       %计算加上处理危险任务的额外能耗
                    g = g + 1;
                    h=1;
                    Mm = Mm + 1; 
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                    ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
                end
            end
        end
        p = newP;          
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
    end
    %计算利润
    profit = Profit - cost - g * Cw - Mm * Cr;
    fit(i,1) = sqrt(F) * 3 / 2;
    fit(i,2) = (-profit) * 2 / 3;
    fit(i,3) = (ec+EC) * 3 / 2;
    PL_Fit(i,1) = sum(fit(i,:))/3;   %加权平均值
end
end



function [F,Level] = non_domination_sort_1(np,fit)
%non_domination_sort 初始种群的非支配排序
%初始化pareto等级为1
np = np*2;
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



%% 领域搜索
function [combine_Chrom] = linyusearch(combine_Chrom,n,P,DA,HA,Td,C,CA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh)
on_inferior_label=find(combine_Chrom(:,1)==1);
on_inferior=combine_Chrom(on_inferior_label,5:end);    %参与邻域搜索的个体集合-----不带目标值
on_inferior_last=combine_Chrom(on_inferior_label,2:end); %包括了函数值在内的非支配个体------带目标值
for i = 1:size(on_inferior,1)
    fitvalue=on_inferior_last(on_inferior_label(i),1:3);  %当前个体的函数值
    Chrom =on_inferior(on_inferior_label(i),:); %第几个个体-------不带目标值
    for j =1:3
        newChrom=linyu(Chrom,n,P,DA,HA,j);  %第几种邻域结构
        [Chrom,fit,PL_Fit] = Cal(newChrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值
        chazhi = fitvalue - fit;
        if length(find(chazhi>=0))==3  %表示解得到改善
            Chrom = newChrom;
            fitvalue = fit;
        end
    end
    combine_Chrom(on_inferior_label(i),2:4) = fitvalue;
    combine_Chrom(on_inferior_label(i),5:end) = Chrom;
end
end

%% 真正的领域搜索
function newChrom=linyu(Chrom,n,P,DA,HA,j)
if j==1      %部分重排
    %Chrom=linyu1(Chrom,n,P);
    Chrom=linyu2(Chrom,n,P,DA,HA);
elseif j==2  %重新选择是否拆卸
    Chrom=linyu2(Chrom,n,P,DA,HA);
else        %部分重排+重新选择是否拆卸
    %Chrom=linyu1(Chrom,n,P);
    Chrom=linyu2(Chrom,n,P,DA,HA);
    Chrom=linyu2(Chrom,n,P,DA,HA);
end
newChrom = Chrom;
end

%% 领域搜索1
function Chrom=linyu1(Chrom,n,P)
% 随机生成2个变异点，从而生成两个片段，中间那个片段的基因即为变异片段，且该片段基因数应大于等于2
% 再将该段基因重新选择
% 举例 个体1 5 2 6 3 8 7 4    产生两个变异点4，7  即1 5 2 6和7 4，变异片段为3 8
% 1 5 2 6和7 4为已选基因，消除这些操作的影响得到3 8中的无约束操作，即可选操作，再依次重新选择操作
% 个人理解：其实就是变异片段里如果有多个无约束操作，即可实现变异效果
A1=[];
Chrom_son=Chrom;
[np,~] = size(Chrom);
% 找到or约束
P0 = P;
%生成or集合
a = 1;
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
for k = 1:n  %遍历行
    if ismember(-1,P0(k,:)) == 1 %存在or约束,记录行数（前驱零件）
        A1(1,a) = k;
        a = a + 1;
    end
end
for k = 1 : n  %遍历列
    if ismember(-1,P0(:,k)) == 1 %存在or约束,记录列数（后驱零件）
        A1(1,a) = k;
        a = a + 1;
    end
end
for i = 1:np
    P0 = P;
    p01=[];p03=[];p2=[];
    p = Chrom(i,:);
    %随机生成2个变异点
    a = randi(n - 3);
    b = randi(n);
    while b - a <= 2   %保证变异段基因大于等于2
        b = randi(n);
    end
    %两个片段之间的基因为变异片段
    p01 = p(1,1:a);
    p03 = p(1,b:n);
    p1 = zeros(1,length(p01) + length(p03));
    p1(1,1:length(p01)) = p01;
    p1(1,length(p01) + 1:length(p01) + length(p03)) = p03;%需要消除影响的零件矩阵
    p2 = p(1,a + 1:b - 1);          %得到变异片段
    p3 = [];
    %消除变异部分外零件的影响
    p00 = [];
    for j = 1:length(p01) + length(p03)
        val = p1(1,j);
        if ismember(val,A1) == 1 && val~=12 && val ~= 22 %若选取的工件在or集合中
            c = find(P0(val,:) == -1); % 找到val or约束的位置
            [~,n2] = size(c);         % or约束对象个数
            %消除or关系影响
            for k = 1:n2
                p00 = P0(:,c(1,k));
                p00(p00 == -1) = 0;
                P0(:,c(1,k)) = p00;
            end
            %消除and关系影响
            p00 = P0(val,:);  %将val对其后继零件影响消除
            p00 = 0;
            P0(val,:) = p00;
        else
            p00 = P0(val,:);    %消除and影响
            p00 = 0;
            P0(val,:) = p00;
        end
    end
    %变异部分重新排序
    for j = 1:length(p2)
        A2 = [];
        %生成优先任务集合
        c  = 1;
        for k = 1:n
            if any(P0(:,k)) == 0 %k零件没有任何前驱约束
                A2(1,c) = k;
                c = c + 1;
            end
        end
        [~,n0] = size(A2);
        %A2只保留p2中的操作
        A20 = intersect(A2, p2);
        [~,n1] = size(A20);
        f = 1;%保证选出的新任务没有安排过
        while f == 1
            val = A20(1,randi(n1));
            if ismember(val,p3) == 0   %在p3中不存在，也就是没有被安排过
                f = 0;
            end
        end
        p00 = [];
        if ismember(val,A1) == 1%若选取的工件在or集合中
            %消除or关系影响
            c = find(P0(val,:) == -1);
            [~,n0] = size(c);%约束对象个数
            for k = 1:n0
                p00 = P0(:,c(1,k));
                p00(find(p1 == -1)) = 0;
                P0(:,c(1,k)) = p00;
            end
            p00 = P0(val,:);    %消除and影响
            p00 = 0;
            P0(val,:) = p00;
        else
            p00 = P0(val,:);    %消除and影响
            p00 = 0;
            P0(val,:) = p00;
        end
        p3(1,j) = val;
    end
    Chrom_son(i,a+1:b-1) = p3(1,:);
end
%根据前面存储的Temp，重新置0，1
for i = 1:np
    ChromOne = Chrom_son(i,:);
    % 获取第一行不等于0的元素的逻辑索引
    nonZeroIndices = Temp(i, :) ~= 0;
    % 使用逻辑索引获取需要拆卸的步骤
    dis = Temp(i, nonZeroIndices);
    [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
    Chrom_son(i,n+idx1) = 1;
end
Chrom = Chrom_son;
end

%% 领域搜索2
function Chrom=linyu2(Chrom,n,P,DA,HA)
%重新选择是否拆卸
[np,~] = size(Chrom);
P0 = P;
indices_ca = find(DA == 1);  % 需求操作
indices_ha = find(HA == 1);  % 危险操作
indices = [indices_ca,indices_ha];    %所有有需求和危险的操作
len = length(indices);
A = indices;
dis = [];    %拆卸集合
A2 = [0]; %保证进入循环
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
    if len > 0
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
        ChromOne = Chrom(j,1:n);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
        [~, idx2] = ismember(random_dis_fm, ChromOne);  %随机的拆卸步骤
        Chrom(j,n+idx2) = 1;
    else
        ChromOne = Chrom(j,1:n);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
    end
end
end

%%非支配解分层排序
function [frontlabel]=quickrank(lastpopution,frontvalue)
lastpopution=[zeros((size(lastpopution,1)),1),lastpopution];
frontlabel=[];
ceng=0;
while ~isempty(lastpopution)
    %      size(lastpopution,1)
    %      pause(2)
    label=[];
    ceng=ceng+1;
    if size(lastpopution,1)>1
        for i1=1:size(lastpopution,1)
            n=0;
            pnum=[];
            qnum=[];
            for j1=1:size(lastpopution,1)
                if i1==j1
                    continue;%结束当前循环,不往下执行
                end
                chazhi=frontvalue(i1,:)-frontvalue(j1,:);
                if (length(find(chazhi(1,:)>=0))==3)&(length(find(chazhi(1,:)>0))>=1)    %%被别人支配
                    n=1;
                    break;
                end
                if  ((length(find(chazhi(1,:)<=0))==3)&(length(find(chazhi(1,:)<0))>=1))||(((length(find(chazhi(1,:)>0))>=1))&(length(find(chazhi(1,:)<0))>=1))
                    plabel=j1;%%%支配解
                    pnum=[pnum; plabel]; %%%支配解标签
                end
                if (length(find(chazhi(1,:)==0))==3)  %%相同点
                    qlabel=j1;  %%%相同点标签
                    qnum=[qnum;qlabel];  %%%相同点-------------------标记相同点的
                end
            end
            pnum1=length(pnum);  %66
            qnum1=length(qnum);
            zonggeshu=pnum1+qnum1;
            if  (n==0)&(qnum1==0)
                lastpopution(i1,1)=ceng;
                label1=i1;
                label=[label;label1];
            end
            if  (n==0)&(zonggeshu==(size(lastpopution,1)-1))&(qnum1>0)
                %qnum=[qnum;i1];
                lastpopution(i1,1)=ceng;
                label1=i1;
                label=[label;label1];
            end
        end
        zhongjian1=lastpopution(label,:);
        frontlabel=[frontlabel;zhongjian1];
        lastpopution(label,:)=[];
        frontvalue(label,:)=[];
    else
        lastpopution(1,1)=ceng;
        zhongjian1=lastpopution(1,:);
        frontlabel=[frontlabel;zhongjian1];
        lastpopution(1,:)=[];
        frontvalue(label,:)=[];
    end
end
end

%% 只是为了得到非支配解
function on_inferior=feizhipei(lastpopution,fitvalue)
%一个解不受任何解支配就是非支配解
%目前是：最小时间，最大利润，最小能耗
%值相减，有正有负，则互不支配，存在都是正的就是存在被支配的解，均为0则去除重复点。
%目前对目标函数2进行了取倒数，处理，就都是求最小了。
simi_m=[];  %存储一系列相同值
on_inferior=[];  %存储非支配解
for i=1:size(lastpopution,1)
    if length(find(simi_m==i))==0
        label=1:size(lastpopution,1);
        fit_value=fitvalue(i,:);  %函数值
        label(i)=[];
        fit_value_1=repmat(fit_value,length(label),1);
        chazhi=fit_value_1-fitvalue(label,:); %差值
        on_inferior_label=1;
        simi_position=[]; %存储相同值的位置
        for j=1:length(label)
            if chazhi(1)<0.000000001&&chazhi(1)>-0.000000001
                chazhi(1)=0;
            end
            if chazhi(2)<0.000000001&&chazhi(2)>-0.000000001
                chazhi(2)=0;
            end
            if chazhi(3)<0.000000001&&chazhi(3)>-0.000000001
                chazhi(3)=0;
            end
            if length(find(chazhi(j,:)==0))==3
                simi_position=[simi_position,label(j)];
            end
            if (length(find(chazhi(j,:)>=0))==3)&&(length(find(chazhi(j,:)>0))>1)
                on_inferior_label=0;
                break
            end
        end
        if  on_inferior_label==1
            if length(simi_position)>0
                simi_m=[simi_m,simi_position];
            end
            on_inferior=[on_inferior,i];
        end
    end
end
on_inferior=lastpopution(on_inferior,:);
end

%% 拥挤度的计算
function population=crowded(on_inferior,popunumber)
functionvalue=on_inferior(:,1:3);
distancevalue=zeros(size(functionvalue,1),1);
fmax=max(on_inferior(:,1:3),[],1);                       %三个目标上的最大值
fmin=min(on_inferior(:,1:3),[],1);                          %三个目标上的最小值
%对数据进行归一化
for i=1:size(functionvalue,1)
    functionvalue(i,1)=(functionvalue(i,1)-fmin(1))/(fmax(1)-fmin(1));
    functionvalue(i,2)=(functionvalue(i,2)-fmin(2))/(fmax(2)-fmin(2));
    functionvalue(i,3)=(functionvalue(i,3)-fmin(3))/(fmax(3)-fmin(3));
end
for l=1:3     %分目标计算每个目标上popu各个体的拥挤距离
    [~,newsite]=sortrows(functionvalue(:,l));
    distancevalue(newsite(1))=inf;
    distancevalue(newsite(end))=inf;
    for m=2:(size(functionvalue,1)-1)
        distancevalue(newsite(m))=distancevalue(newsite(m))+(functionvalue(newsite(m+1),l)-functionvalue(newsite(m-1),l))/(fmax(l)-fmin(l));
    end
end
%         distancevalue
distance_value=zeros(size(functionvalue,1),1);
for l=1:3
    [~,newsite]=sortrows(functionvalue(:,l));
    for m=2:(size(functionvalue,1)-1)
        distance_value(newsite(m))= distance_value(newsite(m))+(abs((functionvalue(newsite(m+1),l)-functionvalue(newsite(m-1),l)))-distancevalue(newsite(m)))^2;
    end
end
distance_value=distance_value/3;
[~,newsite1]=sortrows(distance_value);  %升序
newsite1=fliplr(newsite1');  %反转
newsite1=newsite1(1:popunumber);
population=on_inferior(newsite1,4:end);
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

%% Selch=select1('rws',ChromNumber,Fitnv);  %根据Fitnv适应度值轮盘赌进行选择
function Chrom=select(NIND,frontvalue,FChrom,Objv)
%---计算拥挤距离/选出下一代个体==============================================zxs
fnum=0;                                                                %当前前沿面
while numel(frontvalue,frontvalue<=fnum+1)<NIND              %判断前多少个面的个体能完全放入下一代种群
    fnum=fnum+1;
end

newnum=numel(frontvalue,frontvalue<=fnum);                     %前fnum个面的个体数
Chrom(1:newnum,:)=FChrom(frontvalue<=fnum,:);                  %将前fnum个面的个体复制入下一代
popu=find(frontvalue==fnum+1);                                 %popu记录第fnum+1个面上的个体编号
distancevalue=zeros(size(popu));                               %popu各个体的拥挤距离
fmax=max(Objv(popu,:),[],1);                                   %popu每维上的最大值
fmin=min(Objv(popu,:),[],1);                                   %popu每维上的最小值
for i=1:size(Objv,2)                                           %分目标计算每个目标上popu各个体的拥挤距离
    [~,newsite]=sortrows(Objv(popu,i));                        %popu里对第一维排序之后的位置
    distancevalue(newsite(1))=inf;
    distancevalue(newsite(end))=inf;
    for j=2:length(popu)-1
        distancevalue(newsite(j))=distancevalue(newsite(j))+(Objv(popu(newsite(j+1)),i)-Objv(popu(newsite(j-1)),i))/(fmax(i)-fmin(i));
    end
end
popu=-sortrows(-[distancevalue;popu]')';                             %按拥挤距离降序排序第fnum+1个面上的个体
Chrom(newnum+1:NIND,:)=FChrom(popu(2,1:NIND-newnum),:);          %将第fnum+1个面上拥挤距离较大的前popnum-newnum个个体复制入下一代
end