function Foutput = KOA(P,DA,CA,HA,Td,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,time,C)
% 四目标：站的数量，平滑程度，能耗，成本
% 规定参数
tic %开始计时
np = 50;  %行星数
Tmax = 2; %最大迭代次数
dim = 3; %目标数
[n,A] = size(P); % n代表零件数
favg = 0;   %存储平均适应度值
Q = [];
pareto_frontier = [10 -3000 40];  %计算IGD值真实点
ub = 1;
lb = 0;
Sun_Pos = zeros(1,n*2); % 一个矢量，包含迄今为止最好的解决方案，代表太阳
Sun_Score = inf(1,dim); % 一个标量变量，包含迄今为止最好的分数
Convergence_curve = zeros(1,Tmax); %存储最优解
%---------------------控制参数------------------------%
%
Tc=3;
M0=0.1;
lambda=15;
% 轨道偏心率(e)
orbital=rand(1,np); %% Eq.(4)   生成与种群数个数相同的偏心率（0-1）
%%轨道周期(T)
T=abs(randn(1,np)); %% Eq.(5)   根据正态分布生成与种群数个数相同的偏心率
% 种群初始化

[Chrom,A1] = Shengcheng(np,n,P,DA,HA);
positions=rand(1,n); % 初始化行星的位置
% 对Positions排序
[A,B] = sort(positions);
for j = 1:np
    Chromone = Chrom(j,1:n);
    for k = 1:n
        Positions(j,Chromone(k)) = A(1,k);
    end
end

%把第一代个体的IGD值求和
t = 0;   %函数求值计数器
[Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值
for j = 1:np
    % %计算目标值
    if calculate_igd(pareto_frontier,fit(j,:)) <= calculate_igd(pareto_frontier,Sun_Score)
        Sun_Score=fit(j,:); % 更新到目前为止最好的分数
        Sun_Pos=Positions(j,:); % 更新目前为止最好的解决方案
        Sun_Pos(1,n+1) = fit(j,1);
        Sun_Pos(1,n+2) = fit(j,2);
        Sun_Pos(1,n+3) = fit(j,3);
    end
end
%进化
a_time=toc;   %记录时间
while a_time < time
    [Order] = sort(PL_Fit);
    % 函数求值t时的最差适应度值
    worstFitness = Order(np); %% Eq.(11)
    M = M0*(exp(-lambda*(t/Tmax))); %% Eq. (12) %M是一个随时间（t）呈指数下降的函数，用于控制搜索精度
    % 计算机R表示目前最佳解和第i个解之间的欧氏距离
    for j = 1:np
        R(j) = 0;
        for k = 1:dim
            R(j)=R(j)+(Sun_Pos(n+k)-fit(j,k))^2; %% Eq.(7)    %欧式距离应用到拆卸线得想一下修改
        end
        R(j)=sqrt(R(j));
    end
    for j = 1:np
        Sum = 0;
        for k = 1:np
            Sum=Sum+(PL_Fit(k)-worstFitness);
        end
        MS(j) = rand*(sum(Sun_Score)/4-worstFitness)/(Sum); % Eq.(8)
        m(j)=(PL_Fit(j)-worstFitness)/(Sum); %% Eq.(9)
    end
    %Step 2: 定义引力(F)
    %根据万有引力定律计算太阳和第i颗行星之间的引力
    for j = 1:np
        Rnorm(j)=(R(j)-min(R))/(max(R)-min(R)); %% The normalized R (Eq.(24))XS和Xi之间的欧几里德距离归一化
        MSnorm(j)=(MS(j)-min(MS))/(max(MS)-min(MS)); %% The normalized MS MS质量归一化
        Mnorm(j)=(m(j)-min(m))/(max(m)-min(m)); %% The normalized m      m质量归一化
        Fg(j)=orbital(j)*M*((MSnorm(j)*Mnorm(j))/(Rnorm(j)*Rnorm(j)+eps))+(rand); %% Eq.(6)  计算引力
    end
    % A1表示物体I在时刻t的椭圆轨道的半长轴，
    for j = 1:np
        a1(j)=rand*(T(j)^2*(M*(MS(j)+m(j))/(4*pi*pi)))^(1/3); %% Eq.(23)
    end
    for j =1:np
        % a2是一个从-1到-2逐渐减小的循环控制参数
        a2=-1+-1*(rem(t,Tmax/Tc)/(Tmax/Tc)); %% Eq.(29)
        % 而n1是从1到−2的线性递减因子
        n1=(a2-1)*rand+1; %% Eq.(28)
        a=randi(np);  % 取出一个随机解
        b=randi(np); % 取出一个随机解
        rd=rand(1,n); %根据正态分布生成的向量 r5
        r=rand; %r是0到1的一个随机数 r4
        U1=rd<r; %% Eq.(21)    %好像反了需要最后效果求证
        O_P=Chrom(j,:); %% 存储第i个解
        O_Positions = Positions(j,:);%% 存储原本位置
        % Step 6: 更新行星
        if rand<rand
            % h是用于控制在时间t时太阳与当前行星之间的距离的自适应因子
            h=(1/(exp(n.*randn))); %% Eq.(27)
            %基于三个解决方案的平均向量:当前解决方案、迄今最佳解决方案和随机选择的解决方案
            Xm=(Positions(a,:)+Sun_Pos(1,1:n)+Positions(j,:))/3.0;
            Positions(j,:)=Positions(j,:).*U1+(Xm+h.*(Xm-Positions(b,:))).*(1-U1); %% Eq.(26)
        else
            % Step 3: 计算物体速度
            % 与当前行星的搜索方向相反或偏离的标志
            if rand<0.5 %% Eq.(18)
                f=1;
            else
                f=-1;
            end
            L=(M*(MS(j)+m(j))*abs((2/(R(j)+eps))-(1/(a1(j)+eps))))^(0.5); %% Eq.(15)  eps为一个极小值数
            U=rd>rand(1,n); %% 二值向量Eq.(17)
            if Rnorm(j)<0.5 %% Eq.(13)
                M=(rand.*(1-r)+r); %% Eq.(16)
                l=L*M*U; %% Eq.(14)
                Mv=(rand*(1-rd)+rd); %% Eq.(20)
                l1=L.*Mv.*(1-U);%% Eq.(19)
                V(j,:)=l.*(2*rand*Positions(j,:)-Positions(a,:))+l1.*(Positions(b,:)-Positions(a,:))+(1-Rnorm(j))*f*U1.*rand(1,n).*(ub-lb); %% Eq.(13a)
            else
                U2=rand>rand; %% Eq. (22)
                V(j,:)=rand.*L.*(Positions(a,:)-Positions(j,:))+(1-Rnorm(j))*f*U2*rand(1,n).*(rand*ub-lb);  %% Eq.(13b)
            end %% End IF

            % Step 4: 逃离局部最优
            % 更新标志f到相反的方向或离开当前行星的搜索方向
            if rand<0.5 %% Eq.(18)
                f=1;
            else
                f=-1;
            end
            % Step 5 跟新位置
            Positions(j,:)=((Positions(j,:)+V(j,:).*f)+(Fg(j)+abs(randn))*U.*(Sun_Pos(1,1:n)-Positions(j,:))); %% Eq.(25)
        end
        %根据更新后的位置生成新种群,并检查是否合理,并且重新生成是否拆卸反序列
        Positionsone = Positions(j,:);
        [A,Chromone] = sort(Positionsone);
        % 对更新个体合理性进行判断
        for k = 1:n
            a = Chromone(1,k);
            A1 = find(P(:, a) == -1);  %前置或集合
            A2 = find(P(:, a) == 1);   %前置与集合
            if ~isempty(A1)
                if ~any(ismember(Chromone(1,1:k-1),A1))   % 或集合不满足
                    Chromone(1,:) = O_P(1,1:n);              % 以优秀父代个体作为子代
                    Positions(j,:) = O_Positions;
                    break;
                end
            end
            if ~isempty(A2)
                if ~all(ismember(A2, Chromone(1,1:k-1)))  % 与集合不满足
                    Chromone(1,:) = O_P(1,1:n);              % 以优秀父代个体作为子代
                    Positions(j,:) = O_Positions;
                    break;
                end
            end
        end
        %存放个体
        newChrom(j,1:n) = Chromone(1,:);
        %确定是否需要重新选择是否拆卸
        if ~isequal(Chromone, O_P(1:numel(Chromone)))
            newChrom(j,n+1:n*2) = 0;
            %重新选择
            P0 = P;
            indices_ca = find(DA == 1);  % 需求操作
            indices_ha = find(HA == 1);  % 危险操作
            indices = [indices_ca,indices_ha];    %所有有需求和危险的操作
            len = length(indices);
            A = indices;
            dis = [];    %拆卸集合
            while ~isempty(A2)
                len = length(indices);
                for h = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                    temp = indices(h);   %操作
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
                    for h = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                        temp = random_dis(h);   %操作
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
                [~, idx1] = ismember(dis, Chromone);     %必须拆卸步骤
                newChrom(j,n+idx1) = 1;
                [~, idx2] = ismember(random_dis_fm, Chromone);  %随机的拆卸步骤
                newChrom(j,n+idx2) = 1;
            else
                [~, idx1] = ismember(dis, Chromone);     %必须拆卸步骤
                newChrom(j,n+idx1) = 1;
            end
        else
            newChrom(j,n+1:n*2) = O_P(1,n+1:n*2);
        end
        %计算函数值
        [Chrom_nouse,fit_one,PL_Fit1] = Cal(newChrom(j,:),n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值;
        if calculate_igd(pareto_frontier,fit_one) <= calculate_igd(pareto_frontier,fit(j,:)) % 保留优解
            fit(j,:)=fit_one; % 保留优解
            % 更新全局最优解
            if calculate_igd(pareto_frontier,fit(j,:))<calculate_igd(pareto_frontier,Sun_Score) % 目前解比全局最优解好
                Sun_Score=fit(j,:); % 更新全局最优解适应度值
                Sun_Pos=Positions(j,:); % 更新全局最优解方案
                Sun_Pos(1,n+1) = fit(j,1);
                Sun_Pos(1,n+2) = fit(j,2);
                Sun_Pos(1,n+3) = fit(j,3);
            end
        else
            newChrom(j,:)=O_P;
        end %% End IF
    end
    t=t+1; %% 迭代次数加1
    if t>Tmax %% 检查终止条件
        break;
    end %% End IF
    Chrom = [newChrom;Chrom];
    [Chrom,fit,PL_Fit1] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
    Chrom = [fit,Chrom];
    frontlabel=quickrank(Chrom,fit);
    %去除重复值
    %frontlabel=xipai(frontlabel,5);
    if size(frontlabel,1)>2*np
        frontlabel=frontlabel(1:2*np,:);
    else
        while ~(size(frontlabel,1)==2*np)  %费时3秒
            if size(frontlabel,1)<2*np
                [buchongpop,A1]=Shengcheng(2*np-size(frontlabel,1),n,P,DA,HA);
                [~,fitvalue,~]=Cal(buchongpop,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
                buchongpop=[fitvalue,buchongpop];
                lastpopution=[frontlabel(:,2:end);buchongpop];  %合成种群
                frontvalue=lastpopution(:,1:3);  %%%三个目标值
                frontlabel=quickrank(lastpopution,frontvalue);
                %frontlabel=xipai(frontlabel,3);
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
    a_time=toc;
end
[Chrom,fit,PL_Fit1] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
frontvalue = Newranking(fit);
Chrom(:,2*n+1) = frontvalue;
Chrom(:,2*n+2) = fit(:,1);
Chrom(:,2*n+3) = fit(:,2);
Chrom(:,2*n+4) = fit(:,3);
Chrom=select(np,frontvalue,Chrom,fit);
% [F,Level] = non_domination_sort(np,fit);
% Chrom(:,2*n+1) = Level;
% Chrom(:,2*n+2) = fit(:,1);
% Chrom(:,2*n+3) = fit(:,2);
% Chrom(:,2*n+4) = fit(:,3);
% Chrom(:,2*n+5) = fit(:,4);
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
para = unique_para(:,n*2+2:n*2+4);
Foutput = para;

% para = unique_para(:,n*2+2:n*2+5);
% Foutput = para;

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
        if ca~= 1 && ha ~= 1     
            if j ==1
                p = 1;
                h=1;
                Mm = Mm + 1;
            end
            newP = 1;
            if newP ~= p                 
                g = g + 1;  
                h=1;
                Mm = Mm + 1;  
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t;  
                ec = ec + Ert*t;      
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t;
                ec = ec + Ert*t;       
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
            if newP ~= p                
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
                    ec = ec - Ept*t;      
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
            newP = 1;      
            if newP ~= p                   %和前一个工厂不同
                g = g + 1; 
                h=1;
                Mm = Mm + 1;  
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
    fit(i,1) = sqrt(F)* 3 / 2;
    fit(i,2) = (-profit)* 2 / 3;
    fit(i,3) = (ec+EC)* 3 / 2;
    PL_Fit(i,1) = sum(fit(i,:))/3;   %加权平均值
end
end

%% 计算IGD指标
function igd_value = calculate_igd(pareto_frontier, pareto_set)

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