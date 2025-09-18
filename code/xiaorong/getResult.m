clc;
clear all;
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
gen=20;
shu = 3;  %测试算法个数
result=cell(shu,gen);     %存放解集
Objv=[];                %用于得到真实Pareto解集
jishu=zeros(shu,gen);     %用于记录每一解集中在真实Pareto前沿解集中的个数
igd=zeros(shu,gen);   %存放IGD值
HV=zeros(shu,gen);   %存放HV值
qkpl = cell(1,gen);   %存放结果
kpl = cell(1,gen);   %存放结果
dsKOA = cell(1,gen);  %存放结果
load('Q.mat');
for i=1:gen
       [X1,time] = QDS_KOA(P,DA,CA,HA,Td,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,C,Q);
       X2=KOA(P,DA,CA,HA,Td,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,time,C);
       X3=DS_KOA(P,DA,CA,HA,Td,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,time,C);
       result{1,i}=X1;
       result{2,i}=X2;
       result{3,i}=X3;
       Objv=cat(1,Objv,X1);
       Objv=cat(1,Objv,X2);
       Objv=cat(1,Objv,X3);
       qkpl{1,i} = X1; 
       kpl{1,i} = X2;
       dsKOA{1,i} = X3;
end
save('qkpl.mat', 'qkpl');
save('kpl.mat', 'kpl');
save('dsKOA.mat', 'dsKOA');

frontvalue = Newranking(Objv);
Foutput=sortrows(Objv(frontvalue<=1,:));  %最终结果:种群中非支配解的函数值
[a,b]=size(Foutput);
[zu,ci]=size(result);
for i=1:zu
    for j=1:ci
        x=result{i,j};
        [Gen,~]=size(x);
        for k=1:Gen
            for m=1:a
                if x(k,1)==Foutput(m,1) && x(k,2)==Foutput(m,2)   %判断各解集中是否有存在于真实pareto解之中，有则对应位置加1
                    jishu(i,j)=jishu(i,j)+1;
                end
            end
        end
    end
end
pareto_frontier=Foutput;  %真实Pareto前沿

for i=1:zu
    for j=1:ci
        pareto_set=result{i,j};
        A=calculate_igd(pareto_frontier, pareto_set);
        igd(i,j)=A;
    end
end
referencePoint=pareto_frontier(2,:);
for i=1:zu
    for j=1:ci
        pareto_set=result{i,j};
        A=calculateMultiObjectiveHV(pareto_set, referencePoint);
        HV(i,j)=A;
    end
end

sum_per_row = sum(igd, 2);
for i = 1:zu
    sum_per_row(i,1) = sum_per_row(i,1)/gen;
end

sum_per_HV = sum(HV, 2);
for i = 1:zu
    sum_per_HV(i,1) = sum_per_HV(i,1)/gen;
end

disp(jishu);
disp('IGD平均值：');
disp(sum_per_row);
disp('HV平均值：');
disp(sum_per_HV);
figure

%% 
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



function hv = calculateMultiObjectiveHV(approximationSet, referencePoint)
    % approximationSet: 近似解集，每一行是一个解，每一列是一个目标的值
    % referencePoint: 参考点，一个包含目标函数最大值的向量
    
    % 计算解集的大小
    [numSolutions, numObjectives] = size(approximationSet);
    
    % 初始化HV值为0
    hv = 0;
    % 对每个目标进行计算
    for i = 1:numObjectives
        % 取出目标函数的值
        objectiveValues = approximationSet(:, i);
        
        % 计算该目标下的最大值和最小值
        maxObjective = max(objectiveValues);
        minObjective = min(objectiveValues);
        
        % 如果最大值和最小值不相等，计算归一化的目标函数值
        if maxObjective ~= minObjective
            normalizedValues = (maxObjective - objectiveValues) / (maxObjective - minObjective);
            % 计算超立方体的体积，并累加到HV值中
            hv = hv + prod(referencePoint(i) - normalizedValues);
        end
    end
    
    % 返回HV值
    hv = hv / numSolutions;
end
