function [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,N)
    %% 删除一些差解
    [~,uni] = unique(Population.objs,'rows');
    Population = Population(uni);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);
    N          = min(N,length(Population));

    %% 非支配排序
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% 对每个解计算拥挤度距离
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    %% 基于拥挤度距离选择最后前沿面上的解
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% 获取进入下一次迭代的种群
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
end