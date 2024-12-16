function Population = InitPop(PopSize,Limits)
N = size(Limits,2);
Population(PopSize).X = [];
Population(PopSize).F = [];
Population(PopSize).Var = [];
Population(PopSize).Feff = [];
Population(PopSize).Rank = [];
Population(PopSize).CrowdDis = [];
for i = 1 : PopSize
    Population(i).X = [Limits(1,:) + rand(1,N).*(Limits(2,:) - Limits(1,:))]';
end
