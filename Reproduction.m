function Offsprings = Reproduction(MatingPool,CrossOverRate,MutationRate,Limits)
N =  length(MatingPool);
Offsprings = MatingPool;
%% CrossOver
for i = 1 : 2 : N

    parent1 = MatingPool(i).X;
    parent2 = MatingPool(i+1).X;
    [child1, child2] = SBX(parent1,parent2,CrossOverRate,Limits);
    Offsprings(i).X = child1;
    Offsprings(i+1).X = child2;

    % init Fit
    Offsprings(i).Rank = inf;
    Offsprings(i+1).Rank = inf;
    Offsprings(i).CrowdDis = 0;
    Offsprings(i+1).CrowdDis = 0;
    Offsprings(i).F = [];
    Offsprings(i+1).F = [];
end

%% mutation
for i = 1 : N
    % Eprameters
    parent = Offsprings(i).X;
    child = PolyMutation(parent,MutationRate,Limits);
    Offsprings(i).X = child;
end
