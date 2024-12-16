function PopulationUpdated = EliteFullSorting(Population,Offsprings,N,elitePer)

FitnessValues = [Population(:).F]';
[~,ia,~] = unique(FitnessValues,'rows');
Population = Population(ia);
AllPop = [Population,Offsprings];
AllPop = CalcRankAndDistance(AllPop);
AllPop(isnan([AllPop(:).CrowdDis]))=[];
RankAndDistMatrix = [[AllPop(:).Rank]',[AllPop(:).CrowdDis]'];
[~, I ] = sortrows(RankAndDistMatrix,[1 -2]);

if elitePer == 1

    if length(I)>N
        PopulationUpdated = AllPop(I(1:N));
    else
        PopulationUpdated = AllPop;
        while length(PopulationUpdated) < N
            PopulationUpdated = [PopulationUpdated,AllPop(randi(length(AllPop)))];
        end
    end
else


    PopulationUpdated = Offsprings(1:round((1-elitePer)*N));

    while length(PopulationUpdated) < N
        PopulationUpdated = [PopulationUpdated,AllPop(randi(length(AllPop)))];
    end
end

