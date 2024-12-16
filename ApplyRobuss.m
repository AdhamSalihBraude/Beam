function Population = ApplyRobuss(Population,RobussType,eta,PenValue)
for i = 1 : length(Population)
    
    if RobussType == 1
        Population(i).F = [Population(i).Feff;Population(i).Var];
    elseif RobussType == 2
        Population(i).F = Population(i).F + (Population(i).Var > eta)*PenValue;
    end

end
