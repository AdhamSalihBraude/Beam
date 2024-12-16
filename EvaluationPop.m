function Population = EvaluationPop(Population,ProblemName)
    for i = 1 : length(Population)
        Population(i).F = EvalFuntion(Population(i).X,ProblemName);
    end
end