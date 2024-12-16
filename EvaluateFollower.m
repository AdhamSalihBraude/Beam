function PopulationFollower = EvaluateFollower(PopulationFollower,NonDominatedLeader,FeffType,VarType,ProblemName,Limits)
DH = [NonDominatedLeader(:).X];

for i = 1 : length(PopulationFollower)
    Feff = [];
    x = PopulationFollower(i).X;
    
    Xd = x + DH;
    
    for idxDec = 1 : size(Limits,2)
        Xi = Xd(idxDec,:);
        uLimit = Limits(2,idxDec);
        lLimit = Limits(1,idxDec);
        Xi(Xi>uLimit) = uLimit;
        Xi(Xi<lLimit) = lLimit;
        Xd(idxDec,:) = Xi;
    end
     F = EvalFuntion(Xd,ProblemName);   
    if FeffType == 2
        Feff = max(F,[],2); %max
    elseif FeffType == 3
        Feff = min(F,[],2); %max
    else
        Feff = mean(F,2); %mean
    end

    if VarType == 1
        if size(F,2) > 1
        Var = norm(F-Feff)/norm(F);
        else
            Var = 0;
        end
    else
        if size(F,2) > 1
        Fmax = max(F,[],2);
        Fmin = min(F,[],2);
        Fnorm = (F-Fmin)./(Fmax-Fmin);
        
        Var = sqrt((max(Fnorm(:,1))-min(Fnorm(:,1)))^2 + (max(Fnorm(:,2))-min(Fnorm(:,2)))^2);
        else
            Var = 0;
        end
    end
    if isnan(Var)
        Fmax = max(F,[],2);
        Fmin = min(F,[],2);
        Var = norm(Fmax-Fmin);
    end
    % Feff = [Feff;Var];
    PopulationFollower(i).Feff = -Feff;
    PopulationFollower(i).Var = -Var;
    PopulationFollower(i).F = -Feff;
    PopulationFollower(i).Var;
end