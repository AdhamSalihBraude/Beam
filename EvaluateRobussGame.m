function PopulationLeader = EvaluateRobussGame(PopulationLeader,PopulationFollower,FeffType,VarType,ProblemName,Limits)
DH = [PopulationFollower(:).X];
for i = 1 : length(PopulationLeader)
    Feff = [];
    x = PopulationLeader(i).X;

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
    if FeffType == 1
        Feff = mean(F,2); %mean
    elseif FeffType == 2
        Feff = max(F,[],2); %max
    elseif FeffType == 3
        Feff = min(F,[],2); %min
    else
        Feff = PopulationLeader(i).F; %No uncer
    end
    if size(F,2) > 1


        if VarType == 1
            Var = norm(F-Feff)/norm(F);
        elseif VarType == 2
            Rk = 0;
            for k = 1 : size(F,2)
                Rk = Rk + norm(F(:,k)-PopulationLeader(i).F)/norm(Xd(:,k)-PopulationLeader(i).X);
            end

            Var = Rk/k;

        elseif VarType == 3
            Fmax = max(F,[],2);
            Fmin = min(F,[],2);
            Fnorm = (F-Fmin)./(Fmax-Fmin);
            Var = sqrt((max(Fnorm(:,1))-min(Fnorm(:,1)))^2 + (max(Fnorm(:,2))-min(Fnorm(:,2)))^2);
        else
            Var = [];
        end
    else
        Var = 0;
    end
    if isnan(Var)
        Fmax = max(F,[],2);
        Fmin = min(F,[],2);
        Var = norm(Fmax-Fmin);
    end
    %Feff = [Feff;Var];
    PopulationLeader(i).Feff = Feff;
    PopulationLeader(i).Var = Var;

end