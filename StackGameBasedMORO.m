%This code is related to:
% ------------------------------- Reference --------------------------------
% 
% Salih, A., & Matalon, E. E. (2024). Solving multi-objective robust 
% optimization problems via Stakelberg-based game model. 
% Swarm and Evolutionary Computation, 91, 101734.
%------------------------------- Copyright --------------------------------
% You are free to use this code forresearch purposes.
% Please note that some code was taken from PlatEMO. In addition to 
% acknowledge the use of this code, All publications which
% use this code should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

clc
clear
clf
close all
rng('shuffle')
GenRatioFollower = 10;
PopSizeLeader = 200; %number of individuals

PopSizeFollower = 100;
maxFE = 50000;
MaxGen = ceil(maxFE/(PopSizeLeader)); % maximal number of generations
NpopTosave = 25;
saveEachGen = 10;
%% Problem Def
for FeffType = [0,2]
    ProblemString = ['TP1'; 'TP2'; 'TP3'; 'TP4'; 'TP5'; 'TP6'; 'TP7'; 'TP8'; 'TP9'; 'Tru'; 'bea'];
    for ProbIdx = 1 : size(ProblemString,1)
        ProblemName = ProblemString(ProbIdx,:); % 'TP1', 'TP2', 'TP3', 'TP4', 'TP5', 'TP6', 'TP7', 'TP8', 'TP9', 'Truss', 'beam'
        D = 5;
        delta = 0.05;
        CrossOverRate = 1;
        MutationRate = 0.5;
        PlotFlag = 0;
        PlotEnd = 0;

        NumRuns = 31;
        RobustType = 1;
        % 1- mean, 2- max, 3- min, otherwise- no uncer
        VarType = 3; % 1- Deb and Gupta 2005, 2- Gaspar-Cunha 2014, 3- Diag, 4- no uncer
        FeffTypeFollower = 1; % 1- mean, 2- max, 3- min, otherwise- no uncer
        VarTypeFollower = 3; % 1- Deb and Gupta 2005, 2- Gaspar-Cunha 2014, 3- Diag, 4- no uncer
        eta = 1000;
        PenValue = 1e2;
        elitePer = 1;


        if strcmp(ProblemName,'TP1')
            LimitsLeader = [ zeros(1,D);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP2')

            LimitsLeader = [ zeros(1,D);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP3')

            LimitsLeader = [ zeros(1,D);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP4')

            LimitsLeader = [ zeros(1,D);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP5')

            LimitsLeader = [ zeros(1,D);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP6')

            LimitsLeader = [ 0, -ones(1,D-1);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP7')

            LimitsLeader = [ 0, -ones(1,D-1);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP8')

            LimitsLeader = [ 0, 0,-ones(1,D-2);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP9')

            LimitsLeader = [ 0, 0,-ones(1,D-2);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'TP10')

            LimitsLeader = [ zeros(1,D);...
                ones(1,D)];

            LimitsFollower = [ones(1,D)*(-delta);...
                ones(1,D)*(delta) ];

        elseif strcmp(ProblemName,'Tru')
            ProblemName = 'Truss';
            LimitsLeader = [ 20 20 1000;...
                50 50 3000];

            LimitsFollower = [-1 -1 -50 ;...
                1 1 50 ];

        elseif strcmp(ProblemName,'bea')
            ProblemName = 'beam';
            LimitsLeader = [ 100 100 9 9;...
                800 500 50 50];

            LimitsFollower = [-30 -30 -3 -3 ;...
                30 30 3 3 ];
        else
            error([ProblemName, '- Error in the problem name (unrecognized problem)'])
        end
        FolderName = [ProblemName,'_Ftype_',num2str(FeffType),'_', char(datetime('today','Format','dd_MM_yyy'))];
        mkdir( FolderName)
        if size(LimitsFollower,2) ~= size(LimitsLeader,2)
            error("Error: Delta Length is " + length(Delta) + " For " + size(Limits,2) + " Vars")
        end



        % Saving Run Parameters
        ParamData.ProblemName = ProblemName;
        ParamData.GenRatioFollower = GenRatioFollower;
        ParamData.PopSizeLeader = PopSizeLeader;
        ParamData.MaxGen = MaxGen;
        ParamData.LimitsLeader = LimitsLeader;
        ParamData.PopSizeFollower = PopSizeFollower;
        ParamData.LimitsFollower = LimitsFollower;
        ParamData.CrossOverRate = CrossOverRate;
        ParamData.MutationRate = MutationRate;
        ParamData.NumRuns = NumRuns;
        ParamData.RobustType = RobustType;
        ParamData.FeffType = FeffType; % 1- mean, 2- max, 3- min, otherwise- no uncer
        ParamData.VarType = VarType; % 1- Deb and Gupta 2005, 2- Gaspar-Cunha 2014, 3- Diag, 4- no uncer
        ParamData.FeffTypeFollower = FeffTypeFollower; % 1- mean, 2- max, 3- min, otherwise- no uncer
        ParamData.VarTypeFollower = VarTypeFollower; % 1- Deb and Gupta 2005, 2- Gaspar-Cunha 2014, 3- Diag, 4- no uncer
        ParamData.eta = eta;
        ParamData.PenValue = PenValue;
        ParamData.elitePer = elitePer;

        for RunIdx = 1 : NumRuns
            %% First Generation


            DataGen(NpopTosave).PopulationLeader = [];
            DataGen(NpopTosave).NonDominatedLeader = [];
            DataGen(NpopTosave).PopulationFollower = [];
            NonDominatedLeader = [];
            MatingPoolLeader = [];
            OffspringsLeader = [];
            ParallelCPU = 0;
            NonDominatedValues = [];

            PopulationLeader = InitPop(PopSizeLeader,LimitsLeader);
            PopulationFollower = InitPop(PopSizeFollower,LimitsFollower);

            %% Leader
            PopulationLeader = EvaluationPop(PopulationLeader,ProblemName);
            PopulationLeader = EvaluateRobussGame(PopulationLeader,PopulationFollower,FeffType,VarType,ProblemName,LimitsLeader);
            PopulationLeader = ApplyRobuss(PopulationLeader,RobustType,eta,PenValue);

            GenCounter = 1;
            FitnessValues = [PopulationLeader(:).F]';
            [FrontNo,~] = NDSort(FitnessValues,1);
            NonDominatedLeader = PopulationLeader(FrontNo==1);
            % Update follower
            NonDominatedValues = [NonDominatedLeader(:).F]';
            [~,ia,~] = unique(NonDominatedValues,'rows');
            NonDominatedLeader = NonDominatedLeader(ia);
            NonDominatedValues = [NonDominatedLeader(:).F]';
            PopulationFollower = EvaluateFollower(PopulationFollower,NonDominatedLeader,FeffTypeFollower,VarTypeFollower,ProblemName,LimitsLeader);
            PopulationFollower = ApplyRobuss(PopulationFollower,RobustType,eta,PenValue);

            if PlotFlag
                if size(NonDominatedValues,2) == 3
                    figure(1)
                    plot3(NonDominatedValues(:,1),NonDominatedValues(:,2),NonDominatedValues(:,3),'xr')
                    title(['Generation number  ', num2str(GenCounter)])
                    grid on
                    drawnow

                    % figure(2)
                    % parallelcoords([NonDominated.X]')
                    % title(['Generation number  ', num2str(GenCounter)])
                    % grid on
                    % drawnow

                else

                    figure(1)
                    plot(NonDominatedValues(:,1),NonDominatedValues(:,2),'xr')
                    title(['Generation number  ', num2str(GenCounter)])
                    grid on
                    drawnow

                    % figure(2)
                    % parallelcoords([NonDominated.X]')
                    % title(['Generation number  ', num2str(GenCounter)])
                    % grid on
                    % drawnow
                end
            end
            % end
            %% Gen Loop
            while GenCounter <= MaxGen
                % Leader
                PopulationLeader = CalcRankAndDistance(PopulationLeader);
                MatingPoolLeader = SelectionSol(PopulationLeader,PopSizeLeader);
                OffspringsLeader = Reproduction(MatingPoolLeader,CrossOverRate,MutationRate,LimitsLeader);
                OffspringsLeader = EvaluationPop(OffspringsLeader,ProblemName);
                OffspringsLeader = EvaluateRobussGame(OffspringsLeader,PopulationFollower,FeffType,VarType,ProblemName,LimitsLeader);
                OffspringsLeader = ApplyRobuss(OffspringsLeader,RobustType,eta,PenValue);
                NonDominatedLeader = EvaluationPop(NonDominatedLeader,ProblemName);
                NonDominatedLeader = EvaluateRobussGame(NonDominatedLeader,PopulationFollower,FeffType,VarType,ProblemName,LimitsLeader);
                NonDominatedLeader = ApplyRobuss(NonDominatedLeader,RobustType,eta,PenValue);

                PopulationLeader = EliteFullSorting(NonDominatedLeader,OffspringsLeader,PopSizeLeader,elitePer);
                FitnessValues = [PopulationLeader(:).F]';
                [FrontNo,~] = NDSort(FitnessValues,1);
                NonDominatedLeader = PopulationLeader(FrontNo==1);
                NonDominatedValues = [NonDominatedLeader(:).F]';
                [~,ia,~] = unique(NonDominatedValues,'rows');
                NonDominatedLeader = NonDominatedLeader(ia);

                % follower
                if rem(GenCounter,GenRatioFollower) == 0 || GenCounter == 1
                    PopulationFollower = EvaluateFollower(PopulationFollower,NonDominatedLeader,FeffTypeFollower,VarTypeFollower,ProblemName,LimitsLeader);
                    PopulationFollower = ApplyRobuss(PopulationFollower,RobustType,eta,PenValue);
                    PopulationFollower = CalcRankAndDistance(PopulationFollower);
                    MatingPoolFollower = SelectionSol(PopulationFollower,PopSizeFollower);
                    OffspringsFollower = Reproduction(MatingPoolFollower,CrossOverRate,MutationRate,LimitsFollower);
                    OffspringsFollower = EvaluateFollower(OffspringsFollower,NonDominatedLeader,FeffTypeFollower,VarTypeFollower,ProblemName,LimitsLeader);
                    OffspringsFollower = ApplyRobuss(OffspringsFollower,RobustType,eta,PenValue);

                    PopulationFollower = EliteFullSorting(PopulationFollower,OffspringsFollower,PopSizeFollower,elitePer);
                    % PopulationFollower = OffspringsFollower;
                end
                if rem(GenCounter, saveEachGen) == 0
                    DataGen(GenCounter/saveEachGen).PopulationLeader = PopulationLeader;
                    DataGen(GenCounter/saveEachGen).NonDominatedLeader = NonDominatedLeader;
                    DataGen(GenCounter/saveEachGen).PopulationFollower = PopulationFollower;
                    DataGen(GenCounter/saveEachGen).GenCounter = GenCounter;
                end

                GenCounter = GenCounter + 1;
                NonDominatedValues = [NonDominatedLeader(:).F]';
                if PlotFlag && rem(GenCounter,20) == 0
                    if size(NonDominatedValues,2) == 3
                        figure(1)
                        plot3(NonDominatedValues(:,1),NonDominatedValues(:,2),NonDominatedValues(:,3),'xr',FitnessValues(:,1),FitnessValues(:,2),FitnessValues(:,3),'.k')
                        title(['Generation number  ', num2str(GenCounter)])
                        grid on
                        drawnow

                        % figure(2)
                        % parallelcoords([NonDominatedLeader.X]')
                        % title(['Generation number  ', num2str(GenCounter)])
                        % grid on
                        % drawnow

                    else

                        figure(1)
                        plot(NonDominatedValues(:,1),NonDominatedValues(:,2),'xr',FitnessValues(:,1),FitnessValues(:,2),'.k')
                        title(['Generation number  ', num2str(GenCounter)])
                        grid on
                        drawnow

                        % figure(2)
                        % parallelcoords([NonDominatedLeader.X]')
                        % title(['Generation number  ', num2str(GenCounter)])
                        % grid on
                        % drawnow
                    end
                end
                disp(['Problem ',ProblemName,' Run ', num2str(RunIdx),'\',num2str(NumRuns),' Generation ', num2str(GenCounter),'\',num2str(MaxGen)])
            end

            save([FolderName,'\GameBasedRNSGAII_',ProblemName, '_RUN_', num2str(RunIdx)], "DataGen","ParamData")
            disp(['----------Problem ',ProblemName,' Run ', num2str(RunIdx),'\',num2str(NumRuns),' Done!!!------'])
        end
        disp('=================================================')
        disp(['=============== Problem ',ProblemName,' Done!!! ============='])
        disp('=================================================')
        %%
    end


    if NumRuns == 1 && PlotEnd == 1
        PopulationFollower = EvaluateFollower(PopulationFollower,NonDominatedLeader,FeffTypeFollower,VarTypeFollower,ProblemName,LimitsLeader);
        PopulationFollower = ApplyRobuss(PopulationFollower,RobustType,eta,PenValue);

        NonDominatedLeader = EvaluationPop(NonDominatedLeader,ProblemName);
        NonDominatedFValues = [NonDominatedLeader(:).F]';
        NonDominatedLeader = EvaluateRobussGame(NonDominatedLeader,PopulationFollower,FeffType,VarType,ProblemName,LimitsLeader);
        NonDominatedLeader = ApplyRobuss(NonDominatedLeader,RobustType,eta,PenValue);
        NonDominatedValues = [NonDominatedLeader(:).F]';
        if size(NonDominatedValues,2) == 3
            figure(1)
            plot3(NonDominatedValues(:,1),NonDominatedValues(:,2),NonDominatedValues(:,3),'xr')
            title("Final O.S.")
            xlabel("F_1")
            ylabel("F_2")
            zlabel("Var")
            grid on
            drawnow
            figure(11)
            plot(NonDominatedValues(:,1),NonDominatedValues(:,2),'xr'...
                ,NonDominatedFValues(:,1),NonDominatedFValues(:,2),'oc')
            title("Final O.S.")
            legend('Feff', 'Fnom')
            xlabel("F_1")
            ylabel("F_2")
            grid on
            drawnow
            figure(2)
            parallelcoords([NonDominatedLeader.X]')
            title("Final Leader D.S.")
            grid on
            drawnow

            figure(3)
            parallelcoords([PopulationFollower.X]')
            title("Final Follower D.S.")
            grid on
            drawnow

        else

            figure(1)
            plot(NonDominatedValues(:,1),NonDominatedValues(:,2),'xr')
            title("Final O.S.")
            grid on
            drawnow
            figure(11)
            plot(NonDominatedValues(:,1),NonDominatedValues(:,2),'xr'...
                ,NonDominatedFValues(:,1),NonDominatedFValues(:,2),'oc')
            title("Final O.S.")
            legend('Feff', 'Fnom')
            drawnow
            figure(2)
            parallelcoords([NonDominatedLeader.X]')
            title("Final Leader D.S.")
            grid on
            drawnow

            figure(3)
            parallelcoords([PopulationFollower.X]')
            title("Final Follower D.S.")
            grid on
            drawnow
        end

    end
end