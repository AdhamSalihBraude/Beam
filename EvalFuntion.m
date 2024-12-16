function PopObj = EvalFuntion(PopDec,FuncStr)
%Taken from PLATEMO see Copyright below - cited in the paper
PopDec = PopDec';

if strcmp(FuncStr,'TP1')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------

    PopObj(:,1) = PopDec(:,1);
    g = mean(PopDec(:,2:end),2);
    h = -((PopObj(:,1)-0.6).^3-0.4^3)/(0.6^3+0.4^3);
    s = 1./(PopDec(:,1)+0.2);
    PopObj(:,2) = h + g.*s;


elseif strcmp(FuncStr,'TP2')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = cos(pi/2*PopDec(:,1));
    g = 1 + 10*mean(PopDec(:,2:end),2);
    PopObj(:,2) = sin(pi/2*PopDec(:,1)).*g;
elseif strcmp(FuncStr,'TP3')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = 1 - PopDec(:,1).^2;
    g = 1 + 10*mean(PopDec(:,2:end),2);
    PopObj(:,2) = sin(pi/2*PopDec(:,1)).*g;
elseif strcmp(FuncStr,'TP4')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------


    PopObj(:,1) = (exp(PopDec(:,1))-1)/(exp(1)-1);
    g = 1 + 10*mean(PopDec(:,2:end),2);
    h = sin(4*pi*PopDec(:,1))/15 - PopDec(:,1) + 1;
    PopObj(:,2) = h.*g;

elseif strcmp(FuncStr,'TP5')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = PopDec(:,1);
    g = 1 + 10*mean(PopDec(:,2:end),2);
    h = sin(4*pi*PopDec(:,1))/15 - PopDec(:,1) + 1;
    PopObj(:,2) = h.*g;

elseif strcmp(FuncStr,'TP6')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = PopDec(:,1);
    h = 1 - PopDec(:,1).^2;
    g = sum(10-10*cos(4*pi*(PopDec(:,2:end)))+PopDec(:,2:end).^2,2);
    S = 1./(0.2+PopDec(:,1)) + PopDec(:,1).^2;
    PopObj(:,2) = h + g.*S;
elseif strcmp(FuncStr,'TP7')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = PopDec(:,1);
    h = 1 - PopDec(:,1).^2;
    g = sum(10-10*cos(4*pi*(PopDec(:,2:end)))+PopDec(:,2:end).^2,2);
    S = 1./(0.2+PopDec(:,1)) + 10*PopDec(:,1).^2;
    PopObj(:,2) = h + g.*S;
elseif strcmp(FuncStr,'TP8')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = PopDec(:,1);
    h = 2 - 0.8*exp(-((PopDec(:,2)-0.35)/0.25).^2) - exp(-((PopDec(:,2)-0.85)/0.03).^2);
    g = 50*sum(PopDec(:,3:end).^2,2);
    S = 1 - sqrt(PopObj(:,1));
    PopObj(:,2) = h.*(g+S);
elseif strcmp(FuncStr,'TP9')
    % <multi> <real> <large/none> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = PopDec(:,1);
    h = 2 - PopObj(:,1) - 0.8*exp(-((PopObj(:,1)+PopDec(:,2)-0.35)/0.25).^2) - exp(-((PopObj(:,1)+PopDec(:,2)-0.85)/0.03).^2);
    g = 50*sum(PopDec(:,3:end).^2,2);
    S = 1 - sqrt(PopObj(:,1));
    PopObj(:,2) = h.*(g+S);

elseif strcmp(FuncStr,'TP10')
    % <multi> <real> <large/none> <constrained> <robust>
    % Test problem for robust multi-objective optimization
    % delta --- 0.05 --- Maximum disturbance degree
    % H     ---   50 --- Number of disturbances

    %------------------------------- Reference --------------------------------
    % A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
    % analysis for multi-objective optimization: benchmark problems, Structural
    % and Multidisciplinary Optimization, 2014, 49: 771-793.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    PopObj(:,1) = PopDec(:,1).*sqrt(16+PopDec(:,3).^2) + PopDec(:,2).*sqrt(1+PopDec(:,3).^2);
    PopObj(:,2) = 20*sqrt(16+PopDec(:,3).^2)./PopDec(:,1)./PopDec(:,3);

elseif strcmp(FuncStr,'Truss')
    %from "A Multi-Objective Genetic Algorithm for Robust Design Optimization"
    L1=4000;L2=1000;P=100000;
    n=size(PopDec,2);

    A1=PopDec(:,1);A2=PopDec(:,2);y=PopDec(:,3);
    f1=A1.*sqrt(L1^2+y.^2)+A2.*sqrt(L2^2+y.^2);
    F1=(L2/(L1+L2)*P)*sqrt(L1^2+y.^2)./y;
    F2=(L1/(L1+L2)*P)*sqrt(L2^2+y.^2)./y;
    Sigma1=F1./A1*1e6;
    Sigma2=F2./A2*1e6;
    sigmaMax = 100e6;
    f2=max([Sigma1 Sigma2],[],2);
    alpha = ones(length(f2),1);
    %alpha(f2>sigmaMax) = 10^5;
    PopObj=[f1,f2].*alpha;

elseif strcmp(FuncStr,'beam')
     %%%%  from "Multiobjective robust optimization for crashworthiness design
    %%%%  of foamfilled thin-walled structures with random and interval uncertainties"
    %%%  for the beam, the extraparam=[E, P1, L] ;  X=[h,w,t,b]
    param=[600e3, 50e3, 210e9 ,0.2];%%  P1=600kN, Q=50kN, E=210GPa,   L=200cm
    P=param(1);Q=param(2);E=param(3);L=param(4);
    C=P*L^3/(48*E);
    n=size(PopDec,2);
    PopDec=PopDec/1000;
    h=PopDec(:,1);w=PopDec(:,2);t=PopDec(:,3);b=PopDec(:,4);
    Iz=(w.*h.^3-(w-t).*(h-2*b).^3)/12;
    Iy=(2*b.*w.^3+(h-2*b).*t.^3)/12;
    
    f2=2*w.*b+t.*(h-2*b);

    Sigma=0.5*L*(Q.*w./Iy+P.*h./Iz);
    f1=Sigma;
    Sigmamax=500e6;
    alpha = ones(length(f2),1);
    % alpha(Sigma>Sigmamax) = 100^5;
    PopObj=[f1/10^6,f2].*alpha;
end

PopObj = PopObj';





end