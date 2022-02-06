function []=simSA(seedid)
% The function to simulate the SA algorithm
warning('off','all')
rng(seedid)

for nsa = 1:50
% nsa: How many times simulate SA code(50 times)
    ddoutputnow = inf(1,9);
    while (ddoutputnow(5) > 1 || ddoutputnow(6) > 1 || ddoutputnow(7) > 1 || ddoutputnow(8) > 1) 
    % If fitting costs does not become lowered than their thresholds, then do not go to LD fitting and repeat DD fitting again
        while ddoutputnow(2) > 1
            % Select intial parameter until the fluctuation cost becomes lowered than threshold
            ddnow=[10^(5*rand-4),1,10^(2*rand-1),10^(2*rand-1),10^(5*rand-4),10^(log10(5)*rand),10^(2*rand-1),10^(2*rand-1),10^(2*rand-1),10^(2*rand-1),10^(2*rand-1),10^(2*rand-1),10^(4*rand-5),10^(4*rand-5),10^(4*rand-5),10^(4*rand-5),10^(4*rand),10^(4*rand),10^(4*rand)];ddnow(19)=ddnow(18);            
            [ddcostnow,fluccostnow,ddfitcostnow,dddiffcostnow,ampcostnow,percostnow,levcostnow,relampcostnow,tranTcostnow,propertyampnow,propertyrelampnow,propertyperiodnow]=newddsa(ddnow(1),1,ddnow(3),ddnow(4),ddnow(5),ddnow(6),ddnow(7),ddnow(8),ddnow(9),ddnow(10),ddnow(11),ddnow(12),ddnow(13),ddnow(14),ddnow(15),ddnow(16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19));
            ddoutputnow=[ddcostnow,fluccostnow,ddfitcostnow,dddiffcostnow,ampcostnow,percostnow,levcostnow,relampcostnow,tranTcostnow];
            ddpropertiesnow=[propertyampnow,propertyrelampnow,propertyperiodnow];
        end
        
        result=[];
        disp('SA Start');
        tempdd = 1;num = 1;
        while tempdd > 0.006 && ( ddoutputnow(5) > 1 || ddoutputnow(6) > 1 || ddoutputnow(7) > 1 || ddoutputnow(8) > 1) 
        % Terminate iteration if it iterates about 10000 times, or fitting costs become lowered than their thresholds
            tempdd = 0.9995*tempdd;num=num+1;
            ddprev=[ddnow(1)*exp(0.4*rand-0.2),1,ddnow(3)*exp(0.4*rand-0.2),ddnow(4)*exp(0.4*rand-0.2),ddnow(5)*exp(0.4*rand-0.2),max(1,ddnow(6)*exp(0.4*rand-0.2)),...
                ddnow(7)*exp(0.4*rand-0.2),ddnow(8)*exp(0.4*rand-0.2),ddnow(9)*exp(0.4*rand-0.2),ddnow(10)*exp(0.4*rand-0.2),ddnow(11)*exp(0.4*rand-0.2),max(0.1,min(10,ddnow(12)*exp(0.4*rand-0.2))),...
                ddnow(13)*exp(rand-0.5),ddnow(14)*exp(rand-0.5),ddnow(15)*exp(rand-0.5),ddnow(16)*exp(rand-0.5),max(1,ddnow(17)*exp(rand-0.5)),max(1,ddnow(18)*exp(rand-0.5)),ddnow(19)];
            ddprev(19)=ddprev(18);
            % Perturb parameters
            
            [ddcostprev,fluccostprev,ddfitcostprev,dddiffcostprev,ampcostprev,percostprev,levcostprev,relampcostprev,tranTcostprev,propertyampprev,propertyrelampprev,propertyperiodprev]=newddsa(ddprev(1),ddprev(2),ddprev(3),ddprev(4),ddprev(5),ddprev(6),ddprev(7),ddprev(8),ddprev(9),ddprev(10),ddprev(11),ddprev(12),ddprev(13),ddprev(14),ddprev(15),ddprev(16),ddprev(14)*ddprev(17),ddprev(15)*ddprev(18),ddprev(16)*ddprev(19));
            ddoutputprev=[ddcostprev,fluccostprev,ddfitcostprev,dddiffcostprev,ampcostprev,percostprev,levcostprev,relampcostprev,tranTcostprev];
            ddpropertiesprev=[propertyampprev,propertyrelampprev,propertyperiodprev];
            
            if ddoutputnow(1) > ddoutputprev(1)
                ddnow=ddprev;ddoutputnow=ddoutputprev;ddpropertiesnow=ddpropertiesprev;
                % Accept parameters if the total cost is reduced
            else
                if 1/(exp(abs(ddoutputnow(1) - ddoutputprev(1))/tempdd)) > rand
                    ddnow=ddprev;ddoutputnow=ddoutputprev;ddpropertiesnow=ddpropertiesprev;
                    % Accept parameter randomly
                end
            end
            if mod(num,500)==1
                result=[result; [num,100,ddoutputnow,zeros(1,9),...
                    ddnow(1:16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
                    zeros(1,10),ddpropertiesnow,zeros(1,20)]];
                dlmwrite(strcat(num2str(nsa),'Result', num2str(seedid), '.csv'),result,'precision',20);
                % For every 500 iterations, record the result to csv file
            end
        end
    end
    
    dlmwrite(strcat(num2str(nsa),'Result', num2str(seedid), '.csv'),result,'precision',20);
    
    [~,~,periodlnWTdd,~,~,~,~]=nddmeasure(ddnow(1),1,1,ddnow(3),ddnow(4),1,ddnow(12),ddnow(13),ddnow(14),ddnow(15),ddnow(16));
    numEnt = 1;tempEnt=1;iterEnt=1;
    ldnow=[0.1*rand+0.95,ddnow(4),(2.9*rand+1.1)*ddnow(4),ddnow(10),(2.9*rand+1.1)*ddnow(10),0.19*rand+0.01,0.19*rand+0.01,0.19*rand+0.01,0.19*rand+0.01,rand/12];
    while goodL(ldnow(6),ldnow(7))==0
        ldnow(6)=0.19*rand+0.01;ldnow(7)=0.19*rand+0.01;
    end
    while goodL(ldnow(8),ldnow(9))==0
        ldnow(8)=0.19*rand+0.01;ldnow(9)=0.19*rand+0.01;
    end
    periodnow=ldnow(1)*periodlnWTdd;
    [ldcostnow,stabcostnow,ldfitcostnow,lddiffcostnow,lampcostnow,lpercostnow,lrelampcostnow,lphasecostnow,ltranTcostnow,propertylampnow,propertylrelampnow,propertylperiodnow,propertylpeakphasenow,propertylthphasenow]=...
        newldsa(ddnow(1),1,ddnow(3),ddnow(4),ddnow(5),ddnow(6),ddnow(7),ddnow(8),ddnow(9),ddnow(10),ddnow(11),ddnow(12),ddnow(13),ddnow(14),ddnow(15),ddnow(16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
        periodnow,ldnow(2),ldnow(3)-ldnow(2),ldnow(4),ldnow(5)-ldnow(4),ldnow(6)*periodnow,ldnow(7)*periodnow,ldnow(8)*periodnow,ldnow(9)*periodnow,ldnow(10)*periodnow);
    ldoutputnow=[ldcostnow,stabcostnow,ldfitcostnow,lddiffcostnow,lampcostnow,lpercostnow,lrelampcostnow,lphasecostnow,ltranTcostnow];
    ldpropertiesnow=[propertylampnow,propertylrelampnow,propertylperiodnow,propertylpeakphasenow,propertylthphasenow];
    disp([ddoutputnow ldoutputnow]);
    
    while ldoutputnow(2) > 1 && iterEnt <= 1000
        % Select intial parameter until the entrain cost becomes lowered than threshold
        ldnow=[0.1*rand+0.95,ddnow(4),(2.9*rand+1.1)*ddnow(4),ddnow(10),(2.9*rand+1.1)*ddnow(10),0.19*rand+0.01,0.19*rand+0.01,0.19*rand+0.01,0.19*rand+0.01,rand/12];
        while goodL(ldnow(6),ldnow(7))==0
            ldnow(6)=0.19*rand+0.01;ldnow(7)=0.19*rand+0.01;
        end
        while goodL(ldnow(8),ldnow(9))==0
            ldnow(8)=0.19*rand+0.01;ldnow(9)=0.19*rand+0.01;
        end
        periodnow=ldnow(1)*periodlnWTdd;
        [ldcostnow,stabcostnow,ldfitcostnow,lddiffcostnow,lampcostnow,lpercostnow,lrelampcostnow,lphasecostnow,ltranTcostnow,propertylampnow,propertylrelampnow,propertylperiodnow,propertylpeakphasenow,propertylthphasenow]=...
            newldsa(ddnow(1),1,ddnow(3),ddnow(4),ddnow(5),ddnow(6),ddnow(7),ddnow(8),ddnow(9),ddnow(10),ddnow(11),ddnow(12),ddnow(13),ddnow(14),ddnow(15),ddnow(16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
            periodnow,ldnow(2),ldnow(3)-ldnow(2),ldnow(4),ldnow(5)-ldnow(4),ldnow(6)*periodnow,ldnow(7)*periodnow,ldnow(8)*periodnow,ldnow(9)*periodnow,ldnow(10)*periodnow);
        ldoutputnow=[ldcostnow,stabcostnow,ldfitcostnow,lddiffcostnow,lampcostnow,lpercostnow,lrelampcostnow,lphasecostnow,ltranTcostnow];
        ldpropertiesnow=[propertylampnow,propertylrelampnow,propertylperiodnow,propertylpeakphasenow,propertylthphasenow];
        disp([ddoutputnow ldoutputnow]);
        iterEnt=iterEnt+1;
    end
    result=[result; [0,ddoutputnow(1)+ldoutputnow(1),ddoutputnow,ldoutputnow,...
        ddnow(1:16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
        ldnow(1:5),ldnow(6)*periodnow,ldnow(7)*periodnow,ldnow(8)*periodnow,ldnow(9)*periodnow,ldnow(10)*periodnow,ddpropertiesnow,ldpropertiesnow]];
    dlmwrite(strcat(num2str(nsa),'Result', num2str(seedid), '.csv'),result,'precision',20);
     % Record the result to csv file
    
    disp('LD SA Start');numld=1;temp=1;
    while temp > 0.001 &&  (ddoutputnow(5) > 1 || ddoutputnow(6) > 1 || ddoutputnow(7) > 1 || ddoutputnow(7) > 1 || ddoutputnow(8) > 1 ||ldoutputnow(2) > 1|| ldoutputnow(5) > 1 || ldoutputnow(6) > 1 || ldoutputnow(7) > 1 || ldoutputnow(8) > 1 ...
        ||isrhythmM(ddnow(1),ddnow(2),1,ddnow(3),ddnow(4),1,ddnow(12),ddnow(13),ddnow(14),ddnow(15),ddnow(16),periodnow,ldnow(3)-ldnow(2),ldnow(2),ldnow(6),ldnow(7),ldnow(10)) ...
        ||isrhythmM(ddnow(5),ddnow(7),ddnow(8),ddnow(9),ddnow(10),ddnow(11),ddnow(12),ddnow(13),ddnow(14),ddnow(15),ddnow(16),periodnow,ldnow(5)-ldnow(4),ldnow(4),ldnow(8),ldnow(9),ldnow(10)))
        % Terminate iteration if it iterates about 10000 times, or fitting costs become lowered than their thresholds
        numld=numld+1;temp = 0.9995*temp;
        ddprev=[ddnow(1)*exp(0.1*rand-0.05),1,ddnow(3)*exp(0.1*rand-0.05),ddnow(4)*exp(0.1*rand-0.05),ddnow(5)*exp(0.1*rand-0.05),max(1,ddnow(6)*exp(0.1*rand-0.05)),...
            ddnow(7)*exp(0.1*rand-0.05),ddnow(8)*exp(0.1*rand-0.05),ddnow(9)*exp(0.1*rand-0.05),ddnow(10)*exp(0.1*rand-0.05),ddnow(11)*exp(0.1*rand-0.05),min(10,max(0.1,ddnow(12)*exp(0.1*rand-0.05))),...
            ddnow(13)*exp(0.2*rand-0.1),ddnow(14)*exp(0.2*rand-0.1),ddnow(15)*exp(0.2*rand-0.1),ddnow(16)*exp(0.2*rand-0.1),max(1,ddnow(17)*10^(0.2*rand-0.1)),max(1,ddnow(18)*10^(0.2*rand-0.1)),ddnow(19)];
        ddprev(19)=ddprev(18);
        ldprev=[min(max(0.95,ldnow(1)*exp(0.4*rand-0.2)),1.05),ldnow(2)*exp(0.4*rand-0.2),ldnow(3)*exp(0.4*rand-0.2),...
            ldnow(4)*exp(0.4*rand-0.2),ldnow(5)*exp(0.4*rand-0.2),max(0.01,min(0.2,ldnow(6)*exp(0.4*rand-0.2))),...
            max(0.01,min(0.2,ldnow(7)*exp(0.4*rand-0.2))),max(0.01,min(0.2,ldnow(8)*exp(0.4*rand-0.2))),max(0.01,min(0.2,ldnow(9)*exp(0.4*rand-0.2))),min(max(0.001,ldnow(10)*exp(0.4*rand-0.2)),1/12)];
        ldprev(3)=max([ddprev(4),1.1*ldprev(2),ldprev(3)]);
        ldprev(5)=max([ddprev(10),1.1*ldprev(4),ldprev(5)]);       
        while goodL(ldprev(6),ldprev(7))==0
            ldprev(6)=max(0.01,min(0.2,ldnow(6)*exp(0.4*rand-0.2)));ldprev(7)=max(0.01,min(0.2,ldnow(7)*exp(0.4*rand-0.2)));
        end
        while goodL(ldprev(8),ldprev(9))==0
            ldprev(8)=max(0.01,min(0.2,ldnow(8)*exp(0.4*rand-0.2)));ldprev(9)=max(0.01,min(0.2,ldnow(9)*exp(0.4*rand-0.2)));
        end
        % Perturb parameters
        
        [~,~,periodlnWTdd,~,~,~,~]=nddmeasure(ddprev(1),1,1,ddprev(3),ddprev(4),1,ddprev(12),ddprev(13),ddprev(14),ddprev(15),ddprev(16));
        periodprev=ldprev(1)*periodlnWTdd;
        [ddcostprev,fluccostprev,ddfitcostprev,dddiffcostprev,ampcostprev,percostprev,levcostprev,relampcostprev,tranTcostprev,propertyampprev,propertyrelampprev,propertyperiodprev]=newddsa(ddprev(1),ddprev(2),ddprev(3),ddprev(4),ddprev(5),ddprev(6),ddprev(7),ddprev(8),ddprev(9),ddprev(10),ddprev(11),ddprev(12),ddprev(13),ddprev(14),ddprev(15),ddprev(16),ddprev(14)*ddprev(17),ddprev(15)*ddprev(18),ddprev(16)*ddprev(19));
        ddoutputprev=[ddcostprev,fluccostprev,ddfitcostprev,dddiffcostprev,ampcostprev,percostprev,levcostprev,relampcostprev,tranTcostprev];ddpropertiesprev=[propertyampprev,propertyrelampprev,propertyperiodprev];
        [ldcostprev,stabcostprev,ldfitcostprev,lddiffcostprev,lampcostprev,lpercostprev,lrelampcostprev,lphasecostprev,ltranTcostprev,propertylampprev,propertylrelampprev,propertylperiodprev,propertylpeakphaseprev,propertylthphaseprev]=...
            newldsa(ddprev(1),ddprev(2),ddprev(3),ddprev(4),ddprev(5),ddprev(6),ddprev(7),ddprev(8),ddprev(9),ddprev(10),ddprev(11),ddprev(12),ddprev(13),ddprev(14),ddprev(15),ddprev(16),ddprev(14)*ddprev(17),ddprev(15)*ddprev(18),ddprev(16)*ddprev(19),periodprev,ldprev(2),ldprev(3)-ldprev(2),ldprev(4),ldprev(5)-ldprev(4),ldprev(6)*periodprev,ldprev(7)*periodprev,ldprev(8)*periodprev,ldprev(9)*periodprev,ldprev(10)*periodprev);
        ldoutputprev=[ldcostprev,stabcostprev,ldfitcostprev,lddiffcostprev,lampcostprev,lpercostprev,lrelampcostprev,lphasecostprev,ltranTcostprev];
        ldpropertiesprev=[propertylampprev,propertylrelampprev,propertylperiodprev,propertylpeakphaseprev,propertylthphaseprev];
        
        % Accept parameter and export to excel file if cost is reduced
        if ddoutputnow(1)+ldoutputnow(1) > ddoutputprev(1)+ldoutputprev(1)
            ddnow=ddprev;ddoutputnow=ddoutputprev;ldnow=ldprev;ldoutputnow=ldoutputprev;ddpropertiesnow=ddpropertiesprev;ldpropertiesnow=ldpropertiesprev;
            result=[result; [numld,ddoutputnow(1)+ldoutputnow(1),ddoutputnow,ldoutputnow,...
                ddnow(1:16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
                ldnow(1:5),ldprev(6)*periodprev,ldprev(7)*periodprev,ldprev(8)*periodprev,ldprev(9)*periodprev,ldprev(10)*periodprev,ddpropertiesnow,ldpropertiesnow]];
            dlmwrite(strcat(num2str(nsa),'Result', num2str(seedid), '.csv'),result,'precision',20);
            
        else
            if 1/(exp(abs(ddoutputnow(1)+ldoutputnow(1)-ddoutputprev(1)-ldoutputprev(1))/temp)) > rand
                ddnow=ddprev;ddoutputnow=ddoutputprev;ldnow=ldprev;ldoutputnow=ldoutputprev;ddpropertiesnow=ddpropertiesprev;ldpropertiesnow=ldpropertiesprev;
                result=[result; [numld,ddoutputnow(1)+ldoutputnow(1),ddoutputnow,ldoutputnow,...
                    ddnow(1:16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
                    ldnow(1:5),ldprev(6)*periodprev,ldprev(7)*periodprev,ldprev(8)*periodprev,ldprev(9)*periodprev,ldprev(10)*periodprev,ddpropertiesnow,ldpropertiesnow]];
                dlmwrite(strcat(num2str(nsa),'Result', num2str(seedid), '.csv'),result,'precision',20);
            end
        end
    end
    result=[result; [numld,ddoutputnow(1)+ldoutputnow(1),ddoutputnow,ldoutputnow,...
        ddnow(1:16),ddnow(14)*ddnow(17),ddnow(15)*ddnow(18),ddnow(16)*ddnow(19),...
        ldnow(1:5),ldprev(6)*periodprev,ldprev(7)*periodprev,ldprev(8)*periodprev,ldprev(9)*periodprev,ldprev(10)*periodprev,ddpropertiesnow,ldpropertiesnow];zeros(1,81)];
    dlmwrite(strcat(num2str(nsa),'Result', num2str(seedid), '.csv'),result,'precision',20);
    if ddoutputnow(5) > 1 || ddoutputnow(6) > 1 ||  ddoutputnow(7) > 1 || ddoutputnow(8) > 1 ||ldoutputnow(2) > 1|| ldoutputnow(5) > 1 || ldoutputnow(6) > 1 || ldoutputnow(7) > 1 || ldoutputnow(8) > 1
        nsa=nsa-1;
    end
end
end