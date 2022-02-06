function [ddcost,fluccost,ddfitcost,dddiffcost,ampcost,percost,levcost,relampcost,tranTcost,propertyamp, propertyrelamp, propertyperiod]=newddsa(lnAT,lnK,lnb1,lnb2,dnAT,diffmutAT,dnK,dna3,dnb1,dnb2,dnb3,tranRatio,KPA,wtKA,wtKd,wtKP,mutKA,mutKd,mutKP)
   % The function calculating costs with parameters under DD
   reallnwtdd = [768869.3188, 264042.5173, 126016.8675, 728639.5924, 846200.8203, 854025.8136];
   reallnmutdd = [182167.7517, 130466.4109, 267972.8808, 248968.4884, 183000.8868, 192979.5297];
   realdnwtdd = [583159.2737, 599346.2874, 470768.6051, 509385.2408, 832937.6661, 649438.0032];
   realdnmutdd = [391891.2035,267775.6933, 316607.1311, 312674.3787, 600702.4079, 645335.4512];
   norm = max([max(reallnwtdd),max(reallnmutdd),max(realdnwtdd),max(realdnmutdd)]);
   realamplnwtdd=max(reallnwtdd/norm)-min(reallnwtdd/norm);realamplnmutdd=max(reallnmutdd/norm)-min(reallnmutdd/norm);realampdnwtdd=max(realdnwtdd/norm)-min(realdnwtdd/norm);realampdnmutdd=max(realdnmutdd/norm)-min(realdnmutdd/norm);
   realrelamplnwtdd=1-min(reallnwtdd/norm)/max(reallnwtdd/norm);realrelampdnwtdd=1-min(realdnwtdd/norm)/max(realdnwtdd/norm);realrelampdnmutdd=1-min(realdnmutdd/norm)/max(realdnmutdd/norm);
   reallevlnwtdd=(max(reallnwtdd/norm)+min(reallnwtdd/norm))/2;reallevlnmutdd=(max(reallnmutdd/norm)+min(reallnmutdd/norm))/2;reallevdnwtdd=(max(realdnwtdd/norm)+min(realdnwtdd/norm))/2;reallevdnmutdd=(max(realdnmutdd/norm)+min(realdnmutdd/norm))/2;
   
   mutlnAT=lnAT*diffmutAT;mutdnAT=dnAT*diffmutAT;
   [trscTdnMUTdd,amplitudednMUTdd,perioddnMUTdd,relampdnMUTdd,leveldnMUTdd,maxdnMUTdd,costf]=nddmeasure(mutdnAT,dnK,dna3,dnb1,dnb2,dnb3,tranRatio,KPA,mutKA,mutKd,mutKP);
   fluccost = costf;
   [trscTlnWTdd,amplitudelnWTdd,periodlnWTdd,relamplnWTdd,levellnWTdd,maxlnWTdd,costf]=nddmeasure(lnAT,lnK,1,lnb1,lnb2,1,tranRatio,KPA,wtKA,wtKd,wtKP);
   fluccost = fluccost + costf;
   [trscTdnWTdd,amplitudednWTdd,perioddnWTdd,relampdnWTdd,leveldnWTdd,maxdnWTdd,costf]=nddmeasure(dnAT,dnK,dna3,dnb1,dnb2,dnb3,tranRatio,KPA,wtKA,wtKd,wtKP);
   fluccost = (fluccost + costf);
   [trscTlnMUTdd,amplitudelnMUTdd,periodlnMUTdd,relamplnMUTdd,levellnMUTdd,maxlnMUTdd,costf]=nddmeasure(mutlnAT,lnK,1,lnb1,lnb2,1,tranRatio,KPA,mutKA,mutKd,mutKP);
   % With 'nddmeasure' function, calcuate amplitude, period, level, ... of the four models
   fluccost=fluccost/2;
                
   ampcost=5*(min(abs(1-(amplitudednWTdd/max(amplitudednMUTdd,0.0001))/(realampdnwtdd/realampdnmutdd)),1)+min(abs(1-(amplitudelnWTdd/max(amplitudednMUTdd,0.0001))/(realamplnwtdd/realampdnmutdd)),1)+abs((amplitudelnMUTdd/max(amplitudednMUTdd,0.0001)))); 
   relampcost=2.5*(min(abs(1-relamplnWTdd/realrelamplnwtdd),1)+min(abs(1-relampdnWTdd/realrelampdnwtdd),1)+min(1,abs(1-relampdnMUTdd/realrelampdnmutdd))+abs(relamplnMUTdd));
   percost=5*(min(abs(1-perioddnWTdd/max(periodlnWTdd, 0.0001)),1)+min(abs(1-perioddnMUTdd/max(perioddnWTdd, 0.0001)),1)+min(abs(1-periodlnWTdd/max(perioddnMUTdd, 0.0001)),1)); 
   levcost=2.5*(min(abs(1-(levellnMUTdd/levellnWTdd)/(reallevlnmutdd/reallevlnwtdd)),1)+min(abs(1-(leveldnWTdd/levellnWTdd)/(reallevdnwtdd/reallevlnwtdd)),1)+min(abs(1-(leveldnMUTdd/levellnWTdd)/(reallevdnmutdd/reallevlnwtdd)),1));
   tranTcost=0;
   ddfitcost=max(min(ampcost,10),1)*ampcost+max(min(relampcost,10),1)*relampcost+max(min(levcost,10),1)*levcost+max(min(percost,10),1)*percost;
   dddiffcost=tanh(abs(log10(dnAT/lnAT))/2) + tanh(abs(log10(dnK))/2) + tanh(abs(log10(dna3))/2) + tanh(abs(log10(dnb1/lnb1))/2) + tanh(abs(log10(dnb2/lnb2))/2) + tanh(abs(log10(dnb3))/2);
   ddcost=fluccost+ddfitcost+dddiffcost;
   % Calculate the fitting cost with the measured amplitudes, periods, ...

   propertyamp=[amplitudelnWTdd amplitudelnMUTdd amplitudednWTdd amplitudednMUTdd];
   propertyrelamp=[relamplnWTdd relamplnMUTdd relampdnWTdd relampdnMUTdd];
   propertyperiod=[periodlnWTdd periodlnMUTdd perioddnWTdd perioddnMUTdd];
end