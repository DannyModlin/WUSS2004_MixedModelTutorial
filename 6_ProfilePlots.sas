proc format;  
  value $drug 
         'a'='a-Standard'
		 'c'='c-New'
		 'p'='p-Placebo'
;
proc print data=aglm.fev1mult; 
   format drug $drug.;
proc print data=aglm.fev1uni; 
   format drug $drug.;
run;

proc sgplot data=aglm.fev1uni;
   vline hour / response=fev1 stat=mean group=drug ;
   format drug $drug.;
title 'Average FEV1 vs. Hour by Drug';
run;
title;


proc mixed data=aglm.fev1uni;
   class drug patient hour;
   model fev1=drug basefev1 drug*basefev1 hour drug*hour / ddfm=kr2;
   repeated hour / type=TOEP subject=patient(drug);
run;    

*** comparing mean models using maximum likelihood;
proc mixed data=aglm.fev1uni method=ml;
   class drug patient hour;
   model fev1=drug basefev1 drug*basefev1 hour drug*hour / ddfm=kr2 ;
   repeated hour / type=TOEP subject=patient(drug);
   ods output FitStatistics=FitFull(rename=(value=Full));
   title 'Full Model Using ML';

run;     

proc mixed data=aglm.fev1uni method=ml;
   class drug patient hour;
   model fev1=drug basefev1 hour drug*hour / ddfm=kr2 ;
   repeated hour / type=TOEP subject=patient(drug);
   ods output FitStatistics=FitReduced(rename=(value=Reduced));
   title 'Reduced Model Using ML';
run;    
title; 

data fit;
   merge FitFull FitReduced;
run;

proc print data=fit;
run;  

*** Fit the final model using REML;
proc mixed data=aglm.fev1uni;
   class drug patient hour;
   model fev1=drug basefev1 hour drug*hour / ddfm=kr2;
   repeated hour / type=TOEP subject=patient(drug);
   title 'Final Model Using REML';
run;   
