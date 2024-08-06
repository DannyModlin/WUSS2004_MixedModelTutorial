proc sgplot data=aglm.wafer4;
   reg y=deposit x=thick / group=temp;
run; 

proc mixed data=aglm.wafer4; 
   class temp wafer;
   model deposit=temp thick thick*temp / ddfm=kr2 solution;
   random wafer(temp) ;
run; 

proc mixed data=aglm.wafer4; 
   class temp wafer;
   model deposit=temp thick / solution ddfm=kr2;
   random wafer(temp);
run;                                         

ods select estimates;
proc mixed data=aglm.wafer4; 
   class temp wafer;
   model deposit=temp thick / solution ddfm=kr2;
   random wafer(temp);
   estimate 'Intercept for temp 900' intercept 1 temp 1 0 0;
   estimate 'Intercept for temp 1000' intercept 1 temp 0 1 0;
   estimate 'Intercept for temp 1100' intercept 1 temp 0 0 1;
run;

               

