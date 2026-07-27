proc sgplot data=aglm.toy;
   series x=toy y=pressure / group=adhesive;
   title 'Pressure vs Toy Number';
run;
title;

proc sgplot data=aglm.toy;
   vbox pressure / category=adhesive;
run;

proc mixed data=aglm.toy;
   class adhesive toy;
   model pressure=adhesive / e3;
   random toy;
run; 
