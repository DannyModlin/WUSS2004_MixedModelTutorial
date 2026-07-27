proc sgplot data=aglm.grass;
   vline variety / group=method stat=mean response=yield;
   title 'Yield vs. Variety for Each Method';
run;
title;
