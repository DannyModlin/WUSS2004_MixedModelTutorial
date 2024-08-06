proc mixed data=aglm.grass ;
   class method variety;
   model yield=method / ddfm=kr2;
   random variety method*variety ;
   contrast 'A vs B and C' method 2 -1 -1;
   estimate 'A vs B and C' method 2 -1 -1 / divisor=2 cl alpha=0.02;
   estimate 'Method A mean' intercept 1 method 1 0 0;
run;	 			 	                 
