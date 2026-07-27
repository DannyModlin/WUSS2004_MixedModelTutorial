/* cap input rows for the captured run */
options obs=100;

/* data setup: aglm.toy from 1_DataCreation.sas */
libname aglm (work);
data aglm.toy;     
   length toy $1;
   input toy $ adhesive $ pressure @@;
datalines;
1 c 67.0 1 b 71.9 1 a 72.2
2 c 67.5 2 b 68.8 2 a 66.4
3 c 76.0 3 b 82.6 3 a 74.5
4 c 72.7 4 b 78.1 4 a 67.3
5 c 73.1 5 b 74.2 5 a 73.2
6 c 65.8 6 b 70.8 6 a 68.7
7 c 75.6 7 b 84.9 7 a 69.0
;
run;
proc sort data=aglm.toy;
  by adhesive toy;
run;
