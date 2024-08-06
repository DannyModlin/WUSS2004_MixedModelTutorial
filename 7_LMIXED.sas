libname cat '.';
%let narray  = 6;
%let ndye    = 2;
%let nrow    = 4;
%let ngene   = 500;
%let ntrt    = 6;
%let npin    = 4;
%let ndip    = 4;
%let no      = %eval(&ndye*&nrow*&ngene);
%let tno     = %eval(&narray*&no);

data cat.microarray;
	keep Gene MArray Dye Trt Pin Dip log2i;
	array PinDist{&tno};
	array DipDist{&tno};
	array GeneDist{&tno};
	array ArrayEffect{&narray};
	array ArrayGeneEffect{%eval(&narray*&ngene)};
	array ArrayDipEffect{%eval(&narray*&ndip)};
	array ArrayPinEffect{%eval(&narray*&npin)};

	do i=1 to &tno;
		PinDist{i}=1 + int(&npin*ranuni(12345));
		DipDist{i}=1 + int(&ndip*ranuni(12345));
		GeneDist{i}=1 + int(&ngene*ranuni(12345));
	end;
	igene=0;
	idip=0;
	ipin=0;

	do i=1 to &narray;
		ArrayEffect{i}=sqrt(0.014)*rannor(12345);

		do j=1 to &ngene;
			igene=igene+1;
			ArrayGeneEffect{igene}=sqrt(0.0017)*rannor(12345);
		end;

		do j=1 to &ndip;
			idip=idip + 1;
			ArrayDipEffect{idip}=sqrt(0.0033)*rannor(12345);
		end;

		do j=1 to &npin;
			ipin=ipin + 1;
			ArrayPinEffect{ipin}=sqrt(0.037)*rannor(12345);
		end;
	end;
	i=0;

	do MArray=1 to &narray;

		do Dye=1 to &ndye;

			do Row=1 to &nrow;

				do k=1 to &ngene;

					if MArray=1 and Dye=1 then
						do;
							Trt=0;
							trtc=0;
						end;
					else
						do;

							if trtc >=&no then
								trtc=0;

							if trtc=0 then
								do;
									Trt=Trt + 1;

									if Trt >=&ntrt then
										do;
											Trt=0;
											trtc=0;
										end;
								end;
							trtc=trtc + 1;
						end;
					i=i + 1;
					Pin=PinDist{i};
					Dip=DipDist{i};
					Gene=GeneDist{i};
					a=ArrayEffect{MArray};
					ag=ArrayGeneEffect{(MArray-1)*&ngene+Gene};
					ad=ArrayDipEffect{(MArray-1)*&ndip+Dip};
					ap=ArrayPinEffect{(MArray-1)*&npin+Pin};
					log2i=1 +
                        + Dye
                        + Trt
                        + Gene/1000.0
                        + Dye*Gene/1000.0
                        + Trt*Gene/1000.0
                        + Pin
                        + a
                        + ag
                        + ad
                        + ap
                        + sqrt(0.02)*rannor(12345);
					output;
				end;
			end;
		end;
	end;
run;

cas mySession;
libname Caslib cas;

proc casutil;
	load data=cat.microarray outcaslib="casuser" casout="microArray" replace;
run;

***Sparse data: takes 1 hour 30 minutes;

proc mixed data=microarray;
	*dmmethod=sparse;
	class marray dye trt gene pin dip;
	model log2i=dye trt gene dye*gene trt*gene pin;
	random int gene dip pin/subject=marray s;
	* ods output solutionr=BLUPs;
run;

***Sparse data: takes 53 seconds;

proc lmixed data=caslib.microarray dmmethod=sparse;
	class marray dye trt gene pin dip;
	model log2i=dye trt gene dye*gene trt*gene pin;
	random int gene dip pin/subject=marray;
	* ods output solutionr=BLUPs;
run;

*****************************************************;
***Secondary Example*********************************;
*****************************************************;

%let NClinic  = 100;
%let NPatient = %eval(&NClinic*50);
%let NTime    = 3;
%let SigmaC   = 2.0;
%let SigmaP   = 4.0;
%let SigmaE   = 8.0;
%let Seed     = 12345;
libname cat '.';

data WeekSim;
	keep Gender Clinic Patient Time Measurement;
	array PGender{&NPatient};
	array PClinic{&NPatient};
	array PEffect{&NPatient};
	array CEffect{&NClinic};
	array GEffect{2};

	do Clinic=1 to &NClinic;
		CEffect{Clinic}=sqrt(&SigmaC)*rannor(&Seed);
	end;
	GEffect{1}=10*ranuni(&Seed);
	GEffect{2}=10*ranuni(&Seed);

	do Patient=1 to &NPatient;
		PGender{Patient}=1 + int(2 *ranuni(&Seed));
		PClinic{Patient}=1 + int(&NClinic*ranuni(&Seed));
		PEffect{Patient}=sqrt(&SigmaP)*rannor(&Seed);
	end;

	do Patient=1 to &NPatient;
		Gender=PGender{Patient};
		Clinic=PClinic{Patient};
		Mean=1 + GEffect{Gender} + CEffect{Clinic} + PEffect{Patient};

		do Time=1 to &nTime;
			Measurement=Mean + sqrt(&SigmaE)*rannor(&Seed);
			output;
		end;
	end;
run;

cas mySession;
libname Casli cas;

proc casutil;
	load data=cat.weekSim outcaslib="casuser" casout="weekSim" replace;
run;

***runs for 20 seconds ;
proc mixed data=cat.WeekSim;
	class Gender Clinic Patient Time;
	model Measurement=Gender;
	random Clinic / s;
	repeated Time / sub=Patient type=un;
run;

***runs for 41 seconds;
proc lmixed data=caslib.WeekSim;
	class Gender Clinic Patient Time;
	model Measurement=Gender;
	random Clinic / s;
	repeated Time / sub=Patient type=un;
run;