*******************************************************************************
*The SAS program (HEI-2015-PerDay-ASA24-2016.sas):                            *
*									      *
*This SAS program can be used to calculate Healthy Eating Index (HEI)-2015    *
*scores from 24-hour recall data collected using ASA24-2016. This program     *
*calculates HEI-2015 component and total scores for each recall (i.e., each   *
*individual's reported 1-day intake).  Additional code that calculates        *
*HEI-2015 component and total scores for multiple 24HRs is available on the   * 
*ASA24 HEI Resources page.		                                      *
*									      *
*This program has been tested using SAS, version 9.4 and uses the 'Totals'    *
*analysis file from ASA24-2016.  These program files can be downloaded from   *
*the Researcher Web site. The data file should be in CSV format.	      *
*									      *
*Note: Some users have found that the SAS program will drop observations from *
*the analysis if the ID field is not the same length for all observations.  To* 
*prevent this error, the observations with the longest ID length should be    *
*listed first when the data is imported into SAS.			      *
*								              *
*				            			              *
*									      *
*Please see accompanying readme file.   				      *
*******************************************************************************;

/* Filename for the input dataset */
filename Totals 'C:/Totals.csv'; /*Daily Total Nutrient and Pyramid Equivalents file from ASA24-2016*/

/*Name and location of the output file to be exported, containing HEI-2015 component and total scores for each intake day*/

filename res 'C:/hei2015.perday.result.csv';


/*Read in required macro*/

%include 'C:/hei2015.score.macro.sas';


TITLE 'ASA24-2016 HEI-2015 scores - by person per day';


/*Step 1.
Input daily total data and create five additional required variables.  These variables are:  
FWHOLEFRT, MONOPOLY, VTOTALLEG, VDRKGRLEG, PFALLPROTLEG, and PFSEAPLANTLEG
*/

Proc import datafile=Totals
  Out=Totals
  Dbms=csv
  Replace;
  Getnames=yes;
Run;



DATA Totals;
  SET Totals;

  FWHOLEFRT=F_CITMLB+F_OTHER;

  MONOPOLY=MFAT+PFAT;

  VTOTALLEG=V_TOTAL+V_LEGUMES;
  VDRKGRLEG=V_DRKGR+V_LEGUMES;

  PFALLPROTLEG=PF_MPS_TOTAL+PF_EGGS+PF_NUTSDS+PF_SOY+PF_LEGUMES; 
  PFSEAPLANTLEG=PF_SEAFD_HI+PF_SEAFD_LOW+PF_NUTSDS+PF_SOY+PF_LEGUMES;
  
run; 


/*Step 2. 
 Runs the HEI2015 scoring macro which calculates intake density amounts and HEI scores.
*/

%HEI2015 (indat=Totals,
          kcal= KCAL,
	  vtotalleg= VTOTALLEG,
	  vdrkgrleg= VDRKGRLEG,
	  f_total= F_TOTAL,
	  fwholefrt=FWHOLEFRT,
	  g_whole= G_WHOLE,
	  d_total= D_TOTAL,
          pfallprotleg= PFALLPROTLEG,
	  pfseaplantleg= PFSEAPLANTLEG,
	  monopoly=MONOPOLY,
	  satfat=SFAT,
	  sodium=SODI,
	  g_refined=G_REFINED,
	  add_sugars=ADD_SUGARS,
	  outdat=hei2015);
 
run;

/*Step 3.
 Displays and saves the results.
*/ 

Data hei2015r (keep=UserName UserID RecallNo kcal HEI2015C1_TOTALVEG HEI2015C2_GREEN_AND_BEAN HEI2015C3_TOTALFRUIT
      HEI2015C4_WHOLEFRUIT HEI2015C5_WHOLEGRAIN HEI2015C6_TOTALDAIRY HEI2015C7_TOTPROT HEI2015C8_SEAPLANT_PROT 
      HEI2015C9_FATTYACID HEI2015C10_SODIUM HEI2015C11_REFINEDGRAIN HEI2015C12_SFAT HEI2015C13_ADDSUG HEI2015_TOTAL_SCORE);
  Set hei2015;
  Run;

proc means n nmiss min max mean data=hei2015r;
run;

proc export data= hei2015r
  file=res
  dbms=csv
  replace;
run;


