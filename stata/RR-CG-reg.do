/*   This script replicates Tables (3) and (5) of the paper:
     Caldara, Dario and Edward Herbst (2018), "Monetary Policy, Real Activity, and Credit Spreads:
     Evidence from Bayesian Proxy SVAR", American Economic Journal: Macroeconomics.
*/
  
     /* Set options: */
     #delimit;
     set more off;
     set linesize 255;
     capture log close; 
      
      
     /* Output file: */
     capture log using RR-CG.log, replace; 
     mata: mata set matafavor speed;
      
     /* Read-in the data: */
     use RR-CG-data.dta, clear;
     gen tm = ym(year,month);     

     ////////////////////////////////////////
     // Regression described in Equation (16)
     // Table (3)
     ///////////////////////////////////////

     // First Column (no corporate credit spreads)    
     ivreg2 dtarg oldtarg graym gray0 gray1 gray2 igrym igry0 igry1
     igry2 gradm grad0 grad1 grad2 igrdm igrd0 igrd1 igrd2 grau0, r bw(auto) small;
      
      
     // Output Growth
     lincom graym+gray0+gray1+gray2; 
     // Inflation
     lincom gradm+grad0+grad1+grad2; 
     // Output Growth Revision
     lincom igrym+igry0+igry1+igry2; 
     //Inflation Revision
     lincom igrdm+igrd0+igrd1+igrd2;
      
     reg dtarg oldtarg graym gray0 gray1 gray2 igrym igry0 igry1 igry2
         gradm grad0 grad1 grad2 igrdm igrd0 igrd1 igrd2 grau0, b noheader notable;
         fitstat;
      
     //Store residuals
     cap drop res_*;
     local res_count = 1;
     predict double res_`res_count', residuals;
     local ++res_count;
      
     // Second Column (with corporate credit spreads)    
     ivreg2 dtarg oldtarg graym gray0 gray1 gray2 igrym igry0 igry1
     igry2 gradm grad0 grad1 grad2 igrdm igrd0 igrd1 igrd2 grau0 baa_yld1000_avg5d, r bw(auto) small;
      
     // Output Growth
     lincom graym+gray0+gray1+gray2; 
     // Inflation
     lincom gradm+grad0+grad1+grad2; 
     // Output Growth Revision
     lincom igrym+igry0+igry1+igry2; 
     //Inflation Revision
     lincom igrdm+igrd0+igrd1+igrd2;
      
     reg dtarg oldtarg graym gray0 gray1 gray2 igrym igry0
     igry1 igry2 gradm grad0 grad1 grad2 igrdm igrd0 igrd1 igrd2 grau0 baa_yld1000_avg5d, noheader notable b;
     fitstat;

     //Store residuals
     predict double res_`res_count' , residuals;
     local ++res_count;

     // Robustness (with average corporate credit spreads calculated from the first day of the month when the FOMC
     // meeting takes place to the day prior to the meeting.)

     ivreg2 dtarg oldtarg graym gray0 gray1 gray2 igrym igry0 igry1
     igry2 gradm grad0 grad1 grad2 igrdm igrd0 igrd1 igrd2 grau0 baa_yld1000_avg1st, r bw(auto) small;
      
     // Output Growth
     lincom graym+gray0+gray1+gray2; 
     // Inflation
     lincom gradm+grad0+grad1+grad2; 
     // Output Growth Revision
     lincom igrym+igry0+igry1+igry2; 
     //Inflation Revision
     lincom igrdm+igrd0+igrd1+igrd2;
      
     ////////////////////////////////////////
     // Regression described in Equation (18)
     // Table (5)
     ///////////////////////////////////////

     // First Column (no corporate credit spreads)    
     ivreg2 ffr_cg avggrad12 ygap0 gray0 ffr_cg_l1 ffr_cg_l2, r bw(auto) small;

     //Interest Rate Smoothing
     lincom ffr_cg_l1+ffr_cg_l2;

     reg ffr_cg avggrad12 gray0 ygap0 ffr_cg_l1 ffr_cg_l2, b noheader notable;
     fitstat;

     //Store residuals
     predict double res_`res_count' , residuals;
     local ++res_count;

     // Second Column (with corporate credit spreads)    

     ivreg2 ffr_cg baa_yld1000_avg5d  avggrad12 gray0 ygap0 ffr_cg_l1 ffr_cg_l2, r bw(auto) small;

     //Interest Rate Smoothing
     lincom ffr_cg_l1+ffr_cg_l2;

     reg ffr_cg baa_yld1000_avg5d avggrad12 gray0 ygap0 ffr_cg_l1 ffr_cg_l2, b noheader notable;
     fitstat;

     predict double res_`res_count' , residuals;
     local ++res_count;

     local res_count = `res_count' - 1;
     forvalue i = 1/`res_count' {;
             replace res_`i' = 0 if res_`i' ==.;
     };
     tostring month year, replace;
     gen day = "01";
     gen export_date = month + "/" + day + "/" + year;
     destring month year, replace;
     sort day year  month ;
     keep res_* export_date;
     order export_date;
     export	delimited using "RR-CG-residuals", replace;
      
     drop res_*;
