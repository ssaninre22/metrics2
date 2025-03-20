* use the function below to import data directly from excel spreadsheet
* import excel "DataAP.xlsx", sheet("IndividChoiceData") firstrow
* or use .dta file with the same data
*
use "DataAP.dta"
*
* The formated output below is used to produce Table 2 and parts of Tables A1 and A2 in the paper
*
* Logit 4-outcome specification: EUT safe AC, EUT risky BD, fan-out AD and fan-in BC 
* Logit coefficients are reported in Table A2
*
mlogit y4 lpi lp1i i l s o phpl, b(0) noci
*
* marginal effects are reported in Table 2 and A2 (analogous parts for other models below)
*
margins, dydx(*)  predict(pr outcome(0)) noci
margins, dydx(*)  predict(pr outcome(1)) noci
margins, dydx(*)  predict(pr outcome(2)) noci
margins, dydx(*)  predict(pr outcome(3)) noci
*
* this part is used to obtain clustered standard errors part (analogous parts for other models below)
*
mlogit y4 lpi lp1i i l s o phpl, b(0) noci vce(cluster g)
margins, dydx(*)  predict(pr outcome(0)) noci
margins, dydx(*)  predict(pr outcome(1)) noci
margins, dydx(*)  predict(pr outcome(2)) noci
margins, dydx(*)  predict(pr outcome(3)) noci
*
* Logit 2-outcome specification: EUT vs non-EUT in Tables A1 and A2
*
mlogit y2 lpi lp1i i l s o phpl, b(0) noci
margins, dydx(*)  predict(pr outcome(1)) noci
mlogit y2 lpi lp1i i l s o phpl, b(0) noci vce(cluster g)
margins, dydx(*)  predict(pr outcome(1)) noci
*
* Linear probability model with 2 outcomes reported in Table A1
*
regress y2 lpi lp1i i l s o phpl, noci
regress y2 lpi lp1i i l s o phpl, noci vce(cluster g)
*
* Logit 2-outcome specification: EUT vs non-EUT excluding BC fan-in outcome reported in Tables A1 and A2
*
mlogit y2 lpi lp1i i l s o phpl if y4<3, b(0) noci
margins, dydx(*)  predict(pr outcome(1)) noci
mlogit y2 lpi lp1i i l s o phpl if y4<3, b(0) noci vce(cluster g)
margins, dydx(*)  predict(pr outcome(1)) noci
*
* Logit 3-outcome specification: EUT, fan-out AD and fan-in BC outcome reported in Tables A1 and A2
*
mlogit y3 lpi lp1i i l s o phpl, b(0) noci
margins, dydx(*)  predict(pr outcome(0)) noci
margins, dydx(*)  predict(pr outcome(1)) noci
margins, dydx(*)  predict(pr outcome(2)) noci
mlogit y3 lpi lp1i i l s o phpl, b(0) noci vce(cluster g)
margins, dydx(*)  predict(pr outcome(0)) noci
margins, dydx(*)  predict(pr outcome(1)) noci
margins, dydx(*)  predict(pr outcome(2)) noci
*
* Ordered logit 3-outcome specification in the following order: fan-in BC, EUT, fan-out AD reported in Tables A1 and A2
*
ologit y3o lpi lp1i i l s o phpl, noci
margins, dydx(*)  predict(pr outcome(0)) noci
margins, dydx(*)  predict(pr outcome(1)) noci
margins, dydx(*)  predict(pr outcome(2)) noci
ologit y3o lpi lp1i i l s o phpl, noci vce(cluster g)
margins, dydx(*)  predict(pr outcome(0)) noci
margins, dydx(*)  predict(pr outcome(1)) noci
margins, dydx(*)  predict(pr outcome(2)) noci