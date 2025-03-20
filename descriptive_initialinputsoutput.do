version 18.0
set more off
clear all
set matsize 11000
set maxvar  120000
set seed 0
global bootstraps 1000

// set environment variables
global projects: env projects
global storageb : env storageb

// set general locations
global dataready = "$storageb/dc_data/"
// output
global output    = "$projects/dc/output"

cd  $dataready
use ihdp_data, clear
rename _all, lower
keep if twin == 0

replace bw = bw/1000
// list relevant variables
global baseline_child  sex black hispanic  anga bw mdi_12cor
global baseline_mother mage maternal_ppvt_r meduc works married tot_siblings employed_adult 
global outputs         iqcage stndscor mdi_24cor
global inputs          cum_avg_daycare_36m_sum

// mark sample baseline sample
reg $baseline_child $baseline_mother $outputs
keep if e(sample) == 1

// parenting latents, ages 1 and 3
foreach var of varlist alt_subscale_1_12-alt_subscale_6_12 alt_subscale_1_36-alt_subscale_8_36 {
	summ `var'
	gen  `var'_std = (`var' - r(mean))/r(sd)
}
# delimit
sem (_cons@0 X -> alt_subscale_1_12_std) (_cons@0 X -> alt_subscale_2_12_std) (_cons@0 X -> alt_subscale_3_12_std) 
	(_cons@0 X -> alt_subscale_4_12_std) (_cons@0 X -> alt_subscale_5_12_std) (_cons@0 X -> alt_subscale_6_12_std), method(adf);
predict parenting_age1, latent;
sem (_cons@0 X -> alt_subscale_1_36_std) (_cons@0 X -> alt_subscale_2_36_std) (_cons@0 X -> alt_subscale_3_36_std) (_cons@0 X -> alt_subscale_4_36_std)
	(_cons@0 X -> alt_subscale_5_36_std) (_cons@0 X -> alt_subscale_6_36_std) (_cons@0 X -> alt_subscale_7_36_std) (_cons@0 X -> alt_subscale_8_36_std), cov(e.alt_subscale_1_36_std*e.alt_subscale_2_36_std)
																																						 cov(e.alt_subscale_3_36_std*e.alt_subscale_4_36_std)
																																						 cov(e.alt_subscale_5_36_std*e.alt_subscale_6_36_std)
																																						 cov(e.alt_subscale_7_36_std*e.alt_subscale_8_36_std)
																																						 cov(e.alt_subscale_4_36_std*e.alt_subscale_7_36_std)
																																						 cov(e.alt_subscale_2_36_std*e.alt_subscale_7_36_std) method(adf);		predict parenting_age3, latent;
# delimit cr
egen    parenting_ages13 = rowmean(parenting_age1 parenting_age3)
keep if parenting_ages13 !=.

// initial human capital latent
foreach var of varlist $baseline_child maternal_ppvt_r {
	summ    `var'
	gen     `var'_std = (`var' - r(mean))/r(sd)
}
sem  (_cons@0 X->mdi_12cor_std) (_cons@0 X->bw_std) (_cons@0 X->anga_std) (_cons@0 X->maternal_ppvt_r_std), cov(e.bw_std*e.anga_std) method(adf)
predict baseline_humancapital, latent
drop *_std

// ages 2-3 latent
foreach var of varlist iqcage stndscor mdi_24cor {
	summ    `var' if tg == 0
	gen     `var'_std = (`var' - r(mean))/r(sd)
}
sem  (_cons@0 X-> iqcage_std) (_cons@0 X->stndscor_std) (_cons@0 X->mdi_24cor_std), method(adf)
predict early_humancapital, latent

// standardize latents in sample
foreach var of varlist baseline_humancapital parenting_ages13 early_humancapital {
	summ    `var' if tg == 0
	gen     `var'_std = (`var' - r(mean))/r(sd)
	replace `var'_std = `var'_std + 100
}

summ    cum_avg_daycare_36m_sum      if tg == 0 
gen     cum_avg_daycare_36m_sum_std  = (cum_avg_daycare_36m_sum        - r(mean))/r(sd)
replace cum_avg_daycare_36m_sum_std  = cum_avg_daycare_36m_sum_std     + r(mean)

global inputs_std   baseline_humancapital_std cum_avg_daycare_36m_sum cum_avg_daycare_36m_sum_std parenting_ages13_std
global outputs_std  early_humancapital_std
keep   ihdp site tg $baseline_child $baseline_mother $inputs_std $outputs $outputs_std

// distribution of inputs and outputs
foreach var of varlist $inputs_std $outputs_std {
	reg      `var' tg, robust 
	matrix  mc`var' = e(b)[1,2]
	matrix  md`var' = e(b)[1,1] 
	matrix se`var' = sqrt(e(V)[1,1])
}

// labels
foreach var of varlist $inputs_std $outputs_std {
	foreach stat in mc md se {
		local `stat'`var' = `stat'`var'[1,1]
		local `stat'`var' = string(``stat'`var'',"%9.3f")
	}
}

// baseline human capital
# delimit
global box  text(.375 97.5
	 "{bf: Control Mean}"
	 "{bf: `mcbaseline_humancapital_std'}"
	 " "
	 "{bf: Treatment-Control Mean Difference}" 
	 "{bf: `mdbaseline_humancapital_std' (s.e. = `sebaseline_humancapital_std')}"
, place(c) box size(small) just(c) margin(l+1 b+1 t+1 r+1) fcolor(none)); 
# delimit cr

#delimit
twoway  (kdensity baseline_humancapital_std if tg == 0, lwidth(dash) lcolor(gs0) lpattern(dash) lwidth(thick) yline(.6, lcolor(gs12)))
	    (kdensity baseline_humancapital_std if tg == 1, lcolor(gs6)  lwidth(thick)),	  
		  

		legend(label(1 "Control") label(2 "Treatment") position(6) size(medsmall) order(1 2) rows(1) region(lcolor(gs0)))
		  xlabel(96[2]104, labsize(medium)        grid   glcolor(gs12) glpattern(solid)) 
		  ylabel(0[.15].45,   labsize(medium) angle(h) glcolor(gs12) glpattern(solid)) 
		  ytitle("Density", size(medium))
		  xtitle("Factor Score")
		  graphregion(color(white)) plotregion(color(white)) $box;
#delimit cr
cd $output
graph export baseline_humancapital_std.eps, replace

// parenting
# delimit
global box  text(.375 97.5
	 "{bf: Control Mean}"
	 "{bf: `mcparenting_ages13_std'}"
	 " "
	 "{bf: Treatment-Control Mean Difference}" 
	 "{bf: `mdparenting_ages13_std' (s.e. = `separenting_ages13_std')}"
, place(c) box size(small) just(c) margin(l+1 b+1 t+1 r+1) fcolor(none)); 
# delimit cr

#delimit
twoway  (kdensity parenting_ages13_std if tg == 0, lwidth(dash) lcolor(gs0) lpattern(dash) lwidth(thick) yline(.6, lcolor(gs12)))
	    (kdensity parenting_ages13_std if tg == 1, lcolor(gs6)  lwidth(thick)),	  
		  

		legend(label(1 "Control") label(2 "Treatment") position(6) size(medsmall) order(1 2) rows(1) region(lcolor(gs0)))
		  xlabel(96[2]104, labsize(medium)        grid   glcolor(gs12) glpattern(solid)) 
		  ylabel(0[.15].45,   labsize(medium) angle(h) glcolor(gs12) glpattern(solid)) 
		  ytitle("Density", size(medium))
		  xtitle("Factor Score")
		  graphregion(color(white)) plotregion(color(white)) $box;
#delimit cr
cd $output
graph export parenting_ages13_std.eps, replace

// early human capital
# delimit
global box  text(.375 97.5
	 "{bf: Control Mean}"
	 "{bf: `mcearly_humancapital_std'}"
	 " "
	 "{bf: Treatment-Control Mean Difference}" 
	 "{bf: `mdearly_humancapital_std' (s.e. = `seearly_humancapital_std')}"
, place(c) box size(small) just(c) margin(l+1 b+1 t+1 r+1) fcolor(none)); 
# delimit cr

#delimit
twoway  (kdensity early_humancapital_std if tg == 0, lwidth(dash) lcolor(gs0) lpattern(dash) lwidth(thick) xline(96, lpattern(solid) lcolor(gs12)))
	    (kdensity early_humancapital_std if tg == 1, lcolor(gs6)  lwidth(thick)),	  
		  

		legend(label(1 "Control") label(2 "Treatment") position(6) size(medsmall) order(1 2) rows(1) region(lcolor(gs0)))
		  xlabel(96[2]104, labsize(medium)        grid   glcolor(gs12) glpattern(solid)) 
		  ylabel(0[.15].45,   labsize(medium) angle(h) glcolor(gs12) glpattern(solid)) 
		  ytitle("Density", size(medium))
		  xtitle("Factor Score")
		  graphregion(color(white)) plotregion(color(white)) $box;
#delimit cr
cd $output
graph export early_humancapital_std.eps, replace

// childcare
# delimit
global box  text(.65 65
	 "{bf: Control Mean}"
	 "{bf: `mccum_avg_daycare_36m_sum'}"
	 " "
	 "{bf: Treatment-Control Mean Difference}" 
	 "{bf: `mdcum_avg_daycare_36m_sum' (s.e. = `secum_avg_daycare_36m_sum')}"
, place(c) box size(small) just(c) margin(l+1 b+1 t+1 r+1) fcolor(none)); 
# delimit cr

#delimit
twoway  (histogram cum_avg_daycare_36m_sum if tg == 1,  width(10) discrete fraction bcolor(gs6)  							      gap(5))
		(histogram cum_avg_daycare_36m_sum if tg == 0,  width(10) discrete fraction bfcolor(none) blcolor(gs0)  blwidth(medthick) gap(5)),	  
		  

		legend(label(2 "Control") label(1 "Treatment") position(6) size(medsmall) order(2 1) rows(1) region(lcolor(gs0)))
		  xlabel(0   "0-9" 10 "10-19" 20 "20-29" 30 "30-39" 40 "40-49" 
				50 "50-59" 60 "60-69" 70 "70-79" 80 "80-89", labsize(medium)        grid   glcolor(gs12) glpattern(solid)) 
		  ylabel(0[.26].78,   labsize(medium) angle(h) glcolor(gs12) glpattern(solid)) 
		  ytitle("Fraction", size(medium))
		  xtitle("Average Hours of Childcare per Week")
		  graphregion(color(white)) plotregion(color(white)) $box;
#delimit cr
cd $output
graph export cum_avg_daycare_36m_sum.eps, replace

# delimit
global box  text(.7 11.5
	 "{bf: Control Mean}"
	 "{bf: `mccum_avg_daycare_36m_sum_std'}"
	 " "
	 "{bf: Treatment-Control Mean Difference}" 
	 "{bf: `mdcum_avg_daycare_36m_sum_std' (s.e. = `secum_avg_daycare_36m_sum_std')}"
, place(c) box size(small) just(c) margin(l+1 b+1 t+1 r+1) fcolor(none)); 
# delimit cr

#delimit
twoway  (histogram cum_avg_daycare_36m_sum_std if tg == 1,  width(1)  fraction bcolor(gs6)  							      gap(5))
		(histogram cum_avg_daycare_36m_sum_std if tg == 0,  width(1)  fraction bfcolor(none) blcolor(gs0)  blwidth(medthick)  gap(5)),	  
		  

		legend(label(2 "Control") label(1 "Treatment") position(6) size(medsmall) order(2 1) rows(1) region(lcolor(gs0)))
		  xlabel(5.5[1]12.5, labsize(medium)        grid   glcolor(gs12) glpattern(solid)) 
		  ylabel(0[.28].85,   labsize(medium) angle(h) glcolor(gs12) glpattern(solid)) 
		  ytitle("Fraction", size(medium))
		  xtitle("Standardized Average Hours of Childcare per Week")
		  graphregion(color(white)) plotregion(color(white)) $box;
#delimit cr
cd $output
graph export cum_avg_daycare_36m_sum_std.eps, replace
