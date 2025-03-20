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
global baseline_child  sex black hispanic anga bw mdi_12cor
global baseline_mother mage maternal_ppvt_r meduc works married tot_siblings employed_adult 
global outputs         mdi_24cor iqcage stndscor 
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

global inputs_std   baseline_humancapital_std cum_avg_daycare_36m_sum parenting_ages13_std
global outputs_std  early_humancapital_std
keep   ihdp site tg $baseline_child $baseline_mother $inputs_std $outputs $outputs_std

matrix desc = J(1,4,.)
foreach var of varlist $baseline_child $baseline_mother cum_avg_daycare_36m_sum parenting_ages13_std $outputs {
		reg      `var' tg, robust
		matrix  mc`var' = e(b)[1,2]
		matrix  md`var' = e(b)[1,1]
		matrix  mt`var' = md`var' + mc`var'
		matrix se`var'  = sqrt(e(V)[1,1])
		matrix   `var'  = [mt`var',mc`var',md`var',se`var']
		matrix     desc  = [desc \ `var']
}
matrix desc = desc[2...,1...]

// descriptives table
matrix all = [J(5,4,.)  \ desc[1..6,1...] \ J(1,4,.) \ desc[7..13,1...] \ J(1,4,.) \ desc[14..15,1...] \ J(1,4,.) \ desc[16...,1...]]
clear
svmat all
gen   n = _n

// format
foreach var of varlist all* {
	format  `var' %15.3fc
	replace `var' = round(`var',.001) 
	gen    `var'and = "&" if n !=1
}

// header
gen all0and = "&" if n != 1
gen       all1_0 = "\begin{tabular}{ L{7cm} C{2.25cm} C{2.25cm} C{2.25cm} C{2.25cm} } \toprule"  if n == 1

replace all1_0   = " "                    	      		    if n == 2
replace all0and  = " & \multicolumn{1}{c}{\textbf{(1)}}"    if n == 2
replace all1and  = " & \multicolumn{1}{c}{\textbf{(2)}}"    if n == 2
replace all2and  = " & \multicolumn{1}{c}{\textbf{(3)}}"    if n == 2
replace all3and  = " & \multicolumn{1}{c}{\textbf{(4)}}"    if n == 2

replace all1_0   = " "                    	      		    	    if n == 3
replace all0and  = " & \multicolumn{1}{c}{\textbf{Treatment}}"      if n == 3
replace all1and  = " & \multicolumn{1}{c}{\textbf{Control}}"        if n == 3
replace all2and  = " & \multicolumn{1}{c}{\textbf{Difference}}"     if n == 3
replace all3and  = " & \multicolumn{1}{c}{\textbf{Standard Error}}" if n == 3

replace all1_0   = " "                    	      		    if n == 4
replace all0and  = " & \multicolumn{1}{c}{(mean)}"          if n == 4
replace all1and  = " & \multicolumn{1}{c}{(mean)}"          if n == 4
replace all2and  = " & \multicolumn{1}{c}{ (1) $-$ (2)}"    if n == 4
replace all3and  = " & \multicolumn{1}{c}{(of difference)}" if n == 4

replace all0and  = "  \multicolumn{5}{l}{\underline{\textbf{Panel a. Child at Baseline}}}"     if n == 5
replace all1and  = "   "           						      		                		   if n == 5
replace all2and  = "   "           						      		                		   if n == 5
replace all3and  = "   "           						      		               			   if n == 5

replace all0and  = "  \multicolumn{5}{l}{\underline{\textbf{Panel b. Mother and Household at Baseline}}}"      if n == 12
replace all1and  = "   "           						      		                  						   if n == 12
replace all2and  = "   "           						      		                  				           if n == 12
replace all3and  = "   "           						      		                  				           if n == 12

replace all0and  = "  \multicolumn{5}{l}{\underline{\textbf{Panel c. Childcare and Parenting (Ages 1 to 3)}}}"  if n == 20
replace all1and  = "   "           						      		                   			  				if n == 20
replace all2and  = "   "           						      		                   			  				if n == 20
replace all3and  = "   "           						      		                   			 			    if n == 20

replace all0and  = "  \multicolumn{5}{l}{\underline{\textbf{Panel d. Measures of Early Life Skills (Ages 2 to 3)}}}"    if n == 23
replace all1and  = "   "           						      		                		   			  				if n == 23
replace all2and  = "   "           						      		                		   			  				if n == 23
replace all3and  = "   "           						      		                		   			  				if n == 23


replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Male}}}" 		 		                 if n == 6
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Black}}}"                              if n == 7
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Hispanic}}}"                    		 if n == 8
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Gestational Age (weeks)}}}"     		 if n == 9
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Birth Weight (kilograms)}}}"    		 if n == 10
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Bayley Mental Development Index}}}"    if n == 11
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Age}}}" 		 		          							 if n == 13
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{IQ-Peabody Picture Vocabulary Test}}}"                     if n == 14
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Education}}}"          		  		 if n == 15
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Works}}}"              	      		 if n == 16
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Married}}}"  					  		 if n == 17
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Other Children in Household}}}" 		 if n == 18
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Employed Adults in Household}}}"  	 if n == 19
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Childcare Hours per Week}}}"    		 if n == 21
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Parenting}}}"      			  		 if n == 22
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{Bayley Mental Development Index}}}"    if n == 24
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{IQ-Peabody Picture Vocabulary Test}}}" if n == 25
replace all1_0  = "   \multicolumn{1}{l}{\hspace{2mm} \textbf{\textit{IQ-Stanford-Binet Test}}}"     	  	 if n == 26


// order
global orderlist all1_0 all0and
foreach num of numlist 1(1)4 {
	global orderlist $orderlist all`num' all`num'and
}
order $orderlist
// other formating
replace all4and = "\\   "
replace all4and = "\\ \\" 							if n == 11  | n == 19 | n == 22
replace all4and = " \\ \midrule " 					if n == 2  
replace all4and = " \\ \bottomrule \end{tabular}"   if n == 26
replace all4and = " " 								if n == 1
replace all4and = "\\ \cmidrule(l{.15cm}r{.15cm}){2-2} \cmidrule(l{.15cm}r{.15cm}){3-3} \cmidrule(l{.15cm}r{.15cm}){4-4} \cmidrule(l{.15cm}r{.15cm}){5-5}" if n == 4

// put together in string
gen all = " "
foreach num of numlist 1(1)4 {
	tostring all`num' , replace force format(%15.3fc)
	replace  all`num'  = "" if all`num'  == "."
}

global orderin all1_0 all0and
foreach num of numlist 1(1)4 {
	global orderin $orderin all`num' all`num'and 
}
global orderin $orderin all`num'
order $orderin n

foreach var of varlist all1_0-all4and {	
	replace  all = all + `var'
	replace  all`num' = " " if all`num' == "."
}

// save in tex format
keep all
cd $output
outsheet using descriptive.tex, noquote nonames replace
