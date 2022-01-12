/******************************************************************************/
/* Project: Bartik Instruments - An Applied Introduction					  */
/* Author:	M. Breuer														  */
/* Date:	01/27/2021														  */
/******************************************************************************/
/* Description:																  */
/* Illustration of Bartik instrument construction using simulated data		  */
/* Panel A: unorrelated factor (b_l_c=0, b_e_c=0, b_y_c=1, freq_c=0.1)		  */
/* Panel B: positively correlated factor (b_e_c=0.1) 						  */
/* Panel C:	negatively correlated factor (b_e_c=-0.1)						  */
/******************************************************************************/

/* Version */
version 15.1

/* Necessary programs */
ssc install estout, replace
ssc install ftools, replace
ssc install reghdfe, replace
ssc install ivreghdfe, replace
ssc install ranktest, replace

/* Directory */
local directory = "..." // Insert working directory
cd "`directory'"

/* Delete prior version */
capture rm Influence.dta
capture rm Regressions.smcl

/* Loop: impact of confounder on employment share */
forvalues s = -0.25(0.25)0.25 {

/* Loop: impact of confounder on employment rate */
forvalues r = -0.25(0.05)0.25 {

/* Preliminaries */
clear all
set more off
set seed 1234


/******************************************************************************/
/* Step 1: Generate data													  */
/******************************************************************************/

/* Parameters */

	/* Counties */
	local counties = 1000
	
	/* Industries */
	local industries = 30
	
	/* Years */
	local years = 15

	/* Confounder */
	
		/* Impact on employment share [units: log employment] */
		local b_l_c = round(`s', 0.01)
		
		/* Impact on employment rate [units: employment rate] */
		local b_e_c = round(`r', 0.01)
		
		/* Impact on wage [units: log wage] */
		local b_y_c = 1
		
		/* Frequency */
		local freq_c = 0.1
		
	/* Employment */
	
		/* Employment share */
				
			/* Mean [log employment] */
			local m_l = 5
			
			/* Standard deviation [fraction of mean] */
			local sd_l = 0.1
					
			/* Persistence */
			local rho_l = 0.7

		/* Employment rate */
					
			/* Mean [employment rate] */
			local m_e = 0.85
			
			/* Standard deviation [fraction of mean] */
			local sd_e = 0.05
			
			/* Persistence */
			local rho_e = 0.7

			/* Commonality (across counties) */
			local com_e = 0.8	
			
			/* Lower bound */
			local low_e = 0
			
			/* Upper bound */
			local up_e = 1
			
	/* Wage */
	
		/* Mean [log average monthly wage] */
		local m_y = 1
		
		/* Residual variation [fraction of mean] */
		local sd_y = 0.1
	
		/* Employment impact */
		local b_e_y = 1
		
/* Panel */

	/* Observations */
	local observations = `counties'*`industries'*`years'
	set obs `observations'

	/* County ID */
	gen county = ceil(_n/(`industries'*`years'))
	label var county "County"

	/* Industry ID */
	gen industry = ceil(_n/(`years')) - (county - 1)*`industries'
	label var industry "Industry"
	
	/* County-Industry ID */
	egen ci = group(county industry)
	label var ci "County-Industry"
	
	/* Year ID */
	gen year = _n - (ci - 1)*`years'
	label var year "Year"
	
	/* Set panel */
	xtset ci year
	
	
/* Observed Wage and Unemployment Data */
	
	/* County-industry level */
		
		/* Confounder */
		
			/* Technology shock */
			gen confounder = (runiform() < `freq_c')
			label var confounder "Confounder"
			
		/* Employment share */
		
			/* AR(1) model of employment */
				
				/* Intercept */
				gen a_ik = rnormal(`m_l', `sd_l'*`m_l') if year == 1
				by ci: replace a_ik = a_ik[1]
				
				/* First year */
				gen l_itk = a_ik + rnormal(0, `sd_l'*a_ik) if year == 1
				
				/* Over time */
				replace l_itk = (1-`rho_l')*a_ik + `rho_l'*l.l_itk + rnormal(0, `sd_l'*a_ik) if year > 1

				/* Confounding influence */
				replace l_itk = l_itk + `b_l_c'*confounder
				
				/* Convert from log employment to employment */
				replace l_itk = exp(l_itk)
				label var l_itk "Employment"

				
			/* Share */
			
				/* Total */
				egen l_it = total(l_itk), by(county year) missing

				/* Employment share */
				gen w_itk = l_itk/l_it
				label var w_itk "Employment share"
			
				/* Drop */
				drop a_* l_it

		
		/* Employment rate */
						
			/* AR(1) model with commonality */
		
				/* Intercept */
				gen a_ik = rnormal(`m_e', `sd_e'*`m_e') if year == 1
				by ci: replace a_ik = a_ik[1]				
		
				/* Common trend */
		
					/* First year (and county) */
					gen c_tk = rnormal(0, `sd_e'*a_ik) if year == 1 & county == 1
					
					/* Over time */
					replace c_tk = `rho_e'*l.c_tk + rnormal(0, `sd_e'*a_ik) if year > 1
						
					/* Across counties */
					egen mean = mean(c_tk), by(industry year)
					replace c_tk = mean
					drop mean					
				
				/* Idiosyncratic trend (with confounder) */
				gen e_itk = rnormal(0, `sd_e'*a_ik)

				/* Combined */
				replace e_itk = a_ik + (1-`com_e')*e_itk + `com_e'*c_tk
					
				/* Confounding influence */
				replace e_itk = e_itk + `b_e_c'*confounder
					
				/* Respect bounds */
				replace e_itk = `low_e' if e_itk < `low_e'
				replace e_itk = `up_e' if e_itk > `up_e'
				label var e_itk "Employment rate"
				
				/* Drop */
				drop a_* c_tk				

		/* Wage */
		
			/* Intercept */
			gen a_ik = rnormal(`m_y', `sd_y'*`m_y') if year == 1
			by ci: replace a_ik = a_ik[1]		
		
			/* Wage */
			gen y_itk = a_ik + `b_e_y'*e_itk + `b_y_c'*confounder + rnormal(0, `sd_y'*a_ik)
			label var y_itk "Wage (log)"
			
			/* Drop */
			drop a_* confounder	
			
		
	/* County level */
	
		/* Employment rate */
		egen x_it = total(w_itk*e_itk), by(county year) missing
		label var x_it "Employment rate (county)"
		
		/* Wage */
		egen y_it = total(w_itk*y_itk), by(county year) missing 
		label var y_it "Wage (log; county)"
	
	
/******************************************************************************/
/* Step 2: Construct Bartik instrument										  */
/******************************************************************************/

/* County-industry level */
		
	/* Pre-determined employment share */
	sort county industry year
	by county industry: gen w_ik = w_itk[1]
	label var w_ik "Employment share (county; industry)"
	
	/* Nationwide employment rates */
	
		/* All county-industries */
		egen l_tk = total(l_itk), by(industry year) missing
		egen e_tk = total(l_itk/l_tk*e_itk), by(industry year) missing
	
		/* Adjusted for own county-industry [finite sample adjustment] */
		replace e_tk = l_tk/(l_tk - l_itk)*(e_tk - (l_itk/l_tk)*e_itk)
		label var e_tk "Employment rate (national; industry)"
		
		/* Drop */
		drop l_tk
		
/* County level */

	/* Bartik instrument */
	egen z_it = total(w_ik*e_tk), by(county year)
	label var z_it "Bartik instrument"
	
	
/******************************************************************************/
/* Step 3: Run regressions													  */
/******************************************************************************/

/* County level */

	/* Preserve */
	preserve

		/* Duplicates */
		duplicates drop county year, force

		/* OLS */
		qui reghdfe y_it x_it, a(county year) cluster(county)
		gen b_ols = _b[x_it]
		gen se_ols = _se[x_it]
		est store OLS
		
		/* Reduced Form */
		qui reghdfe y_it z_it, a(county year) cluster(county)
		est store RF	
		
		/* First Stage */
		qui reghdfe x_it z_it, a(county year) cluster(county)
		est store FS		
		
		/* IV */
		qui ivreghdfe y_it (x_it = z_it), a(county year) cluster(county)
		gen b_iv = _b[x_it]
		gen se_iv = _se[x_it]		
		est store IV
		
		/* Create output */ 
		
			/* Keep relevant variables */
			keep b_* se_*
			
			/* Keep relevant observation */
			keep if _n == 1
			
			/* Confounding influence */
			gen s = round(`s', 0.01)
			gen r = round(`r', 0.01)
			
			/* Append */
			capture append using Influence
			
			/* Save */
			save Influence, replace
			
		/* Condition for Table 2 */
		if round(`s',0.01) == 0 & (round(`r',0.01) == -0.1 | round(`r',0.01) == 0 | round(`r',0.01) == 0.1) {
				
			/* Log: open */
			log using Regressions, append smcl 
			
			/* Output */
			estout OLS RF FS IV, ///
				drop(_cons) ///
				cells(b(star fmt(3)) t(par fmt(2))) starlevels(* 0.10 ** 0.05 *** 0.01) ///
				title("Table 2: Bartik Instrument Estimation") ///								
				mlabels("OLS: y_it" "RF: y_it" "FS: x_it" "SS: y_it")  modelwidth(15) unstack ///
				stats(N N_clust r2, fmt(0 0 3) label("Obs." "Clusters" "R-Squared"))
			
			/* Log close */
			log close
			
		}
		
	/* Restore */
	restore 

	
/******************************************************************************/
/* Step 4: Create illustrative example 										  */
/******************************************************************************/

/* Condition for Table 1 */
if round(`s', 0.01) == 0 & round(`r', 0.01) == 0 {
	
	/* Keep relevant observations */
	keep if (county == 150 | county == 151) & (industry == 10 | industry == 11) & (year == 1 | year == 2)  

	/* Keep relevant items */
	keep county  year industry w_itk e_itk x_it w_ik e_tk z_it

	/* Order */
	order county year industry w_itk e_itk x_it w_ik e_tk z_it
		
	/* Sort */
	sort county year industry

	/* Round */
	foreach var of varlist w_itk e_itk x_it w_ik e_tk z_it {
	
		/* Replace */
		replace `var' = round(`var', 0.001)
	
	}
	
	/* Save */
	save Example, replace
	
}
}
}

/******************************************************************************/
/* Step 5: Create bias graphs												  */
/******************************************************************************/

/* Estimates data */
use Influence, clear

/* Confidence intervals */
		
	/* OLS */
	gen ci_high_ols = b_ols + 1.96*se_ols
	gen ci_low_ols = b_ols - 1.96*se_ols
		
	/* IV */
	gen ci_high_iv = b_iv + 1.96*se_iv
	gen ci_low_iv = b_iv - 1.96*se_iv
	
/* Graphs */
		
	/* OLS vs. IV */
	graph twoway ///
		(rarea ci_high_ols ci_low_ols r if round(s, 0.01) == 0, color(gs10)) ///
		(rarea ci_high_iv ci_low_iv r if round(s, 0.01) == 0, color(gs10)) ///
		(line b_ols r if round(s, 0.01) == 0, lcolor(black)) ///
		(line b_iv r if round(s, 0.01) == 0, lcolor(black) lpattern(dash)) ///
			, legend(label(1 "CI 95%") label(3 "{&beta}{superscript:OLS}") label(4 "{&beta}{superscript:IV}") rows(3) order(3 4 1) ring(0) position(11) bmargin(medium) symxsize(5)) /// 
			graphregion(color(white)) plotregion(fcolor(white)) ///
			xtitle("Impact of Confounder on {it:e{subscript:i,t,k}}") ///
			xlabel(-0.25(0.1)0.25, format(%9.2f)) ///
			xline(0, lcolor(black) lwidth(thin)) ///			
			yline(0, lcolor(black) lwidth(thin)) ///
			ylabel(-4(1)6, format(%9.0f) angle(0)) /// 
			ytitle("{&beta}", orientation(horizontal)) ///			
			name(Figure_1, replace) saving(Figure_1, replace)	
			graph export Figure_1.png, replace

	/* Moderating Impact on Share */
	graph twoway ///
		(rarea ci_high_ols ci_low_ols r if round(s, 0.01) == 0.25, color(gs10)) ///
		(rarea ci_high_ols ci_low_ols r if round(s, 0.01) == -0.25, color(gs10)) ///
		(line b_ols r if round(s, 0.01) == 0.25, lcolor(black)) ///
		(line b_ols r if round(s, 0.01) == -0.25, lcolor(black) lpattern(dash)) ///
			, legend(label(1 "CI 95%") label(3 "Positive Impact on {it:w{subscript:i,t,k}}") label(4 "Negative Impact on {it:w{subscript:i,t,k}}") rows(3) order(3 4 1) ring(0) position(11) bmargin(medium) symxsize(5)) /// 
			graphregion(color(white)) plotregion(fcolor(white)) ///
			xtitle("Impact of Confounder on {it:e{subscript:i,t,k}}") ///
			xlabel(-0.25(0.1)0.25, format(%9.2f)) ///
			xline(0, lcolor(black) lwidth(thin)) ///			
			yline(0, lcolor(black) lwidth(thin)) ///
			ylabel(-4(1)6, format(%9.0f) angle(0)) /// 
			ytitle("{&beta}{superscript:OLS}", orientation(horizontal)) ///
			name(Figure_2, replace) saving(Figure_2, replace)
			graph export Figure_2.png, replace
