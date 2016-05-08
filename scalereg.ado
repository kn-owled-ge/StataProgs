*! version 1.0 scalereg estimate a scale regression

/* 

Author: Robert Parham, University of Rochester

http://kn.owled.ge

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

*/


capture program drop scalereg
program define scalereg, eclass
	version 12
	syntax varlist(min=4 max=4 numeric) [if] [in] , [MINscale(real -100.0) MAXscale(real 100.0) Step(integer 20) METhod(string) REGMEThod(string) MINBin(integer 10) REPeat(integer 1)]
	marksample touse , novarlist
	
	qui{
	tempname N
	quietly count if `touse'
	local `N' = `r(N)'
	if ``N'' == 0 error 2000
	
	// check which method to use
	tempname met
	local method = strupper("`method'")
	if ("`method'"=="" | "`method'"=="MEDIAN") 	local `met' = 0
	else if ("`method'"=="MEAN") 				local `met' = 1
	else if ("`method'"=="MAD") 				local `met' = 2
	else if ("`method'"=="IQR") 				local `met' = 3
	else {
		di as err "Unknown value provided for Method. Choose one of MEDIAN,MEAN,MAD,IQR."
		exit 197
	}

	// check which regression method to use
	tempname regmet
	local regmethod = strupper("`regmethod'")
	if ("`regmethod'"=="" | "`regmethod'"=="OLS") 	local `regmet' = 0
	else if ("`regmethod'"=="LAD") 					local `regmet' = 1
	else {
		di as err "Unknown value provided for Method. Choose one of OLS,LAD."
		exit 197
	}

	
	// separate varlist
	tempname X bin
	gettoken `X' varlist: varlist
	gettoken `bin' varlist: varlist
	disp "``X''"
	disp "``bin''"
	disp "`minscale'"
	disp "`maxscale'"
	disp "`step'"
	disp "``met''"
	disp "`minbin'"

	// verify X non-negative
	tempname tmpres minX
	tabstat ``X'' , s(min) save
	mat `tmpres' = r(StatTotal)
	local `minX' = `tmpres'[1,1]
	if (``minX''<=0) {
		di as err "Scale variable ``X'' must take only positive values. Aborting."
		exit 197
	}
	
	// generate a log value for X
	tempvar logX
	gen `logX' = log(``X'')
	
	// generate bins and weights and elect a representative bin member
	tempvar wgt binrep
	replace ``bin'' = round(`logX'*`step')/`step'
	bysort ``bin'': egen `wgt' = count(``bin'')
	bysort ``bin'': gen `binrep' = 1 if _n==1 & _N>`minbin'

	// get the Y variable
	tempname Y mebY slope 
	tempvar bY logmebY iter2
	gettoken `Y' varlist: varlist
	gettoken `mebY' varlist: varlist
	disp "``Y''"
	disp "``mebY''"

	// run the binning procedure `repeat' times
	local `slope' = 1
	gen `logmebY' = .
	// process very stable after 2 times
	forvalues `iter2' = 1/`repeat' {
		gen `bY' = ``Y''/exp((`logX'-``bin'')*``slope'')
		drop ``mebY''
		if (``met'' == 0) {
			bysort ``bin'': egen ``mebY''  = median(`bY')
		}
		else if (``met'' == 1) {
			bysort ``bin'': egen ``mebY''  = mean(`bY')
		}
		else if (``met'' == 2) {
			bysort ``bin'': egen ``mebY''  = mad(`bY')
		}
		else {
			bysort ``bin'': egen ``mebY''  = iqr(`bY')
		}
		replace `logmebY' = log(``mebY'')
		drop `bY'
		
		// run the regression on bin centers, with frequency weights
		if (``regmet'' == 0) {
			regress `logmebY' ``bin'' [fweight=`wgt'] if `binrep'==1 & ``bin''>`minscale' & ``bin''<`maxscale'
		}
		else {
			qreg `logmebY' ``bin'' [fweight=`wgt'] if `binrep'==1 & ``bin''>`minscale' & ``bin''<`maxscale'
		}

		// update slope for next iteration
		mat `tmpres' = e(b)
		local `slope' = `tmpres'[1,1]
	}
	drop `logmebY' `logX' `binrep' `wgt'
	
	// save return matrix
	mat `tmpres' = e(b)
	matname `tmpres' "``X'' _cons"  , c(.)
	matname `tmpres' "``Y''"  , r(.)
	ereturn post `tmpres' , depname(``Y'') obs(``N'') esample(`touse')
	
	
	}
end
