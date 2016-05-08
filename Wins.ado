capture program drop Wins
program define Wins, nclass
	version 10.1
	syntax varlist(numeric) [if] [in] [, Perc(real 1) SORTBY(varname)]
	
	marksample touse , novarlist
	
	quietly count if `touse'
	if `r(N)' == 0 error 2000
	
//	gen TMPPOSLOC = _n
	local lperc = 100-`perc'
	
	if "`sortby'"!="" {
		sort `sortby'
	}
	foreach v of varlist `varlist' {
		if "`sortby'"!="" {
			qui: by `sortby': egen `v'_TMP_pu = pctile(`v') if `touse', p(`lperc')
			qui: by `sortby': egen `v'_TMP_pl = pctile(`v') if `touse', p(`perc')
		}
		else {
			qui: egen `v'_TMP_pu = pctile(`v') if `touse', p(`lperc')
			qui: egen `v'_TMP_pl = pctile(`v') if `touse', p(`perc')
		}
		qui: replace `v'=`v'_TMP_pu if `v'>`v'_TMP_pu & `v'<. & `v'_TMP_pu<.
		qui: replace `v'=`v'_TMP_pl if `v'<`v'_TMP_pl & `v'<. & `v'_TMP_pl<.
		drop `v'_TMP_pu `v'_TMP_pl
	}
//	sort TMPPOSLOC
//	drop TMPPOSLOC
end
