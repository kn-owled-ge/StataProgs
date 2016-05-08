**** Prepare GDP deflator ****
use "${datsrcfold}NGDP.dta", clear
gen fyear = year(date)
gen ngdp2005 = ngdp/13093.7        		   // base NGDP year = 2005 (13093 Billion)
destring rgdp, replace
destring def, replace
keep if rgdp<.
drop date

// find time structure of GDP series
sort fyear
gen t = _n
tsset t
gen lrgdp = log(rgdp)
gen lngdp = log(ngdp)
regress lrgdp t
predict rsdlr,resid
regress lngdp t
predict rsdln,resid
var rsdlr , lags(1)
tabstat rsdlr, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f) save
var rsdln , lags(1)
tabstat rsdln, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
disp "sigepsx: " sqrt(0.061^2 * (1-0.935^2))
keep fyear ngdp2005 ngdp rgdp def
save "${datfold}NGDP.dta", replace



**** Tidy data ****
use "${datsrcfold}COMPUCRSP_141111.dta", clear
ren variabl0 lt
drop xt prstkcc prstkpc scstkc

// make sic and gvkey into numeric
destring sic, replace
destring gvkey, replace

// generate market value
generate mve = csho*prcc_f
keep if mve>0
label variable mve "Market value of equity at FY end"
drop csho prcc_f mkvalt
drop datadate sich exchg lpermno fyr

// correct nonsensical and missing values
replace sppe=. if sppe<0
replace capx=. if capx<0
replace dp=. if dp<0
replace txditc=. if txditc<0
replace dvc=. if dvc<0
replace dvp=. if dvp<0
replace sstk=. if sstk<0
replace prstkc=. if prstkc<0
replace pstk = pstkrv if pstk == .
replace pstk = pstkl if pstk == .
replace pstk=. if pstk<0
replace dltis=. if dltis<0
replace dltr=. if dltr<0
replace xint=. if xint<0
replace xrd=. if xrd<0
replace xad=. if xad<0
replace xsga=. if xsga<0
replace ppegt=. if ppegt<0
replace sale=. if sale<0
replace cogs=. if cogs<0
replace revt=. if revt<0
replace aqc=abs(aqc) //aqc used as measure of structural firm changes
drop pstkrv pstkl

// lose repeat obs and declare panel
bysort gvkey fyear: keep if _n==1
xtset gvkey fyear

// Import the GDP data
merge m:1 fyear using "${datfold}NGDP.dta" , nogen keep(3)

// we deflate by total (nominal) GP, which will make data real and out of a fixed size economy
gen def1 = gp
// remove utilities / financial / software firms, as not included in our eventual sample
replace def1 = . if ((sic>=4900 & sic<=4999) | (sic>=6000 & sic<=6999) | (sic>=9000) | (sic>=7370 & sic<=7372))

bysort fyear: egen Tdef = total(def1)
drop def1
// make deflator base year = 2005 - that was a good year to be me :)
tabstat Tdef if fyear==2005, s(mean) save
mat tmp = r(StatTotal)
replace Tdef = Tdef/tmp[1,1]

// Correct all cash variables
foreach x of varlist * {
	if "`x'" != "gvkey" && "`x'" != "fyear" && "`x'" != "sic" && "`x'" != "Tdef" {
		replace `x' = `x'/ngdp2005
		//USE TDEF AS DEFLATOR: replace `x' = `x'/Tdef
	}
}

/*
// DEBUG - verify size of compustat universe now stable
keep if fyear>1981
drop if ((sic>=4900 & sic<=4999) | (sic>=6000 & sic<=6999) | (sic>=9000) | (sic>=7370 & sic<=7372))
bysort fyear: gen rep = 1 if _n==1
foreach x of varlist sale gp ppegt at mve {
	bysort fyear: egen T`x' = total(`x')
	replace T`x' = log(T`x')
	scatter T`x' fyear if rep==1,name(T`x')
	drop T`x'
}
*/
drop Tdef ngdp rgdp def
save "${datfold}Interim1.dta", replace




**** Interpolate missing values using two obs from each side ****
use "${datfold}Interim1.dta", clear
// Create obs with missing values for firms missing one year in the middle
// Pick the missing spells
sort gvkey fyear
gen tmp=1 if gvkey==gvkey[_n+1] & fyear==fyear[_n+1]-2
keep if tmp==1
drop tmp
replace fyear=fyear+1

// Null all cash variables
foreach x of varlist * {
	if "`x'" != "gvkey" && "`x'" != "fyear" && "`x'" != "sic" {
		replace `x' = .
	}
}

// Add dataset back
append using "${datfold}Interim1.dta"
sort gvkey fyear

// Who are legit candidates for interpolation?
gen int_me = 1
forvalues i = -2/2 {
	if `i'==0 continue
	replace int_me = 0 if (gvkey[_n + (`i')]!=gvkey)
	replace int_me = 0 if (fyear[_n + (`i')]!=fyear + (`i'))
}

// Do interpolation for all cash variables
foreach x of varlist * {
	if "`x'" != "gvkey" && "`x'" != "fyear" && "`x'" != "sic" {
		replace `x' = (`x'[_n-2] + `x'[_n-1] + `x'[_n+1] + `x'[_n+2])/4 if `x'==. & int_me==1
	}
}
drop int_me

// Create age variable
sort gvkey fyear
by gvkey: egen minyear = min(fyear)
gen age = fyear - minyear + 1
drop minyear

// Derive diff variables needed later
foreach x of varlist che invt at lt seq txditc pstk {
	gen D`x'    = `x'-L.`x'
}
save "${datfold}Interim2.dta", replace



foreach RD_NRD in "RD" "NRD" "BOTH" {
**** Set sample selection criteria ****
clear all
use "${datfold}Interim2.dta", clear

// Year to begin sample - due to lags, 1983 is actual beginning, with good RD data starting at 1982
local BEG_YEAR 			1982
// Year to end sample
local END_YEAR			2012
// Minimal PPEGT to be included in sample
local MIN_PPEGT			1
// Minimal MVE to be included in sample
local MIN_MVE			1
// Minimal V to be included in sample
local MIN_V			    1
// Minimal XRD/V ratio to be included in RD sample
local MIN_XRDV			0.01
// Maximal AQC/AT ratio to be included in sample
local MAX_AQC			0.1
// Maximal SPPE/PPEGT ratio to be included in sample
local MAX_SPP			0.1
// Which type of industry classification?
local FFIND 			48



**** Zero out missing and generate variables ****
count
// drop missing important member of composite variables
drop if ppegt>=.
drop if mve>=. | at>=. | seq>=.
drop if capx>=.
drop if dvc>=. | xint>=.
drop if gp>=.
drop if oibdp>=.
drop if xsga>=.
drop if dp>=.

// anyhting still missing is assumed zero
foreach x of varlist * {
	if "`x'" != "gvkey" && "`x'" != "fyear" && "`x'" != "sic" {
		replace `x' = 0 if `x'>=.
	}
}

gen K = ppegt
label variable K "Physical Capital"
gen SL = sale
label variable SL "Sales"
gen GP = gp
label variable GP "Gross Profit"
gen XS = xsga
label variable XS "XSGA"
gen OI = GP - XS    // == oibdp (99.9%)
label variable OI "Operating Income"
gen RD = xrd
label variable RD "R&D Investment"
gen OI1 = GP - XS + RD 
label variable OI1 "Operating Income plus RD"
gen IK  = capx-sppe
label variable IK "Physical Investment"
gen DP  = dp
label variable IK "Physical Depreciation"
gen D = (dvc + dvp + prstkc) + (xint - (Dat - (Dseq+Dtxditc-Dpstk))) + (Dinvt + Dche)
label variable D "Distributions"
gen V = mve + (at - (seq+txditc-pstk)) - (invt + che)
label variable V "Firm Value"


**** define subsample to use ****
gen touse = 1

// mark firms changing considerably
gen chg1 = aqc/at
gen chg2 = sppe/ppegt

// apply sample criteria
count if touse
// lose obs missing beg-of-year data
replace touse=0 if L.K>=. | L.V>=.		//lose 19K
count if touse
// not in year range
replace touse=0 if fyear < `BEG_YEAR' | fyear > `END_YEAR' | fyear==.  //lose 36K
count if touse
// utilities / financial firms
replace touse=0 if (sic>=4900 & sic<=4999) | (sic>=6000 & sic<=6999) | (sic>=9000) //lose 6K
count
// software firms
replace touse=0 if (sic>=7370 & sic<=7372) //lose 6K
count if touse
// firms changing considerably
replace touse=0 if chg1 > `MAX_AQC' | (L.chg1 > `MAX_AQC' & L.chg1 < .) | chg1==. 	//lose 12K
replace touse=0 if chg2 > `MAX_SPP' | (L.chg2 > `MAX_SPP' & L.chg2 < .) | chg2==. 	//lose 5K
drop chg1 chg2
count if touse
// too low PPEGT (K)
replace touse=0 if K < `MIN_PPEGT' | K==. | L.K < `MIN_PPEGT' | L.K==. | L2.K < `MIN_PPEGT' | L2.K==.		//lose 9K
count if touse
// too low MVE
replace touse=0 if mve < `MIN_MVE' | mve==. 	//lose 0K
count if touse
// too low V
replace touse=0 if V < `MIN_V' | V==. | L.V < `MIN_V' | L.V==. | L2.V < `MIN_V' | L2.V==.	//lose 1K
count if touse

**** generate FF industries ****
ffind(sic), newvar("ind") type(`FFIND')

// too little / too much R&D
gen RDV = RD/L.V

if "`RD_NRD'" == "RD" {
	replace touse=0 if RDV < `MIN_XRDV' 	//lose 53K
}
else if "`RD_NRD'" == "NRD" {
	replace touse=0 if RDV > `MIN_XRDV' 	//lose 33K
}
drop RDV
count if touse

save "${datfold}Interim3`RD_NRD'.dta", replace
}

