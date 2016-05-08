**** General initial setup ****
clear all
set matsize 11000

global plc = 2
if ${plc} == 0 {
	global workfold = "C:\Users\finguy\SkyDrive\Documents\KnowledgeDynamics\\"
}
if ${plc} == 1 {
	global workfold = "D:\SkyDrive\Documents\KnowledgeDynamics\\"
}
if ${plc} == 2 {
	global workfold = "C:\Users\robert.parham\OneDrive\Documents\KnowledgeDynamics\\"
}
global scriptfold = "Prog\StataProgs\"
global resfold = "Results\ToPaper\"
global datsrcfold = "DataSrc\"
global datfold = "Data\"

cd ${workfold}
adopath + "${workfold}${scriptfold}"


**** Generate sample set ****
qui: do "${workfold}${scriptfold}DataPrep.do"


**** Get data into matlab ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count

foreach vr of varlist K V {
	gen uF_`vr' = `vr'
	gen u`vr'   = L.`vr'
	gen uL_`vr' = L2.`vr'
}

foreach vr of varlist GP XS OI IK RD D {
	gen uF_`vr' = F.`vr'
	gen u`vr'   = `vr'
	gen uL_`vr' = L.`vr'
}

local vrs = "uF_K uK uL_K uF_V uV uL_V uF_GP uGP uL_GP uF_XS uXS uL_XS uF_OI uOI uL_OI uF_IK uIK uL_IK uF_RD uRD uL_RD uF_D uD uL_D"
//keep touse gvkey fyear ind age `vrs'
gen X = 0

foreach vr of varlist `vrs' {
	if ("`vr'"=="uF_RD" | "`vr'"=="uRD"  | "`vr'"=="uL_RD") & "`RD_NRD'" ~= "RD" {
		continue
	}
	replace touse=0 if `vr'>=.
	if ("`vr'"=="uF_GP" | "`vr'"=="uGP"  | "`vr'"=="uL_GP") {
		replace touse=0 if `vr'<0.1
	}
	if ("`vr'"=="uF_XS" | "`vr'"=="uXS"  | "`vr'"=="uL_XS") {
		replace touse=0 if `vr'<0.1
	}
	if ("`vr'"=="uF_IK" | "`vr'"=="uIK"  | "`vr'"=="uL_IK") {
		replace touse=0 if `vr'<0.1
	}
}
count if touse

// only for Konly models, to avoid negative OI
//gen OIK = uOI/uK
//replace touse=0 if OIK<0.01

// gvkey fyear ind age Z L_Z X L_X W L_W F_N N L_N vrs
outsheet gvkey fyear ind age X X X X X X X X X `vrs' using data\data_160508`RD_NRD'.csv if touse, comma replace nonames nolabel




































**** Do initial sample statistics ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
cls

// sample base stats
count
tabstat K SL GP XS OI RD IK IX DP D V, s(mean sd median skewness) f(%8.3f)

// negative obs percentage count
qui{
count
mat tot = r(N)
local cnt = 1
matrix results = J(1,10,0)
foreach X of varlist K SL GP XS OI RD IK IX DP D V {
	count if `X'<0
	mat tmp = r(N)
	mat results[1,`cnt'] = (tmp[1,1]/tot[1,1])*100
	local cnt = `cnt' + 1
}
}
mat list results

// firm and year counts
bysort fyear: egen gvcount=count(gvkey)
tabstat gvcount , by(fyear)
bysort gvkey: gen rep=1 if _n==1
egen totgv = total(rep)
list totgv in 1/1

// view scale data
hist logK ,bin(100) xmtick(0(1)15) ymtick(0(0.05)0.2) name(histK) graphregion(color(gs16))
graph export ${resfold}scaleK_`RD_NRD'.pdf, replace name(histK)
hist logV ,bin(100) xmtick(0(1)15) ymtick(0(0.05)0.2) name(histV) graphregion(color(gs16))
graph export ${resfold}scaleV_`RD_NRD'.pdf, replace name(histV)




**** Deflate all vars ****
clear all
local RD_NRD 	"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
xtset , clear

foreach X of varlist K V {
	if ("`X'" == "K") {
		local min`X'    "3"
		local max`X'    "10"
		local minBin 	"50"
		local vars`X'   "SL GP XS OI RD IK IX DP D V"
		//local vars`X'   "dp"
	}
	else {
		local min`X'    "3"
		local max`X'    "10"
		local minBin 	"50"
		local vars`X'   "K SL GP XS OI RD IK IX DP D"
	}
	
	qui: gen bin`X'  = .
	
	// do growth of deflator and forward of deflator
	foreach Y of varlist rF`X' F_`X'  {
		qui: gen iqr`Y' = .
		//bootstrap res1=(_b[`X']) res2=(exp(_b[_cons])) , reps(50) seed(1337) cluster(gvkey) : scalereg `X' bin`X' rF`X' tmpr ,method(iqr) min(`min`X'') max(`max`X'') repeat(2) regmethod(lad) minbin(`minBin')
		scalereg `X' bin`X' `Y' iqr`Y' ,method(iqr) min(`min`X'') max(`max`X'') repeat(2) regmethod(lad) minbin(`minBin')
		disp "Deflate `Y' on `X':"
		disp "Slope: " _b[`X'] " / Icept: " exp(_b[_cons])
		mat res`Y'iqr = e(b)
		qui: bysort bin`X': gen binrep=1 if _n==1  & _N>`minBin'
		qui: gen liqr`Y' = log(iqr`Y')
		qui: gen lp`Y'iqr = res`Y'iqr[1,2] + res`Y'iqr[1,1]*bin`X'
		//twoway (scatter liqr`Y' bin`X' if binrep==1 , msize(tiny) mlcolor(gs7)) (line lp`Y'iqr bin`X' if binrep==1 , lwidth(thin) lcolor(gs0)) ,  xmtick(0(1)15) xline(`min`X'' `max`X'') name("`Y'`X'")
		//graph export ${resfold}scatter`Y'`X'_`RD_NRD'.pdf, replace name(`Y'`X')
		qui: gen `Y'`X' = (`Y'/exp(res`Y'iqr[1,1]*log`X')) if log`X' > `min`X'' & log`X' < `max`X''
		drop iqr`Y' binrep liqr`Y' lp`Y'iqr
	}

	// K SL GP XS OI RD IK IX DP D V
	foreach Y of varlist `vars`X'' {
		qui: gen iqr`Y'  = .

		//bootstrap res1=(_b[`X']) res2=(exp(_b[_cons])) , reps(50) seed(1337) cluster(gvkey) : scalereg `X' bin`X' `Y' iqr`Y' ,method(iqr) min(`min`X'') max(`max`X'') repeat(2) regmethod(lad) minbin(`minBin')
		//scalereg `X' bin`X' `Y' iqr`Y' ,method(iqr) min(`min`X'') max(`max`X'') repeat(2) regmethod(lad) minbin(`minBin')
		scalereg `X' bin`X' `Y' iqr`Y' ,method(median) min(`min`X'') max(`max`X'') repeat(2) regmethod(lad) minbin(`minBin')
		disp "Deflate `Y' on `X':"
		disp "Slope: " _b[`X'] " / Icept: " exp(_b[_cons])
		mat res`Y'iqr = e(b)

		// create scatter
		qui: bysort bin`X': gen binrep=1 if _n==1  & _N>`minBin'
		qui: gen liqr`Y' = log(iqr`Y')
		qui: gen lp`Y'iqr = res`Y'iqr[1,2] + res`Y'iqr[1,1]*bin`X'
		twoway (scatter log`Y'  log`X' , msize(vtiny) mlcolor(gs12)) (scatter liqr`Y' bin`X' if binrep==1 , msize(tiny) mlcolor(gs7)) (line lp`Y'iqr bin`X' if binrep==1 , lwidth(thin) lcolor(gs0)) ,  xmtick(0(1)15) xline(`min`X'' `max`X'') name("`Y'`X'")
		//graph export ${resfold}scatter`Y'`X'_`RD_NRD'.pdf, replace name(`Y'`X')

		// generate deflated variables (and lags,logs,growth)
		qui: gen `Y'`X' = (`Y'/exp(res`Y'iqr[1,1]*log`X')) if log`X' > `min`X'' & log`X' < `max`X''
		qui: gen L_`Y'`X' = (L_`Y'/exp(res`Y'iqr[1,1]*logL_`X')) if log`X' > `min`X'' & log`X' < `max`X''
		qui: gen log`Y'`X' = log(`Y'`X')
		qui: gen logL_`Y'`X' = log(L_`Y'`X')
		qui: gen r`Y'`X' = log`Y'`X' - logL_`Y'`X'
		
		// cleanup
		drop binrep iqr`Y' liqr`Y' lp`Y'iqr
	}
}
save "${datfold}Interim4`RD_NRD'.dta", replace



**** Analyze deflated vars ****
clear all
local RD_NRD 	"BOTH"
use "${datfold}Interim4`RD_NRD'.dta", clear

// concentrate on obs that are within both K and V limits
count
count if rFKK<. & rFVV<.
keep if rFKK<. & rFVV<.

// growth
tabstat rFKK rFVV if rFKK<. & rFVV<., s(mean sd median iqr skewness) f(%8.3f)
//regress rFKK rFKK  if rFKK<. & rFVV<. 
//regress rFVV rFVV  if rFKK<. & rFVV<. 

// base stats
tabstat SLK GPK XSK OIK RDK IKK IXK DK VK if rFKK<. & rFVV<., s(mean sd median iqr skewness) f(%8.3f)
tabstat KV SLV GPV XSV OIV RDV IKV IXV DV if rFKK<. & rFVV<., s(mean sd median iqr skewness) f(%8.3f)

// growth base stats
tabstat rSLK rGPK rXSK rOIK rRDK rIKK rIXK rDK rVK if rFKK<. & rFVV<., s(mean sd median iqr skewness n) f(%8.3f)
tabstat rKV rSLV rGPV rXSV rOIV rRDV rIKV rIXV rDV if rFKK<. & rFVV<., s(mean sd median iqr skewness n) f(%8.3f)

// dynamics - regress log`Y' logL_`Y' if rFKK<. & rFVV<.
qui{
foreach X of varlist K V {
	if ("`X'" == "K") {
		local vars`X'   "SLK GPK XSK OIK RDK IKK IXK DK VK"
	}
	else {
		local vars`X'   "KV SLV GPV XSV OIV RDV IKV IXV DV"
	}

	matrix results`X' = J(8,9,0)
	matrix rownames results`X' = "b0" "b1" "b1_(1-b0)" "MSE" "MSS" "R2" "eMed" "eMode"
	matrix colnames results`X' = `vars`X''

	local curr = 1	
	foreach Y of varlist `vars`X'' {
		regress log`Y' logL_`Y' if rFKK<. & rFVV<. 
		mat tmp1 = e(b)
		mat tmp2 = e(rmse)
		mat tmp3 = e(mss)
		mat tmp4 = e(r2)
		mat tmp5 = e(N)
		
		mat results`X'[1,`curr'] = tmp1[1,1]
		mat results`X'[2,`curr'] = tmp1[1,2]
		mat results`X'[3,`curr'] = tmp1[1,2]/(1-tmp1[1,1])
		mat results`X'[4,`curr'] = tmp2[1,1]^2
		mat results`X'[5,`curr'] = tmp3[1,1]/tmp5[1,1]
		mat results`X'[6,`curr'] = tmp4[1,1]
		
		tabstat log`Y'  if rFKK<. & rFVV<. , s(mean p50 sd) save
		mat tmp = r(StatTotal)
		mat results`X'[7,`curr'] = exp((tmp[1,1] + tmp[2,1])/2)
		mat results`X'[8,`curr'] = results`X'[7,`curr']/exp(tmp[2,1]^2)
		
		local curr = `curr'+1
	}
}
}
mat list resultsK, format(%5.3f)
mat list resultsV, format(%5.3f)

// correlations
corr rFKK rFVV SLK GPK XSK OIK RDK IKK IXK DK VK KV SLV GPV XSV OIV RDV IKV IXV DV if rFKK<. & rFVV<.
corr rFKK rFVV logSLK logGPK logXSK logOIK logRDK logIKK logIXK logDK logVK logKV logSLV logGPV logXSV logOIV logRDV logIKV logIXV logDV if rFKK<. & rFVV<.



**** Get data into matlab ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"NRD"
use "${datfold}Interim3`RD_NRD'.dta", clear
count

foreach vr of varlist K V {
	gen uF_`vr' = r`vr'
	gen u`vr'   = L.r`vr'
	gen uL_`vr' = L2.r`vr'
}

foreach vr of varlist GP XS OI IK RD D {
	gen uF_`vr' = F.r`vr'
	gen u`vr'   = r`vr'
	gen uL_`vr' = L.r`vr'
}

local vrs = "uF_K uK uL_K uF_V uV uL_V uF_GP uGP uL_GP uF_XS uXS uL_XS uF_OI uOI uL_OI uF_IK uIK uL_IK uF_RD uRD uL_RD uF_D uD uL_D"
//keep touse gvkey fyear ind age `vrs'
gen X = 0

foreach vr of varlist `vrs' {
	if ("`vr'"=="uF_RD" | "`vr'"=="uRD"  | "`vr'"=="uL_RD") & "`RD_NRD'" ~= "RD" {
		continue
	}
	replace touse=0 if `vr'>=.
}
count if touse

// gvkey fyear ind age Z L_Z X L_X W L_W F_N N L_N vrs
outsheet gvkey fyear ind age X X X X X X X X X `vrs' using data\data_160412`RD_NRD'.csv if touse, comma replace nonames nolabel



**** Get entry/exit distributions ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count

sort gvkey fyear
by gvkey: gen pos=_n
by gvkey: gen maxpos=_N
gen atminpos = pos==1
gen atmaxpos = pos==maxpos
bysort fyear: gen totN = _N
keep if atminpos | atmaxpos
keep if fyear>1980 & fyear<2010
bysort fyear: egen tminpos = total(atminpos)
bysort fyear: egen tmaxpos = total(atmaxpos)
gen lK = log(K)
gen lV = log(V)
gen lXS = log(XS)
tabstat lK if atminpos, s(mean sd p1 p5 p50 p95 p99)
tabstat lK if atmaxpos, s(mean sd p1 p5 p50 p95 p99)
tabstat lV if atminpos, s(mean sd p1 p5 p50 p95 p99)
tabstat lV if atmaxpos, s(mean sd p1 p5 p50 p95 p99)
gen lQ = lV-lK
gen lO = lXS-lV
tabstat lQ if atminpos, s(mean sd p1 p5 p50 p95 p99)
tabstat lQ if atmaxpos, s(mean sd p1 p5 p50 p95 p99)
hist lQ if atmaxpos,bin(100) name(atmax)
hist lQ if atminpos,bin(100) name(atmin)
tabstat lO if atminpos, s(mean sd p1 p5 p50 p95 p99)
tabstat lO if atmaxpos, s(mean sd p1 p5 p50 p95 p99)

bysort fyear: keep if _n==1
gen tminpos1 = tminpos/totN
gen tmaxpos1 = tmaxpos/totN
tabstat tminpos1 tmaxpos1 , s(mean sd)


bysort fyear: egen avgpos = mean(pos)
sum avgpos
hist avgpos,bin(100)
hist avgpos,bin(10)
list fyear pos avgpos in 1/20
keep if fyear>1980
list fyear pos avgpos in 1/20
gen pos1 = pos==1
bysort fyear tpos1 = total(pos1)
bysort fyear: egen tpos1 = total(pos1)
bysort gvkey: egen posmax = max(pos)
gen posmax = pos==posmax
gen atposmax = pos==posmax
keep if fyear<2010
count
bysort fyear gen totN = _N
bysort fyear: gen totN = _N
bysort fyear: egen tposmax = total(atposmax)
gen percpos1 = tpos1/totN
gen percposmax = tposmax/totN
bysort fyear: keep if _n==1
list fyear tpos1 totN tposmax percpos1 percposmax
tabstat percpos1 percposmax , s(mean)






**** Investment Q regressions ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count

gen VK  = L.rV/L.rK			if touse
gen lVK = log(L.rV/L.rK)	if touse
gen IKK  = rIK/L.rK			if touse
gen lIKK = log(rIK/L.rK)	if touse
gen DPK = rDP/L.rK			if touse
gen lDPK = log(rDP/L.rK)	if touse

regress IKK VK
local b0_1 = _b[_cons]
local b1_1 = _b[VK]

regress lIKK lVK
local b0_2 = _b[_cons]
local b1_2 = _b[lVK]

regress lIKK lVK if lVK>0
local b0_3 = _b[_cons]
local b1_3 = _b[lVK]

gen p1 = (`b0_1' + `b1_1'*VK)
gen p2 = exp(`b0_2' + `b1_2'*lVK)
gen p3 = exp(`b0_3' + `b1_3'*lVK)

gen r1  = IKK - p1
gen r2  = IKK - p2
gen r3  = IKK - p3
gen r12 = r1^2
gen r22 = r2^2
gen r32 = r3^2


tabstat r1 r2 r3, s(p1 p50 p99 mean iqr)
tabstat r12 r22 r32, s(mean)
corr IKK p1 p2 p3

gen lp1  = log(p1)
gen lp2  = log(p2)
gen lp3  = log(p3)

gen lr1  = lIKK - lp1
gen lr2  = lIKK - lp2
gen lr3  = lIKK - lp3
gen lr12 = lr1^2
gen lr22 = lr2^2
gen lr32 = lr3^2

tabstat lr1 lr2 lr3, s(p1 p50 p99 mean iqr)
tabstat lr12 lr22 lr32, s(mean)
corr lIKK lp1 lp2 lp3



**** Depreciation Q regressions ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count

gen VK  = V/K
gen lVK = log(V/K)
gen IKK  = IK/K
gen lIKK = log(IK/K)
gen DPK = DP/K
gen lDPK = log(DP/K)
gen NKK  = (IK-DP)/K
gen GPK = GP/K
gen XSK = XS/K
gen lGPK = log(GPK)
gen lXSK = log(XSK)

egen slGPK = std(lGPK)
egen slDPK = std(lDPK)
egen slXSK = std(lXSK)


keep if lVK>0

egen mlVK = mean(lVK)
gen  rlVK = (lVK-mlVK)^2
egen  SST = total(rlVK)
drop mlVK rlVK

truncreg lVK slGPK slDPK slXSK if lVK>0, vce(r) ll(0)
local b0_1 = _b[_cons]
local b1_1 = _b[slGPK]
local b2_1 = _b[slDPK]
local b3_1 = _b[slXSK]

gen p1  = (`b0_1' + 0.318 + `b1_1'*slGPK + `b2_1'*slDPK + `b3_1'*slXSK)
gen r1  = lVK - p1
gen r12 = r1^2
egen SSR1 = total(r12)
disp (1-SSR1[1]/SST[1])


tabstat lVK r1, s(p1 p50 p99 mean iqr)
tabstat r12 , s(mean)
corr lVK p1 

truncreg lVK lGPK lDPK lXSK if lVK>0, vce(r) ll(0)
local b0_2 = _b[_cons]
local b1_2 = _b[lGPK]
local b2_2 = _b[lDPK]
local b3_2 = _b[lXSK]

gen p2  = (`b0_2' - 1.319 + `b1_2'*slGPK + `b2_2'*slDPK + `b3_2'*slXSK)
gen r2  = lVK - p2
gen r22 = r2^2
egen SSR2 = total(r22)
disp (1-SSR2[1]/SST[1])

tabstat lVK r2, s(p1 p50 p99 mean iqr)
tabstat r22 , s(mean)
corr lVK p2 







regress lIKK lVK
regress lDPK lVK lIKK
regress lVK lDPK

regress IKK VK
local b0_1 = _b[_cons]
local b1_1 = _b[VK]

regress lIKK lVK
local b0_2 = _b[_cons]
local b1_2 = _b[lVK]

regress lIKK lVK if lVK>0
local b0_3 = _b[_cons]
local b1_3 = _b[lVK]

gen p1 = (`b0_1' + `b1_1'*VK)
gen p2 = exp(`b0_2' + `b1_2'*lVK)
gen p3 = exp(`b0_3' + `b1_3'*lVK)

gen r1  = IKK - p1
gen r2  = IKK - p2
gen r3  = IKK - p3
gen r12 = r1^2
gen r22 = r2^2
gen r32 = r3^2


tabstat r1 r2 r3, s(p1 p50 p99 mean iqr)
tabstat r12 r22 r32, s(mean)
corr IKK p1 p2 p3

gen lp1  = log(p1)
gen lp2  = log(p2)
gen lp3  = log(p3)

gen lr1  = lIKK - lp1
gen lr2  = lIKK - lp2
gen lr3  = lIKK - lp3
gen lr12 = lr1^2
gen lr22 = lr2^2
gen lr32 = lr3^2

tabstat lr1 lr2 lr3, s(p1 p50 p99 mean iqr)
tabstat lr12 lr22 lr32, s(mean)
corr lIKK lp1 lp2 lp3








**** Another look at IK ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count

gen IK1 = IK-DP
gen DK  = F_K - K

gen IK1K = IK1/K
gen DKK  = DK/K

scatter IK1K DKK if K>exp(3) & F_K>exp(3) & DKK>-1 & DKK<1 & IK1K>-1 & IK1K<1, msize(vtiny)
















gen XD = xint/(at - (seq+txditc-pstk))
gen DK = (at - (seq+txditc-pstk))/K
gen XSK1 = XS/K
gen GPK1 = GP/K
regress XD DK GPK1 XSK1 I.ind, vce(cl gvkey)

gen IKK1 = IK/K
gen VK1 = L_V/K
gen IKK2 = IK/K
gen VK2 = L.mve/K
gen IKK3 = log(IK/K)
gen VK3 = log(L_V/K)
gen IKK4 = log(IK/K)
gen VK4 = log(L.mve/K)
gen IKK5 = log(IK/K^0.0001)
gen VK5 = log(L_V/K^0.0001)
gen IKK6 = log(IK/K^0.0001)
gen VK6 = log(L.mve/K^0.0001)

regress IKK1 VK1 , vce(cl gvkey)
regress IKK2 VK2 , vce(cl gvkey)
regress IKK3 VK3 , vce(cl gvkey)
regress IKK4 VK4 , vce(cl gvkey)
regress IKK5 VK5 , vce(cl gvkey)
regress IKK6 VK6 , vce(cl gvkey)




**** UNDERSTAND CAPITAL INSTALLATION ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count
gen XI  = xint
gen lGP = log(GP)
gen lXS = log(XS)
gen lDP = log(DP)
gen lXI = log(XI)
gen lK  = log(K)
//qreg lGP lK if K>exp(3) & K<exp(10)
//qreg lXS lK if K>exp(3) & K<exp(10)
//qreg lDP lK if K>exp(3) & K<exp(10)
//qreg lXI lK if K>exp(3) & K<exp(10)
local thetaZ = 0.8185
local thetaX = 0.7340
local thetaD = 0.9277
local thetaI = 0.9601
gen lGPK = log(GP/K^`thetaZ')
gen lXSK = log(XS/K^`thetaX')
gen lDPK = log(DP/K^`thetaD')
gen lXIK = log(XI/K^`thetaI')
//hist lGPK if K>exp(3) & K<exp(10) & GP>1 & XS>1 & DP>1 & XI>1 , name(lGPK)
//hist lXSK if K>exp(3) & K<exp(10) & GP>1 & XS>1 & DP>1 & XI>1 , name(lXSK)
//hist lDPK if K>exp(3) & K<exp(10) & GP>1 & XS>1 & DP>1 & XI>1 , name(lDPK)
//hist lXIK if K>exp(3) & K<exp(10) & GP>1 & XS>1 & DP>1 & XI>1 , name(lXIK)

gen KpK = F_K/K
gen CKpK = (K - DP + IK)/K
gen IKK = IK/K
tabstat CKpK KpK if K>exp(3) & K<exp(10), s(p5 p25 p50 p75 p95 mean)

gen EB = GP - XS - DP - XI
gen IB = ib
gen EBK = EB/K
gen IBK = IB/K
qreg IBK EBK if K>exp(3) & K<exp(10) & EBK>0
qreg IBK EBK if K>exp(3) & K<exp(10) & EBK<0
local tauC = 0.35

//gen EA = (1-`tauC'*(EB>0))*EB + DP + XI
gen EA = (1-`tauC')*EB + DP + XI

gen DVK = (EA - IK)/K
gen DK  = D/K
gen t = DK - DVK
hist t if abs(t)<1 & K>exp(3) & K<exp(10) & DVK>0, bin(100)
hist t if abs(t)<1 & K>exp(3) & K<exp(10) & DVK<0, bin(100)
hist t if abs(t)<1 & K>exp(3) & K<exp(10), bin(100)
tabstat DK DVK t if K>exp(3) & K<exp(10), s(p5 p25 p50 p75 p95 mean)
twoway (scatter t DVK if K>exp(3) & K<exp(10) & abs(t)<2 & abs(DVK)<2 , msize(vtiny)) 
twoway (scatter t DK  if K>exp(3) & K<exp(10) & abs(t)<2 & abs(DK)<2  , msize(vtiny))
twoway (scatter t lK  if K>exp(3) & K<exp(10) & abs(t)<2 , msize(vtiny))
twoway (scatter t EBK if K>exp(3) & K<exp(10) & abs(t)<2 & abs(EBK)<2, msize(vtiny))
twoway (scatter t IKK if K>exp(3) & K<exp(10) & abs(t)<2 & abs(IKK)<1, msize(vtiny))
qreg t DK if K>exp(3) & K<exp(10)
qreg t DVK if K>exp(3) & K<exp(10)

gen noise = rnormal()
replace DVK1 = (DVK-0.082)*exp(0.9*noise) + 0.082
replace DVK1 = (DVK-0.082)*1.4 + 0.082
replace DVK1 = (0.8088*DVK+0.0268)*(DVK>0) + (0.6234*DVK+0.0183)*(DVK<0)
replace DVK1 = (1.0*DVK)*(DVK>0) + (0.65*DVK)*(DVK<0)
replace DVK1 = (DVK)*(EB>0) + (DVK-0.35*EBK)*(EB<0)
replace t1 = DK - DVK1
tabstat DK DVK DVK1 t t1 if K>exp(3) & K<exp(10), s(p5 p25 p50 p75 p95 mean)
twoway (scatter t1 DVK1 if K>exp(3) & K<exp(10) & abs(t1)<2 & abs(DVK1)<2, msize(vtiny)) 
twoway (scatter t1 DK  if K>exp(3) & K<exp(10) & abs(t1)<2 & abs(DK)<2 , msize(vtiny)) //(scatter f3 DK if K>exp(3) & K<exp(10) & abs(t)<2, msize(vtiny))
hist t1 if abs(t1)<1 & K>exp(3) & K<exp(10) & DVK1>0 , bin(100)
hist t1 if abs(t1)<1 & K>exp(3) & K<exp(10) & DVK1<0 , bin(100)
hist t1 if abs(t1)<1 & K>exp(3) & K<exp(10), bin(100)

twoway (scatter DK DVK1 if K>exp(3) & K<exp(10) & abs(DK)<2 & abs(DVK1)<2, msize(vtiny)) (scatter DK DK if K>exp(3) & K<exp(10) & abs(DK)<2, msize(vtiny))




**** Now calculated ones ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
local delK   = 0.1250
local thetaZ = 0.8185
local thetaX = 0.7340
local thetaD = 0.9277
local tauC   = 0.3500

gen DP = `delK'*K^`thetaD'
gen Kp = K-DP+IK
gen EB = GP-XS-DP
gen EA = (1-`tauC'*(EB>0))*EB + DP
gen DV = EA - IK

gen DVK = DV/K
gen DK  = D1/K
tabstat DVK DK if K>exp(3) & K<exp(10), s(p5 p25 p50 p75 p95 mean)

gen IB1 = D1 + IK - DP - XI








gen DP = dp
gen DK = Kp - K + DP
gen DKK1 = DK/K
gen DKK2 = DK/K^0.85
hist DKK1 if DKK1>-2 & DKK1<2 & K>exp(3) & K<exp(10) & IK>0 , name(DKK1)
hist DKK2 if DKK2>-2 & DKK2<2 & K>exp(3) & K<exp(10) & IK>0 , name(DKK2)
gen DKIK = DK/IK - 1
hist DKIK if DKIK>-4 & DKIK<4 & K>exp(3) & K<exp(10) & IK>0 , name(DKIK)
scatter DKIK IK if abs(DKIK)<50 & IK<5000 & K>exp(3) & K<exp(10), msize(vtiny)
scatter DKIK K  if abs(DKIK)<50 & IK<5000 & K>exp(3) & K<exp(10), msize(vtiny)
tabstat DKIK    if abs(DKIK)<50 & IK<5000 & K>exp(3) & K<exp(10), s(p5 p25 p50 p75 p95 mean)

gen lDP = log(DP)
gen lK = log(K)
scatter lDP lK  if K>exp(3) & K<exp(10), msize(vtiny)
regress lDP lK  if K>exp(3) & K<exp(10), vce(r)
gen Cdp1 = 0.1265*K^0.926
gen lCdp1 = log(Cdp1)
scatter lCdp1 lK  if K>exp(3) & K<exp(10), msize(vtiny)
regress lCdp1 lK  if K>exp(3) & K<exp(10)




gen Cdp1 = 0.19*K^0.85
//gen Cdp2 = 0.075*K
gen Cdp2 = 0.080*K
gen Cdpk1 = log(Cdp1/K^0.85)
gen Cdpk2 = log(Cdp2/K^0.85)
gen DPK = log(DP/K^0.85)
tabstat Cdpk1 Cdpk2 DPK if K>exp(3) & K<exp(10) , s(p5 p25 p50 p75 p95 mean)



gen DKK = DK/K 
gen DPK = DP/K 
gen IKK = 1.306*IK/K






**** whats up with dp ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count
gen Cdp1 = 0.19*K^0.85
gen Cdp2 = 0.075*K
gen Cdpk1 = Cdp1/K^0.85
gen Cdpk2 = Cdp2/K^0.85
gen dpk = dp/K^0.85
gen lCdpk1 = log(Cdpk1)
gen lCdpk2 = log(Cdpk2)
gen ldpk = log(dpk)
corr lCdpk1 lCdpk2 ldpk
gen Dldpk1 = ldpk - lCdpk1
gen Dldpk2 = ldpk - lCdpk2
hist Dldpk1 if Dldpk1>-2 & Dldpk1<2, bin(100)
hist Dldpk2 if Dldpk2>-2 & Dldpk2<2, bin(100)

// now we have a calculated dp, what is the investment function?
gen fIK1 = F_K - K + Cdp1
gen fIK2 = F_K - K + Cdp2
gen fIKK1 = fIK1/K^0.85
gen fIKK2 = fIK2/K^0.85
gen IKK = IK/K^0.85
corr fIKK1 fIKK2 IKK
scatter fIKK1 IKK if fIKK1<10 & IKK<10, msize(tiny)
scatter fIKK2 IKK if fIKK2<10 & IKK<10, msize(tiny)
regress fIKK1 IKK
regress fIKK2 IKK

// does that fit the next period capital?
gen CKp1 = K - Cdp1 + 1.053*IK
gen CKp2 = K - Cdp2 + 1.053*IK
gen Kp  = F_K
tabstat CKp1 CKp2 Kp , s(p5 p25 p50 p75 p95 mean)
gen lCKp1 = log(CKp1)
gen lCKp2 = log(CKp2)
gen lKp  = log(Kp)
tabstat lCKp1 lCKp2 lKp , s(p5 p25 p50 p75 p95 mean)
gen DKp1 = lKp - lCKp1
gen DKp2 = lKp - lCKp2
hist DKp1 if DKp1>-0.4 & DKp1<0.4, bin(100)
hist DKp2 if DKp2>-0.4 & DKp2<0.4, bin(100)
tabstat DKp1 DKp2 if K>exp(3) & K<exp(10), s(p5 p25 p50 p75 p95 mean)
count if DKp1>-0.04 & DKp1<0.04 & K>exp(3) & K<exp(10)
count if DKp2>-0.04 & DKp2<0.04 & K>exp(3) & K<exp(10)

// To summarize: A model of dp = 0.075*K and Kp = K - 0.075K + IK is stable




**** whats up with tau ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
count
tabstat sale cogs gp xsga ebitda, s(p25 p50 p75 mean)
gen Cgp = sale-cogs
gen Cebitda = Cgp - xsga
tabstat Cgp Cebitda, s(p25 p50 p75 mean)
tabstat dp ebit xint, s(p25 p50 p75 mean)
gen Cebit = Cebitda-dp
gen Cebt = Cebit-xint
tabstat Cebit Cebt, s(p25 p50 p75 mean)
// prepare deflated values for regression
gen ibk = ib/K^0.75
gen Cebtk = Cebt/K^0.75
// regress above and below zero, without outliers, and find tax rates
regress ibk Cebtk if abs(ibk)<5 & K>100 & Cebtk>0 & Cebtk<6 , nocons
gen tau1 = 1-_b[Cebtk]
regress ibk Cebtk if abs(ibk)<5 & K>100 & Cebtk<0 & Cebtk>-3, nocons
gen tau2 = 1-_b[Cebtk]
twoway  (scatter ibk Cebtk if abs(ibk)<5 & K>100 & Cebtk<6 & Cebtk>-3, msize(vtiny) xline(0)) ///
		(lfit    ibk Cebtk if abs(ibk)<5 & K>100 & Cebtk>0 & Cebtk<6 , estopts(nocons)) ///
		(lfit    ibk Cebtk if abs(ibk)<5 & K>100 & Cebtk<0 & Cebtk>-3, estopts(nocons))
// note tau2 not significantly different from 0, and tau1 precisely at 35%!

		


gen Cebtk2 = Cebtk^2
regress ibk Cebtk Cebtk2 if abs(Cebtk)<5 & abs(ibk)<5 & K>100 , nocons

twoway (scatter ibk Cebtk if abs(Cebtk)<20 & abs(ibk)<20 & K>100, msize(vtiny) xline(0)) (lfit ibk Cebtk if abs(Cebtk)<10 & abs(ibk)<10 & K>100 & Cebtk>0)
twoway (scatter ibk Cebtk if abs(Cebtk)<5 & abs(ibk)<5 & K>100, msize(vtiny) xline(0)) (lfit ibk Cebtk if abs(Cebtk)<5 & abs(ibk)<5 & K>100 & Cebtk>0) (lfit ibk Cebtk if abs(Cebtk)<5 & abs(ibk)<5 & K>100 & Cebtk<0) (qfit ibk Cebtk if abs(Cebtk)<5 & abs(ibk)<5 & K>100)
twoway (scatter ibk Cebtk if abs(Cebtk)<5 & abs(ibk)<5 & K>100, msize(vtiny) xline(0)) (qfit ibk Cebtk if abs(Cebtk)<5 & abs(ibk)<5 & K>100, estopts(nocons))

gen Ctaxes = Cebt*tau1*(Cebt>0)
gen Cib = Cebt-Ctaxes
tabstat ib ni Ctaxes Cib capx sppe, s(p25 p50 p75 mean)

gen Cik = capx - sppe
gen Cdist = Cib+dp+xint-Cik
tabstat Cik Cdist D, s(p25 p50 p75 mean)

gen Cdistk = Cdist/K^0.75
gen DK = D/K^0.75
twoway (scatter DK Cdistk if abs(Cdistk)<5 & abs(DK)<5 & K>1000, msize(vtiny) xline(0)) //(line Cdistk Cdistk if abs(Cdistk)<5 & abs(DK)<5 & K>1000) //(lfit DK Cdistk if abs(Cdistk)<5 & abs(DK)<5 & K>1000) 
regress DK Cdistk if K>100 , nocons



gen D1 = (dvc + dvp + prstkc)
gen D2 = (dvc + dvp + prstkc) + xint
gen D3 = (dvc + dvp + prstkc) + xint - (Dat - Dlt)
gen D4 = (dvc + dvp + prstkc) + xint - (Dat - Dlt) + (Dinvt + Dche)
gen D5 = (dvc + dvp + prstkc) + xint - (Dat - (Dseq+Dtxditc-Dpstk))
gen D6 = (dvc + dvp + prstkc) + xint - (Dat - (Dseq+Dtxditc-Dpstk)) + (Dinvt + Dche)
tabstat Cdist D D1 D2 D3 D4 D5 D6, s(p25 p50 p75 mean)


regress ib Cib, nocons
scatter ib Cib if Cib<10000 & ib<10000 & Cib>-5000 &ib>-5000, msize(vtiny)



gen theta = 0.75
gen DP = dp
gen XI = xint
gen EB = OI - DP - XI
gen NI = ni
gen EBK = EB/K^theta
gen NIK = NI/K^theta
regress NIK EBK if abs(EBK)<10 & abs(NIK)<10 & EBK<-2
regress NIK EBK if abs(EBK)<10 & abs(NIK)<10 & EBK>2  
gen tau1 = .
replace tau1 = 0.00 if EBK<=0
replace tau1 = 0.35 if EBK>0
gen tau2 = 0.35
gen EA  = (1-tau1)*EB
gen EAK = EA/K^theta
//twoway (scatter NIK EBK if abs(EBK)<10 & abs(NIK)<10 , msize(vtiny)) (line EAK EBK if abs(EBK)<10 & abs(NIK)<10 & EBK<0, lwidth(thin) lcolor(gs0)) (line EAK EBK if abs(EBK)<10 & abs(NIK)<10 & EBK>0, lwidth(thin) lcolor(gs0))

gen D2 = (1-tau1)*OI + tau2*DP + tau2*XI - IK
//gen D2 = EA + DP + XI - IK //= (1-tau)*(OI-DP-XI) + DP + XI - IK = (1-tau)*OI + tau*DP + tau*XI - IK
gen D3 = D + tau2*xint //+ tau*DP
gen DK  = D/K^0.75
gen D2K = D2/K^0.75
gen D3K = D3/K^0.75
tabstat D D2 D3 DK D2K D3K , s(mean sd p25 p50 p75)

//hist D2K if abs(D2K)<5, name(D2K)
//hist D3K if abs(D3K)<5 ,name(D3K)




regress EAK EBK if abs(EBK)<10 & abs(EAK)<10 & EBK<-1 , nocons
regress EAK EBK if abs(EBK)<10 & abs(EAK)<10 & EBK>1  , nocons



gen N1K = EAK/EBK - 1
scatter N1K EBK if abs(EBK)<10 & abs(NIK)<10 & abs(N1K)<10 , msize(vtiny)
scatter N1K EBK if abs(EBK)<10 & abs(N1K)<10 , msize(vtiny)
regress N1K EBK if abs(EBK)<10 & abs(N1K)<10 & EBK<-1
gen N2K = ((NI-DP)/K^theta)/EBK - 1
scatter N2K EBK if abs(EBK)<10 & abs(NIK)<10 & abs(N2K)<10 , msize(vtiny)
regress N2K EBK if abs(EBK)<10 & abs(NIK)<10 & abs(N2K)<10 & EBK<-2
regress N2K EBK if abs(EBK)<10 & abs(NIK)<10 & abs(N2K)<10 & EBK>2


regress N1K EBK if abs(EBK)<10 & abs(NIK)<10 &abs(N1K)<10

regress NIK EBK if EBK>=0 & EBK<10 & abs(NIK)<10 , nocons
gen tau = 1-_b[EBK] if EBK>=0
replace tau=tau*0.222
regress NIK EBK if EBK<0 & EBK>-10 & abs(NIK<10) , nocons
replace tau = 1-_b[EBK] if EBK<0
gen EA = (1-tau)*EB
gen D2 = EA + DP + XINT - IK



gen e    = (1-tau)*ebt
gen D2 = e + dp + xint - IK
tabstat ni e D2 D ,s(mean sd p25 p50 p75)
gen D3 = (1-tau)*(gp-xs) + tau*dp -IK 
gen D4 = D + tau*xint
gen D3K = D3/K^0.75
gen D4K = D4/K^0.75
tabstat D3 D4 D3K D4K,s(mean sd p25 p50 p75)
hist D3K if abs(D3K)<5, name(D3K)
hist D4K if abs(D4K)<5 ,name(D4K)




tabstat ebt ni e ebtk nik ek ,s(mean sd p25 p50 p75)

corr ebt ni e ebtk nik ek
gen ep = e + xint
gen epk = ep/K^0.75
gen DK = D/K^0.75
tabstat ni ep D nik epk DK ,s(mean sd p25 p50 p75)
corr epk DK
gen D1 = (dvc + dvp + prstkc) - (Dat - (Dseq+Dtxditc-Dpstk)) + (Dinvt + Dche)
gen D1K = D1/K^0.75
tabstat ni ep D D1 nik epk DK D1K ,s(mean sd p25 p50 p75)
corr nik epk DK D1K

gen ep1  = ep - 0.17*K^0.75
gen epk1 = ep1/K^0.75
tabstat ni ep ep1 D D1 nik epk epk1 DK D1K ,s(mean sd p25 p50 p75)

gen nix = ni+xint
gen nixk = nix/K^0.75
tabstat ni nix ep ep1 D D1 nik nixk epk epk1 DK D1K ,s(mean sd p25 p50 p75)




gen oibxk = oibx/K^0.75
gen ibk = ib/K^0.75
tabstat oibdp dp oib xint oibx ib oibxk ibk ,s(mean sd p25 p50 p75)


gen toib = (1-0.37)*oib
gen oibk = oib/K^0.75
gen ibk = ib/K^0.75
gen DK  = D/K^0.75
tabstat oibdp dp oib ib D toib oibk ibk DK ,s(mean sd p25 p50 p75)
corr oibdp dp oib ib D toib oibk ibk DK
























































***********************
**** THE GRAVEYARD ****
***********************
























**** Get data into matlab ****
// RD sample, NRD sample, or BOTH?
clear all
local RD_NRD 			"BOTH"
use "${datfold}Interim3`RD_NRD'.dta", clear
cls
count

//local vrs = "F_K K L_K F_V V L_V IK L_IK RD L_RD OI L_OI D L_D"
local vrs = "F_K K L_K F_V V L_V IK L_IK GP XS OI L_OI D L_D"
//keep gvkey fyear ind age `vrs'
gen X = 0

foreach vr of varlist `vrs' {
	drop if `vr'>=.
}

count


qui: gen bin  = .
qui: gen iqr  = .
gen IKK = IK/K
gen KpK = F_K/K
gen logKpK = log(KpK)
gen logIKK = log(IKK)
keep if logIKK<.
keep if logIKK<. & logIKK>-7 & logIKK<3 & logKpK>-1 & logKpK<2
replace iqr  = -0.1 + 0.00018*(logIKK+7)^4.28
replace iqr = log(1 - 0.069/0.999 + (0.069^(1-0.999)/0.999)*IKK^0.999 - 0.0001*IK)
replace iqr = log(1 - 0.069/0.999 + (0.069^(1-0.999)/0.999)*IKK^0.999)
replace iqr = log(0.9455 + IKK-0.0005*IK)
//replace iqr = log(0.9455 + IKK - 0.01*IK)
//replace iqr = log(1.7339*IKK^0.8 + 0.939+0.001*K)
twoway (scatter logKpK  logIKK , msize(vtiny) mlcolor(gs12)) (scatter iqr logIKK , msize(vtiny) mlcolor(gs7)) if iqr> -1


scalereg IKK bin KpK iqr ,method(iqr) min(-10) max(5) repeat(1) regmethod(ols) minbin(30)
disp "Deflate `Y' on `X':"
disp "Slope: " _b[IKK] " / Icept: " exp(_b[_cons])
mat resiqr = e(b)

// create scatter
qui: bysort bin: gen binrep=1 if _n==1  & _N>30
qui: gen liqr = log(iqr)
qui: gen lpiqr = resiqr[1,2] + resiqr[1,1]*bin
twoway (scatter logKpK  logIKK , msize(vtiny) mlcolor(gs12)) (scatter liqr bin if binrep==1 , msize(tiny) mlcolor(gs7)) (line lpiqr bin if binrep==1 , lwidth(thin) lcolor(gs0)) ,  xmtick(-10(1)5) xline(-10 5)
//graph export ${resfold}scatter`Y'`X'_`RD_NRD'.pdf, replace name(`Y'`X')




outsheet gvkey fyear ind age X X X X X X `vrs' using data\data_160218`RD_NRD'.csv, comma replace nonames nolabel














// censor, present and save
replace `Y'`X' = . if log`X'<`min`X'' | log`X'>`max`X''
tabstat `Y'`X' , s(p1 p99) save
mat tmpres = r(StatTotal)
replace `Y'`X' = . if `Y'`X'<tmpres[1,1] | `Y'`X'>tmpres[2,1]
outsheet `Y'`X' using "`Y'`X'.csv" , comma replace  nonames nolabel

// specific factors, such as collateral constraint
gen tRDV = RD/V
gen tdpK = dp/K
tabstat tRDV tdpK , s(mean median)
gen tlogRDK = log(RD/K)
gen tlogIK = log(I/K)
gen tlogIRDK = tlogRDK + tlogIK
gen tlogCFK = log(CF/K)
scatter tlogIK tlogCFK if tlogIK>-10 & tlogIK<4, msize(vtiny) xmtick(-10(1)4) xscale(range(-10 4)) name(IK)
scatter tlogRDK tlogCFK if tlogRDK>-10 & tlogRDK<4, msize(vtiny) xmtick(-10(1)4) xscale(range(-10 4)) name(RDK)
scatter tlogIRDK tlogCFK if tlogIRDK>-10 & tlogIRDK<4, msize(vtiny) xmtick(-10(1)4) xscale(range(-10 4)) name(IRDK)

// some simple regressions
regress E V , vce(r)
regress E V , vce(cl fyear)
regress E V , vce(cl gvkey)
regress E V , vce(r) nocons




// export to matlab
save "${datfold}Interim4${RD_NRD}.dta", replace
foreach x of varlist IK RDK OIK EK VK IV RDV OIV EV {
	replace LL_`x' = -9999 if LL_`x'>=.
}
outsheet gvkey fyear ind age L_K L_V `vrs' using data\data_151109${RD_NRD}.csv, comma replace nonames nolabel




// review histograms
clear all
use "${datfold}Interim4.dta", clear
tabstat K, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
replace K = K/250
tabstat K, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
tabstat V, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
replace V = V/250
tabstat V, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
tabstat OI, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
replace OI = OI/250
tabstat OI, s(mean sd p1 p5 p25 median p75 p95 p99) f(%8.3f)
tabstat IK RDK CFK EK VK IV RDV CFV EV, s(mean sd p25 median p75 skewness kurtosis) f(%8.3f)
tabstat rIK rRDK rOIK rEK rVK rIV rRDV rOIV rEV, s(mean sd p25 median p75 skewness kurtosis) f(%8.3f)
hist IK   if IK  >-0.2 & IK  <0.6, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.2 0.6)) yscale(range(0 10)) xmtick(-0.2(0.01)0.6) name(IK)   xtitle("IK   - Data") 
hist rIK  if rIK >-2.5 & rIK <2.5, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rIK)  xtitle("rIK  - Data") 
//hist RDK  if RDK >-0.0 & RDK <2.0, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.0 2.0)) yscale(range(0 10)) xmtick(-0.0(0.01)2.0) name(RDK)  xtitle("RDK  - Data") 
//hist rRDK if rRDK>-2.5 & rRDK<2.5, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rRDK) xtitle("rRDK - Data") 
hist OIK  if OIK >-0.5 & OIK <2.0, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.5 2.0)) yscale(range(0 3 )) xmtick(-0.5(0.10)2.0) name(OIK)  xtitle("OIK  - Data")
hist rOIK if rOIK>-2.5 & rOIK<2.5, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rOIK) xtitle("rOIK - Data")
hist VK   if VK  >-0.0 & VK  <20 , bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.0 20 )) yscale(range(0 .5)) xmtick(-0.0(1.00)20 ) name(VK)   xtitle("VK   - Data") 
hist rVK  if rVK >-2.5 & rIK <2.5, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rVK)  xtitle("rVK  - Data") 


// review histograms - top 1/3 by K
clear all
use "${datfold}Interim4.dta", clear
bysort ind: egen K_u = pctile(K) , p(67)
hist IK   if IK  >-0.2 & IK  <0.6 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.2 0.6)) yscale(range(0 10)) xmtick(-0.2(0.01)0.6) name(IKu)   xtitle("IK   - Data") 
hist rIK  if rIK >-2.5 & rIK <2.5 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rIKu)  xtitle("rIK  - Data") 
hist RDK  if RDK >-0.0 & RDK <2.0 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.0 2.0)) yscale(range(0 10)) xmtick(-0.0(0.01)2.0) name(RDKu)  xtitle("RDK  - Data") 
hist rRDK if rRDK>-2.5 & rRDK<2.5 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rRDKu) xtitle("rRDK - Data") 
hist OIK  if OIK >-0.5 & OIK <2.0 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.5 2.0)) yscale(range(0 3 )) xmtick(-0.5(0.10)2.0) name(OIKu)  xtitle("OIK  - Data")
hist rOIK if rOIK>-2.5 & rOIK<2.5 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rOIKu) xtitle("rOIK - Data")
hist VK   if VK  >-0.0 & VK  <20  & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.0 20 )) yscale(range(0 .5)) xmtick(-0.0(1.00)20 ) name(VKu)   xtitle("VK   - Data") 
hist rVK  if rVK >-2.5 & rIK <2.5 & K>K_u, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rVKu)  xtitle("rVK  - Data") 


// review histograms - bottom 1/3 by K
clear all
use "${datfold}Interim4.dta", clear
bysort ind: egen K_l = pctile(K) , p(33)
hist IK   if IK  >-0.2 & IK  <0.6 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.2 0.6)) yscale(range(0 10)) xmtick(-0.2(0.01)0.6) name(IKl)   xtitle("IK   - Data") 
hist rIK  if rIK >-2.5 & rIK <2.5 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rIKl)  xtitle("rIK  - Data") 
hist RDK  if RDK >-0.0 & RDK <2.0 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.0 2.0)) yscale(range(0 10)) xmtick(-0.0(0.01)2.0) name(RDKl)  xtitle("RDK  - Data") 
hist rRDK if rRDK>-2.5 & rRDK<2.5 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rRDKl) xtitle("rRDK - Data") 
hist OIK  if OIK >-0.5 & OIK <2.0 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.5 2.0)) yscale(range(0 3 )) xmtick(-0.5(0.10)2.0) name(OIKl)  xtitle("OIK  - Data")
hist rOIK if rOIK>-2.5 & rOIK<2.5 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rOIKl) xtitle("rOIK - Data")
hist VK   if VK  >-0.0 & VK  <20  & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.0 20 )) yscale(range(0 .5)) xmtick(-0.0(1.00)20 ) name(VKl)   xtitle("VK   - Data") 
hist rVK  if rVK >-2.5 & rIK <2.5 & K>K_l, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2.5 2.5)) yscale(range(0 4 )) xmtick(-2.5(0.10)2.5) name(rVKl)  xtitle("rVK  - Data") 


// This block shows that firm growth in OIK is positively correlated for positive
// growth firms, but negatively correlated for negative growth firms, so growing firms
// tend to continue to grow, but declining firms regress around a mean growth rate
// The positive correlation result for growing firms only holds in the NRD set,
// whereas in the RD set, coeff is insignificant.
// This also holds when using mean or median LrOIK as the cutoff
clear all
use "${datfold}Interim4.dta", clear
regress rOIK LrOIK if LrOIK>0, vce(cl fyear)
regress rOIK LrOIK if LrOIK<0, vce(cl fyear)






egen K_p33 = pctile(K) ,p(33)
egen K_p66 = pctile(K) ,p(66)

hist IK, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.2 0.6)) xtitle("Physical investment intensity IK/K - entire sample")
graph export ${resfold}hist_IKK.pdf, replace
hist IK if K<K_p33, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.2 0.6)) xtitle("Physical investment intensity IK/K - small tercile")
graph export ${resfold}hist_IKKs.pdf, replace
hist IK if K>K_p66, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-0.2 0.6)) xtitle("Physical investment intensity IK/K - large tercile")
graph export ${resfold}hist_IKKl.pdf, replace

hist RDK if RDK<1.5 , bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(0 1.5)) xtitle("Knowledge investment intensity IN/K - entire sample")
graph export ${resfold}hist_INK.pdf, replace
hist RDK if K<K_p33 & RDK<1.5 , bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(0 1.5)) xtitle("Knowledge investment intensity IN/K - small tercile")
graph export ${resfold}hist_INKs.pdf, replace
hist RDK if K>K_p66 & RDK<1.5 , bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(0 1.5)) xtitle("Knowledge investment intensity IN/K - large tercile")
graph export ${resfold}hist_INKl.pdf, replace

hist OIK if OIK>-2 & OIK<2, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2 2)) xtitle("Operating efficiency OI/K - entire sample")
graph export ${resfold}hist_OIK.pdf, replace
hist OIK if OIK>-2 & OIK<2 & K<K_p33, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2 2)) xtitle("Operating efficiency OI/K - small tercile")
graph export ${resfold}hist_OIKs.pdf, replace
hist OIK if OIK>-2 & OIK<2 & K>K_p66, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-2 2)) xtitle("Operating efficiency OI/K - large tercile")
graph export ${resfold}hist_OIKl.pdf, replace

hist rOIK if rOIK>-1 & rOIK<1, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-1 1)) xtitle("Firm-demeaned operating efficiency OI/K - entire sample")
graph export ${resfold}hist_rOIK.pdf, replace
hist rOIK if rOIK>-1 & rOIK<1 & K<K_p33, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-1 1)) xtitle("Firm-demeaned operating efficiency OI/K - small tercile")
graph export ${resfold}hist_rOIKs.pdf, replace
hist rOIK if rOIK>-1 & rOIK<1 & K>K_p66, bin(100) color(gs8) graphregion(fcolor(gs16)) xscale(axis(1) range(-1 1)) xtitle("Firm-demeaned operating efficiency OI/K - large tercile")
graph export ${resfold}hist_rOIKl.pdf, replace

drop K_p33 K_p66


tabstat IV RDV , s(mean sd median skewness)
corr IV RDV
hist IV
hist RDV
hist OIV

gen logV = log(V)
gen logRDV = log(RDV)
egen dOIV = std(OIV)
egen dRD  = std(RD)
egen dI   = std(I)
egen dIK  = std(IK)
egen dRDV = std(RDV)
xtset gvkey fyear
regress logV l.dRD l.logV
regress logV l.dI l.logV
regress logV l.dOIV l.logV
regress logV l.dOIV l.dI l.dRD l.logV


