**** General initial setup ****
clear all
set matsize 11000

global plc = 1
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

insheet using "${datfold}Data.csv" , comma

ren v1  GVKEY
ren v2  FYEAR
ren v4  A
ren v5  K
ren v7  Z
ren v8  X
ren v9  V
ren v10 IK
ren v11 L_IK
ren v14 CF
ren v15 L_CF
ren v16 D
ren v17 L_D
drop v3 v6 v12 v13

gen logV = log(V)
gen logK = log(K)
gen logIK = log(IK)

tabstat K logK CF IK D V logV, s(mean sd median skewness) f(%8.3f)

qui: gen binK  = .
qui: gen iqrCF = .
qui: gen iqrIK = .

scalereg K binK CF iqrCF ,method(iqr) min(1) max(10) repeat(2) regmethod(lad) minbin(50)
disp "Slope: " _b[K] " / Icept: " exp(_b[_cons])
mat res = e(b)
qui: gen CFK = (CF/exp(res[1,1]*logK)) if logK > 1 & logK < 10

scalereg K binK IK iqrIK ,method(iqr) min(1) max(10) repeat(2) regmethod(lad) minbin(50)
disp "Slope: " _b[K] " / Icept: " exp(_b[_cons])
mat res = e(b)
qui: gen IKK = (IK/exp(res[1,1]*logK)) if logK > 1 & logK < 10

xtset GVKEY FYEAR
regress Z L1.Z
regress X L1.X
