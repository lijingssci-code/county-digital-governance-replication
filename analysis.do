* County digital procurement and firm entry.
* Stata 15 or later. Place data.dta and analysis.do in the same folder.
* Set that folder as the working directory, then run: do analysis.do
* Output files are written to results/.

version 15
clear all
confirm file "data.dta"
set more off
set maxvar 30000
set seed 20261106

global JAE_OUT "results"
capture mkdir "$JAE_OUT"
global JAE_DATA "data.dta"
capture log close
log using "$JAE_OUT/analysis.log",text replace
foreach cmd in ftools reghdfe ppmlhdfe ranktest ivreg2 ivreghdfe boottest {
 capture which `cmd'
 if _rc {
  display as error "Required command not found: `cmd'. Install its SSC or official version and rerun."
  exit 499
 }
}
use "$JAE_DATA",clear
isid county_code year
assert !inlist(floor(county_code/100),5301,3205)
count if sample_main==1
local main_n=r(N)
assert `main_n'>0
capture drop py
egen py=group(province_id year)
capture drop cy
egen cy=group(city_id year)
xtset county_code year
global JAE_S "ln_pop ln_gdppc fiscal_gap ln_retail"
global JAE_C "T $JAE_S"
confirm variable IV_M2000

tempfile reg_results wild_results ar_results controls_results
postfile rr str40 model str24 block str32 outcome str48 term str8 estimator str8 fe str16 cluster double b double se double p double lo double hi double n double counties double clusters double r2within double kp double ar_p using `reg_results',replace
postfile ww str24 model double p double clusters double reps using `wild_results',replace
postfile aa str8 instrument double ar_p0 double ar_p05 str244 ar95_set double n using `ar_results',replace
postfile cb str8 model str32 term double b double se double p using `controls_results',replace
program define jae_record
 syntax, MODEL(string) BLOCK(string) Y(string) TERM(string) [EST(string) FE(string) CL(string)]
 if "`est'"=="" local est "OLS"
 if "`fe'"=="" local fe "py"
 if "`cl'"=="" local cl "city_id"
 local b=_b[`term']
 local se=_se[`term']
 local p=2*ttail(e(df_r),abs(`b'/`se'))
 local cv=invttail(e(df_r),.025)
 local rw=e(r2_within)
 if "`est'"=="PPML" {
  local p=2*normal(-abs(`b'/`se'))
  local cv=invnormal(.975)
  local rw=e(r2_p)
 }
 local kp=.
 local ar=.
 if "`est'"=="IV" {
  * Use the default normal approximation for 2SLS inference.
  local p=2*normal(-abs(`b'/`se'))
  local cv=invnormal(.975)
  local kp=e(widstat)
  local ar=e(arfp)
 }
 tempvar tagcounty
 egen byte `tagcounty'=tag(county_code) if e(sample)
 quietly count if `tagcounty'==1
 local nc=r(N)
 post rr ("`model'") ("`block'") ("`y'") ("`term'") ("`est'") ("`fe'") ("`cl'") (`b') (`se') (`p') (`b'-`cv'*`se') (`b'+`cv'*`se') (e(N)) (`nc') (e(N_clust)) (`rw') (`kp') (`ar')
 if inlist("`model'","B1","B2","B3","B4") {
  foreach v in D T ln_pop ln_gdppc fiscal_gap ln_retail Traw primary_share secondary_share {
   capture local bc=_b[`v']
   if _rc==0 {
    local sc=_se[`v']
    local pc=2*ttail(e(df_r),abs(`bc'/`sc'))
    post cb ("`model'") ("`v'") (`bc') (`sc') (`pc')
   }
  }
 }
 capture drop s_`model'
 gen byte s_`model'=e(sample)
end

* 1. Descriptive statistics, zero procurement and group profiles.
tabstat Y D T NN open_new env_rate city_new prov_new ndcity_new ndprov_new mean_new p90_new procN_new procA_new IV_M2000 if sample_main==1,statistics(n mean sd min p25 p50 p75 max) columns(statistics)
count if sample_main==1 & D==0
display "Share with zero recorded digital procurement = " %8.4f r(N)/`main_n'
foreach group in high_tertiary_share high_fiscal_gap ind_strong dev_zone high_hist_exposure {
 tabstat D T ln_pop ln_gdppc fiscal_gap ln_retail if sample_main==1,by(`group') statistics(n mean sd)
}

* Descriptive statistics.
tempfile desc_results
postfile dd str32 variable double n double mean double sd double min double max using `desc_results',replace
foreach v in Y D T NN ln_pop ln_gdppc fiscal_gap ln_retail countY countT Traw primary_share secondary_share dg_amt_w lngdp_prodsvc_w ln_share_prodsvc sh_en_prodsvc_w ln_en_digital_w ln_en_service_w sh_ex_pollute_w open_new env_rate city_new prov_new mean_new p90_new procN_new procA_new high_tertiary_share high_fiscal_gap ind_strong dev_zone IV_M2000 ndcity_new ndprov_new {
 quietly summarize `v' if sample_main==1
 post dd ("`v'") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max))
}
postclose dd
preserve
use `desc_results',clear
format mean sd min max %21.15g
export delimited "$JAE_OUT/descriptive_statistics.csv",replace
restore

* 2. Baseline regressions controlling for total entry.
* [B1] Baseline with total entry as the only time-varying control.
reghdfe Y D T if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(B1) block(baseline) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [B2] Preferred county and province-year effects; full controls; prefecture clustering.
reghdfe Y D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(B2) block(baseline) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [B3] Add uncapped ln(1 + total annual registrations) from the raw counts.
reghdfe Y D $JAE_C Traw if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(B3) block(baseline) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [B4] Add the checked industry value-added shares on the full baseline sample.
reghdfe Y D $JAE_C primary_share secondary_share if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(B4) block(baseline) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* 3. Fixed-effect and clustering sensitivity.
* [yearFE] County and year fixed effects.
reghdfe Y D $JAE_C if sample_main==1,absorb(county_code year) vce(cluster city_id)
jae_record,model(yearFE) block(inference) y(Y) term(D) est(OLS) fe(year) cl(city_id)

* [cityyearFE] County and prefecture-year fixed effects; retain singleton removal.
reghdfe Y D $JAE_C if sample_main==1,absorb(county_code cy) vce(cluster city_id)
jae_record,model(cityyearFE) block(inference) y(Y) term(D) est(OLS) fe(cy) cl(city_id)

* [provinceCluster] Preferred fixed effects with province-clustered standard errors.
reghdfe Y D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster province_id)
jae_record,model(provinceCluster) block(inference) y(Y) term(D) est(OLS) fe(py) cl(province_id)

* [wild_baseline] Province-level wild bootstrap, 9,999 replications and fixed seed.
reghdfe Y D $JAE_C i.py if sample_main==1,absorb(county_code) vce(cluster province_id)
local province_clusters=e(N_clust)
boottest D,cluster(province_id) reps(9999) seed(20261106) nograph
post ww ("baseline") (r(p)) (`province_clusters') (r(reps))

* 4. Other robustness checks and alternative outcomes.
* [R_totalDigital] Use total digital expenditure in place of per-capita expenditure.
reghdfe Y dg_amt_w $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(R_totalDigital) block(robustness) y(Y) term(dg_amt_w) est(OLS) fe(py) cl(city_id)

* [R_noMuni] Exclude direct-controlled municipalities.
reghdfe Y D $JAE_C if sample_main==1 & municipality==0,absorb(county_code py) vce(cluster city_id)
jae_record,model(R_noMuni) block(robustness) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [R_2017] Restrict estimation to 2017-2024.
reghdfe Y D $JAE_C if sample_main==1 & year>=2017,absorb(county_code py) vce(cluster city_id)
jae_record,model(R_2017) block(robustness) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [R_rawLog] Use uncapped log producer-service entry and uncapped log total entry.
reghdfe Yraw D Traw $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(R_rawLog) block(robustness) y(Yraw) term(D) est(OLS) fe(py) cl(city_id)

* [O_lngdp_prodsvc_w] Outcome: producer-service entry per RMB 100 million of GDP.
reghdfe lngdp_prodsvc_w D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_lngdp_prodsvc_w) block(outcomes) y(lngdp_prodsvc_w) term(D) est(OLS) fe(py) cl(city_id)

* [O_ln_share_prodsvc] Outcome: ln(1 + producer-service entry share).
reghdfe ln_share_prodsvc D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_ln_share_prodsvc) block(outcomes) y(ln_share_prodsvc) term(D) est(OLS) fe(py) cl(city_id)

* [O_sh_en_prodsvc_w] Outcome: producer-service share of all new registrations.
reghdfe sh_en_prodsvc_w D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_sh_en_prodsvc_w) block(outcomes) y(sh_en_prodsvc_w) term(D) est(OLS) fe(py) cl(city_id)

* [O_ln_en_digital_w] Outcome: log digital-service entry.
reghdfe ln_en_digital_w D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_ln_en_digital_w) block(outcomes) y(ln_en_digital_w) term(D) est(OLS) fe(py) cl(city_id)

* [O_ln_en_service_w] Outcome: log total service entry.
reghdfe ln_en_service_w D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_ln_en_service_w) block(outcomes) y(ln_en_service_w) term(D) est(OLS) fe(py) cl(city_id)

* [O_exit] Exit-share equation with one-year-lagged procurement.
reghdfe sh_ex_pollute_w lagD $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_exit) block(exit) y(sh_ex_pollute_w) term(lagD) est(OLS) fe(py) cl(city_id)

* [O_exit_env] Add recorded environmental penalty intensity to the exit equation.
reghdfe sh_ex_pollute_w lagD $JAE_C env_rate if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_exit_env) block(exit) y(sh_ex_pollute_w) term(lagD) est(OLS) fe(py) cl(city_id)

* [O_exit_noTotal] Exit equation without the total-entry control.
reghdfe sh_ex_pollute_w lagD $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(O_exit_noTotal) block(exit) y(sh_ex_pollute_w) term(lagD) est(OLS) fe(py) cl(city_id)

* 5. Broad and narrow procurement measures on a common sample.

quietly reghdfe Y D NN $JAE_C if sample_main==1,absorb(county_code py)
capture drop nested_common
gen byte nested_common=e(sample)
* [NN_broad_common] Broad procurement measure on the common broad/narrow sample.
reghdfe Y D $JAE_C if nested_common==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(NN_broad_common) block(nested_narrow) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [NN_narrow] Narrow governance-keyword measure on the same common sample.
reghdfe Y NN $JAE_C if nested_common==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(NN_narrow) block(nested_narrow) y(Y) term(NN) est(OLS) fe(py) cl(city_id)

tempfile magnitudes
postfile mm str32 measure double b double sd double percent_change double n using `magnitudes',replace
foreach measure in D NN {
 quietly reghdfe Y `measure' $JAE_C if nested_common==1,absorb(county_code py) vce(cluster city_id)
 local coef=_b[`measure']
 local n=e(N)
 quietly summarize `measure' if e(sample)
 local sd=r(sd)
 post mm ("`measure'") (`coef') (`sd') (100*(exp(`coef'*`sd')-1)) (`n')
}
postclose mm
preserve
use `magnitudes',clear
format b sd percent_change %21.15g
export delimited "$JAE_OUT/narrow_magnitudes.csv",replace
restore

* 6. Poisson pseudo-maximum likelihood.
* [P_raw_raw] Original entrant counts with uncapped log total entry.
ppmlhdfe countY D Traw $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id) tolerance(1e-10)
jae_record,model(P_raw_raw) block(PPML) y(countY) term(D) est(PPML) fe(py) cl(city_id)

* 7. Channel outcome equations and measurement comparisons.
* [C_open_new] Preferred outcome equation for open-tender value share.
reghdfe open_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_open_new) block(channel) y(open_new) term(D) est(OLS) fe(py) cl(city_id)

* [C0_open_new] Outcome equation for open-tender value share, excluding the total-entry control.
reghdfe open_new D $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C0_open_new) block(channel_noTotal) y(open_new) term(D) est(OLS) fe(py) cl(city_id)

* [CL_open_new] Earlier open-tender value share measure on the available comparison sample.
reghdfe old_open_amt D $JAE_C if s_C_open_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CL_open_new) block(legacy_common) y(old_open_amt) term(D) est(OLS) fe(py) cl(city_id)

* [CC_open_new] Rebuilt open-tender value share measure on the preceding matched sample.
reghdfe open_new D $JAE_C if s_CL_open_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CC_open_new) block(rebuilt_common) y(open_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_env_rate] Preferred outcome equation for recorded penalty intensity.
reghdfe env_rate D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_env_rate) block(channel) y(env_rate) term(D) est(OLS) fe(py) cl(city_id)

* [C0_env_rate] Outcome equation for recorded penalty intensity, excluding the total-entry control.
reghdfe env_rate D $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C0_env_rate) block(channel_noTotal) y(env_rate) term(D) est(OLS) fe(py) cl(city_id)

* [CL_env_rate] Earlier recorded penalty intensity measure on the available comparison sample.
reghdfe env_rate D $JAE_C if s_C_env_rate==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CL_env_rate) block(legacy_common) y(env_rate) term(D) est(OLS) fe(py) cl(city_id)

* [CC_env_rate] Rebuilt recorded penalty intensity measure on the preceding matched sample.
reghdfe env_rate D $JAE_C if s_CL_env_rate==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CC_env_rate) block(rebuilt_common) y(env_rate) term(D) est(OLS) fe(py) cl(city_id)

* [C_city_new] Preferred outcome equation for supplier-source cities.
reghdfe city_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_city_new) block(channel) y(city_new) term(D) est(OLS) fe(py) cl(city_id)

* [C0_city_new] Outcome equation for supplier-source cities, excluding the total-entry control.
reghdfe city_new D $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C0_city_new) block(channel_noTotal) y(city_new) term(D) est(OLS) fe(py) cl(city_id)

* [CL_city_new] Earlier supplier-source cities measure on the available comparison sample.
reghdfe old_src_city D $JAE_C if s_C_city_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CL_city_new) block(legacy_common) y(old_src_city) term(D) est(OLS) fe(py) cl(city_id)

* [CC_city_new] Rebuilt supplier-source cities measure on the preceding matched sample.
reghdfe city_new D $JAE_C if s_CL_city_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CC_city_new) block(rebuilt_common) y(city_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_prov_new] Preferred outcome equation for supplier-source provinces.
reghdfe prov_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_prov_new) block(channel) y(prov_new) term(D) est(OLS) fe(py) cl(city_id)

* [C0_prov_new] Outcome equation for supplier-source provinces, excluding the total-entry control.
reghdfe prov_new D $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C0_prov_new) block(channel_noTotal) y(prov_new) term(D) est(OLS) fe(py) cl(city_id)

* [CL_prov_new] Earlier supplier-source provinces measure on the available comparison sample.
reghdfe old_src_prov D $JAE_C if s_C_prov_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CL_prov_new) block(legacy_common) y(old_src_prov) term(D) est(OLS) fe(py) cl(city_id)

* [CC_prov_new] Rebuilt supplier-source provinces measure on the preceding matched sample.
reghdfe prov_new D $JAE_C if s_CL_prov_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CC_prov_new) block(rebuilt_common) y(prov_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_mean_new] Preferred outcome equation for mean supplier distance.
reghdfe mean_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_mean_new) block(channel) y(mean_new) term(D) est(OLS) fe(py) cl(city_id)

* [C0_mean_new] Outcome equation for mean supplier distance, excluding the total-entry control.
reghdfe mean_new D $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C0_mean_new) block(channel_noTotal) y(mean_new) term(D) est(OLS) fe(py) cl(city_id)

* [CL_mean_new] Earlier mean supplier distance measure on the available comparison sample.
reghdfe old_dist_mean D $JAE_C if s_C_mean_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CL_mean_new) block(legacy_common) y(old_dist_mean) term(D) est(OLS) fe(py) cl(city_id)

* [CC_mean_new] Rebuilt mean supplier distance measure on the preceding matched sample.
reghdfe mean_new D $JAE_C if s_CL_mean_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CC_mean_new) block(rebuilt_common) y(mean_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_p90_new] Preferred outcome equation for P90 supplier distance.
reghdfe p90_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_p90_new) block(channel) y(p90_new) term(D) est(OLS) fe(py) cl(city_id)

* [C0_p90_new] Outcome equation for P90 supplier distance, excluding the total-entry control.
reghdfe p90_new D $JAE_S if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C0_p90_new) block(channel_noTotal) y(p90_new) term(D) est(OLS) fe(py) cl(city_id)

* [CL_p90_new] Earlier P90 supplier distance measure on the available comparison sample.
reghdfe old_dist_p90 D $JAE_C if s_C_p90_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CL_p90_new) block(legacy_common) y(old_dist_p90) term(D) est(OLS) fe(py) cl(city_id)

* [CC_p90_new] Rebuilt P90 supplier distance measure on the preceding matched sample.
reghdfe p90_new D $JAE_C if s_CL_p90_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(CC_p90_new) block(rebuilt_common) y(p90_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_ndcity] Non-digital supplier-source cities as the outcome.
reghdfe ndcity_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_ndcity) block(nondigital) y(ndcity_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_ndprov] Non-digital supplier-source provinces as the outcome.
reghdfe ndprov_new D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_ndprov) block(nondigital) y(ndprov_new) term(D) est(OLS) fe(py) cl(city_id)

* [C_open_complete] Open-tender share in county-years with complete method records.
reghdfe open_complete D $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_open_complete) block(method_sensitivity) y(open_complete) term(D) est(OLS) fe(py) cl(city_id)

* [C_open_samecomplete] Main openness measure on the preceding complete-method sample.
reghdfe open_new D $JAE_C if s_C_open_complete==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(C_open_samecomplete) block(method_sensitivity) y(open_new) term(D) est(OLS) fe(py) cl(city_id)

* [PS_env_rate] Add procurement contract count and value to the recorded penalty intensity equation.
reghdfe env_rate D $JAE_C procN_new procA_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS_env_rate) block(scale_adjusted) y(env_rate) term(D) est(OLS) fe(py) cl(city_id)

* [PS0_env_rate] Omit procurement-scale controls on the same recorded penalty intensity sample.
reghdfe env_rate D $JAE_C if s_PS_env_rate==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS0_env_rate) block(scale_commonbase) y(env_rate) term(D) est(OLS) fe(py) cl(city_id)

* [PS_p90_new] Add procurement contract count and value to the P90 supplier distance equation.
reghdfe p90_new D $JAE_C procN_new procA_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS_p90_new) block(scale_adjusted) y(p90_new) term(D) est(OLS) fe(py) cl(city_id)

* [PS0_p90_new] Omit procurement-scale controls on the same P90 supplier distance sample.
reghdfe p90_new D $JAE_C if s_PS_p90_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS0_p90_new) block(scale_commonbase) y(p90_new) term(D) est(OLS) fe(py) cl(city_id)

* [PS_ndcity_new] Add procurement contract count and value to the non-digital supplier-source cities equation.
reghdfe ndcity_new D $JAE_C procN_new procA_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS_ndcity_new) block(scale_adjusted) y(ndcity_new) term(D) est(OLS) fe(py) cl(city_id)

* [PS0_ndcity_new] Omit procurement-scale controls on the same non-digital supplier-source cities sample.
reghdfe ndcity_new D $JAE_C if s_PS_ndcity_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS0_ndcity_new) block(scale_commonbase) y(ndcity_new) term(D) est(OLS) fe(py) cl(city_id)

* [PS_ndprov_new] Add procurement contract count and value to the non-digital supplier-source provinces equation.
reghdfe ndprov_new D $JAE_C procN_new procA_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS_ndprov_new) block(scale_adjusted) y(ndprov_new) term(D) est(OLS) fe(py) cl(city_id)

* [PS0_ndprov_new] Omit procurement-scale controls on the same non-digital supplier-source provinces sample.
reghdfe ndprov_new D $JAE_C if s_PS_ndprov_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(PS0_ndprov_new) block(scale_commonbase) y(ndprov_new) term(D) est(OLS) fe(py) cl(city_id)

* 8. Individual and joint channel covariates in the entry equation.
* [M1_open_new] Add open-tender value share to the entry equation; record procurement and channel coefficients.
reghdfe Y D $JAE_C open_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_open_new) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_open_new) block(mediator_term) y(Y) term(open_new)
* [M0_open_new] Baseline on the same sample as the entry equation with open-tender value share.
reghdfe Y D $JAE_C if s_M1_open_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_open_new) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [M1_env_rate] Add recorded penalty intensity to the entry equation; record procurement and channel coefficients.
reghdfe Y D $JAE_C env_rate if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_env_rate) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_env_rate) block(mediator_term) y(Y) term(env_rate)
* [M0_env_rate] Baseline on the same sample as the entry equation with recorded penalty intensity.
reghdfe Y D $JAE_C if s_M1_env_rate==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_env_rate) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [M1_city_new] Add supplier-source cities to the entry equation; record procurement and channel coefficients.
reghdfe Y D $JAE_C city_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_city_new) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_city_new) block(mediator_term) y(Y) term(city_new)
* [M0_city_new] Baseline on the same sample as the entry equation with supplier-source cities.
reghdfe Y D $JAE_C if s_M1_city_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_city_new) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [M1_prov_new] Add supplier-source provinces to the entry equation; record procurement and channel coefficients.
reghdfe Y D $JAE_C prov_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_prov_new) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_prov_new) block(mediator_term) y(Y) term(prov_new)
* [M0_prov_new] Baseline on the same sample as the entry equation with supplier-source provinces.
reghdfe Y D $JAE_C if s_M1_prov_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_prov_new) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [M1_mean_new] Add mean supplier distance to the entry equation; record procurement and channel coefficients.
reghdfe Y D $JAE_C mean_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_mean_new) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_mean_new) block(mediator_term) y(Y) term(mean_new)
* [M0_mean_new] Baseline on the same sample as the entry equation with mean supplier distance.
reghdfe Y D $JAE_C if s_M1_mean_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_mean_new) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [M1_p90_new] Add P90 supplier distance to the entry equation; record procurement and channel coefficients.
reghdfe Y D $JAE_C p90_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_p90_new) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_p90_new) block(mediator_term) y(Y) term(p90_new)
* [M0_p90_new] Baseline on the same sample as the entry equation with P90 supplier distance.
reghdfe Y D $JAE_C if s_M1_p90_new==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_p90_new) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [M1_joint] Include all six channel covariates in the entry equation.
reghdfe Y D $JAE_C open_new env_rate city_new prov_new mean_new p90_new if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M1_joint) block(mediator_added) y(Y) term(D) est(OLS) fe(py) cl(city_id)

jae_record,model(M1_joint) block(mediator_term) y(Y) term(open_new)
jae_record,model(M1_joint) block(mediator_term) y(Y) term(env_rate)
jae_record,model(M1_joint) block(mediator_term) y(Y) term(city_new)
jae_record,model(M1_joint) block(mediator_term) y(Y) term(prov_new)
jae_record,model(M1_joint) block(mediator_term) y(Y) term(mean_new)
jae_record,model(M1_joint) block(mediator_term) y(Y) term(p90_new)
* [M0_joint] Baseline on the joint channel model's estimation sample.
reghdfe Y D $JAE_C if s_M1_joint==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(M0_joint) block(mediator_base) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* 9. Subgroup estimates and pooled interaction tests.
* [H_high_tertiary_share_0] Subgroup estimate for pre-period service foundation: group 0.
reghdfe Y D $JAE_C if sample_main==1 & high_tertiary_share==0,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_high_tertiary_share_0) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [H_high_tertiary_share_1] Subgroup estimate for pre-period service foundation: group 1.
reghdfe Y D $JAE_C if sample_main==1 & high_tertiary_share==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_high_tertiary_share_1) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [I_high_tertiary_share] Pooled procurement interaction with pre-period service foundation; test the group difference.
reghdfe Y c.D##i.high_tertiary_share $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(I_high_tertiary_share) block(interaction) y(Y) term(1.high_tertiary_share#c.D) est(OLS) fe(py) cl(city_id)

test 1.high_tertiary_share#c.D
lincom D+1.high_tertiary_share#c.D
* [H_high_fiscal_gap_0] Subgroup estimate for pre-period fiscal pressure: group 0.
reghdfe Y D $JAE_C if sample_main==1 & high_fiscal_gap==0,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_high_fiscal_gap_0) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [H_high_fiscal_gap_1] Subgroup estimate for pre-period fiscal pressure: group 1.
reghdfe Y D $JAE_C if sample_main==1 & high_fiscal_gap==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_high_fiscal_gap_1) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [I_high_fiscal_gap] Pooled procurement interaction with pre-period fiscal pressure; test the group difference.
reghdfe Y c.D##i.high_fiscal_gap $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(I_high_fiscal_gap) block(interaction) y(Y) term(1.high_fiscal_gap#c.D) est(OLS) fe(py) cl(city_id)

test 1.high_fiscal_gap#c.D
lincom D+1.high_fiscal_gap#c.D
* [H_ind_strong_0] Subgroup estimate for industrial strength: group 0.
reghdfe Y D $JAE_C if sample_main==1 & ind_strong==0,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_ind_strong_0) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [H_ind_strong_1] Subgroup estimate for industrial strength: group 1.
reghdfe Y D $JAE_C if sample_main==1 & ind_strong==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_ind_strong_1) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [I_ind_strong] Pooled procurement interaction with industrial strength; test the group difference.
reghdfe Y c.D##i.ind_strong $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(I_ind_strong) block(interaction) y(Y) term(1.ind_strong#c.D) est(OLS) fe(py) cl(city_id)

test 1.ind_strong#c.D
lincom D+1.ind_strong#c.D
* [H_dev_zone_0] Subgroup estimate for 2018 development-zone status: group 0.
reghdfe Y D $JAE_C if sample_main==1 & dev_zone==0,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_dev_zone_0) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [H_dev_zone_1] Subgroup estimate for 2018 development-zone status: group 1.
reghdfe Y D $JAE_C if sample_main==1 & dev_zone==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(H_dev_zone_1) block(subgroup) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* [I_dev_zone] Pooled procurement interaction with 2018 development-zone status; test the group difference.
reghdfe Y c.D##i.dev_zone $JAE_C if sample_main==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(I_dev_zone) block(interaction) y(Y) term(1.dev_zone#c.D) est(OLS) fe(py) cl(city_id)

test 1.dev_zone#c.D
lincom D+1.dev_zone#c.D
* 10. Single-instrument OLS comparison, 2SLS and first stage.
quietly reghdfe Y D IV_M2000 $JAE_C if sample_main==1,absorb(county_code py)
capture drop commonIV
gen byte commonIV=e(sample)
count if commonIV==1
display "IV estimation sample size = " r(N)
tabstat D IV_M2000 T ln_pop ln_gdppc fiscal_gap ln_retail if sample_main==1,by(commonIV) statistics(n mean sd)

* Observable sample profiles.
tempfile profile_results
postfile gg str32 grouping byte group str32 variable double n double mean double sd using `profile_results',replace
foreach grouping in high_tertiary_share high_fiscal_gap ind_strong dev_zone high_hist_exposure commonIV {
 forvalues group=0/1 {
  foreach v in ln_pop ln_gdppc fiscal_gap ln_retail T D {
   quietly summarize `v' if sample_main==1 & `grouping'==`group'
   post gg ("`grouping'") (`group') ("`v'") (r(N)) (r(mean)) (r(sd))
  }
 }
}
postclose gg
preserve
use `profile_results',clear
format mean sd %21.15g
export delimited "$JAE_OUT/group_profiles.csv",replace
restore

* OLS, first stage and 2SLS use the same IV sample.
* Use the supplied instrument without reconstructing its components.
reghdfe Y D $JAE_C if commonIV==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(IV_OLS_same) block(iv) y(Y) term(D) est(OLS) fe(py) cl(city_id)

* First stage.
reghdfe D IV_M2000 $JAE_C if commonIV==1,absorb(county_code py) vce(cluster city_id)
jae_record,model(IV_first) block(iv_first) y(D) term(IV_M2000) est(OLS) fe(py) cl(city_id)
test IV_M2000

* 2SLS and first-stage diagnostics.
ivreghdfe Y $JAE_C (D=IV_M2000) if commonIV==1,absorb(county_code py) cluster(city_id) first
jae_record,model(IV_2SLS) block(iv) y(Y) term(D) est(IV) fe(py) cl(city_id)

* 11. Cluster-robust Anderson-Rubin confidence set.
* Invert the AR test and verify it at beta = 0 and 0.05.
* AR inference assumes instrument validity.
foreach z in IV_M2000 {
 preserve
 keep if commonIV==1
 reghdfe Y `z' $JAE_C,absorb(county_code py) vce(cluster city_id) residuals(ar_y)
 local cy=_b[`z']
 local sy=_se[`z']
 local df=e(df_r)
 local p0=2*ttail(e(df_r),abs(_b[`z']/_se[`z']))
 reghdfe D `z' $JAE_C,absorb(county_code py) vce(cluster city_id) residuals(ar_x)
 local cx=_b[`z']
 quietly reghdfe `z' $JAE_C,absorb(county_code py) residuals(ar_z)
 gen double ar_z2=ar_z^2
 quietly summarize ar_z2,meanonly
 scalar jae_q=r(sum)
 gen double ar_scy=ar_z*ar_y
 gen double ar_scx=ar_z*ar_x
 bysort city_id: egen double ar_gsy=total(ar_scy)
 bysort city_id: egen double ar_gsx=total(ar_scx)
 by city_id: gen byte ar_tag=(_n==1)
 gen double ar_yy=ar_gsy^2 if ar_tag
 gen double ar_xx=ar_gsx^2 if ar_tag
 gen double ar_yx=ar_gsy*ar_gsx if ar_tag
 quietly summarize ar_yy,meanonly
 scalar jae_vyy=r(sum)/jae_q^2
 quietly summarize ar_xx,meanonly
 scalar jae_vxx=r(sum)/jae_q^2
 quietly summarize ar_yx,meanonly
 scalar jae_vyx=r(sum)/jae_q^2
 scalar jae_scale=(`sy')^2/jae_vyy
 scalar jae_vyy=jae_vyy*jae_scale
 scalar jae_vxx=jae_vxx*jae_scale
 scalar jae_vyx=jae_vyx*jae_scale
 scalar jae_crit=invFtail(1,`df',.05)
 scalar jae_A=(`cx')^2-jae_crit*jae_vxx
 scalar jae_B=-2*`cy'*`cx'+2*jae_crit*jae_vyx
 scalar jae_C=(`cy')^2-jae_crit*jae_vyy
 scalar jae_disc=jae_B^2-4*jae_A*jae_C
 local arset ""
 if abs(jae_A)<1e-20 {
  if abs(jae_B)<1e-20 local arset=cond(jae_C<=0,"(-inf,+inf)","empty")
  else if jae_B>0 local arset="(-inf,"+string(-jae_C/jae_B,"%12.8f")+"]"
  else local arset="["+string(-jae_C/jae_B,"%12.8f")+",+inf)"
 }
 else if jae_disc<0 local arset=cond(jae_A<0,"(-inf,+inf)","empty")
 else {
  scalar jae_r1=(-jae_B-sqrt(jae_disc))/(2*jae_A)
  scalar jae_r2=(-jae_B+sqrt(jae_disc))/(2*jae_A)
  scalar jae_low=min(jae_r1,jae_r2)
  scalar jae_high=max(jae_r1,jae_r2)
  if jae_A>0 local arset="["+string(jae_low,"%12.8f")+","+string(jae_high,"%12.8f")+"]"
  else local arset="(-inf,"+string(jae_low,"%12.8f")+"] U ["+string(jae_high,"%12.8f")+",+inf)"
 }
 scalar jae_ar_p0=Ftail(1,`df',(`cy')^2/jae_vyy)
 scalar jae_ar_p05=Ftail(1,`df',(`cy'-.05*`cx')^2/(jae_vyy-.1*jae_vyx+.05^2*jae_vxx))
 gen double ar_y05=Y-.05*D
 reghdfe ar_y05 `z' $JAE_C,absorb(county_code py) vce(cluster city_id)
 local p05=2*ttail(e(df_r),abs(_b[`z']/_se[`z']))
 assert abs(jae_ar_p0-`p0')<2e-5
 assert abs(jae_ar_p05-`p05')<2e-5
 post aa ("`z'") (jae_ar_p0) (jae_ar_p05) ("`arset'") (e(N))
 gen double first_x=ar_x+`cx'*ar_z
 xtile ar_bin=ar_z,nq(20)
 collapse (mean) ar_z first_x (count) bin_n=first_x,by(ar_bin)
 twoway (scatter first_x ar_z [aw=bin_n]) (lfit first_x ar_z),xtitle("IV (residualized)") ytitle("Digital procurement (residualized)") legend(off) graphregion(color(white)) plotregion(color(white))
 graph export "$JAE_OUT/Figure_S1_stata.png",width(2400) replace
 export delimited "$JAE_OUT/first_stage_bins.csv",replace
 restore
}

* 12. Annual first stages: province fixed effects and prefecture clustering.
tempfile year_results
postfile yr int year double b double se double p double lo double hi double n double counties double clusters double r2 using `year_results',replace
forvalues yr=2016/2024 {
 reghdfe D IV_M2000 $JAE_C if commonIV==1 & year==`yr',absorb(province_id) vce(cluster city_id)
 local b=_b[IV_M2000]
 local se=_se[IV_M2000]
 local pv=2*ttail(e(df_r),abs(`b'/`se'))
 local cv=invttail(e(df_r),.025)
 egen byte tag=tag(county_code) if e(sample)
 quietly count if tag==1
 local nc=r(N)
 drop tag
 post yr (`yr') (`b') (`se') (`pv') (`b'-`cv'*`se') (`b'+`cv'*`se') (e(N)) (`nc') (e(N_clust)) (e(r2_within))
}
postclose yr
preserve
use `year_results',clear
format b se p lo hi %21.15g
export delimited "$JAE_OUT/annual_first_stage.csv",replace
restore

* 13. Export all results, estimation samples and AR confidence sets.
preserve
keep county_code year sample_main commonIV s_*
export delimited "$JAE_OUT/model_samples.csv",replace
restore
postclose rr
postclose ww
postclose aa
postclose cb
use `reg_results',clear
format b se p lo hi %21.15g
export delimited "$JAE_OUT/regressions.csv",replace
use `wild_results',clear
export delimited "$JAE_OUT/wild_tests.csv",replace
use `ar_results',clear
export delimited "$JAE_OUT/ar_confidence_sets.csv",replace
use `controls_results',clear
format b se p %21.15g
export delimited "$JAE_OUT/baseline_controls.csv",replace
use "$JAE_DATA",clear
display "All regressions and diagnostics completed."
display "CLEAN_ANALYSIS_COMPLETE"
log close

