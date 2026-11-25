*==============================================================================
* 3_analyze.do   --   STAGE 3 of the replication pack (analysis stage)
*
* "Fitting dynamic network actor models to thousands of team-matches: pass frequency and receiver
*  choice across match contexts in European football"
*
* All tables and figures of the paper: descriptives, Poisson models of pass frequency, observed shares,
* pass geometry, tempo and goal difference, the one-stage conditional logit, random-effects pooling of
* the team-match estimates with club-clustered variance, moderators, Figures 1 and 2.
*
* Inputs (nothing from stages 1-2 is re-estimated):
*   output/build/ (all csv files)                                     (stage 1)
*   output/fit/ ev_<league>, act_<league>, act_<league>_prior (csv)   (stage 2)
*   output/fit/est_<league>_<tag>.csv                                 (team-match nwdynam fits)
*   output/fit/variants/prior_<variant>/est_<league>_VP.csv           (variant fits)
*   data/raw/matches/matches_<league>.json                            (final scores; only for Table S7)
*   To reuse the stored estimates instead of fitting, copy stored/estimates/ into output/fit/.
*
* Outputs: $OUT = <root>/output/analyze/ : one CSV per paper table, Figure1/2
*          as .png/.eps/.tif;  $OUT/tmp holds intermediate .dta files.
*
* Usage:  stata-mp -b do 3_analyze.do                  (all blocks)
*         stata-mp -b do 3_analyze.do pool moderators  (only the named blocks)
*   blocks: desc rate observed geometry outcome clogit pool moderators figures simplectx
*   Block "simplectx" (Table S11) needs only the variant estimates of 2_fit.do (tags VP_base, VP_simple).
*   Block "figures" and block "moderators" reload results of "rate", "observed",
*   "pool" from $OUT/tmp, so run those first (the default run does everything).
* Root folder: environment variable PACK_ROOT, else the default path below.
*==============================================================================
version 17
clear all
set more off
set linesize 200
args only1 only2 only3 only4 only5 only6 only7 only8 only9
local only "`only1' `only2' `only3' `only4' `only5' `only6' `only7' `only8' `only9'"
local only = trim("`only'")
timer clear
timer on 99

local root : environment PACK_ROOT
if "`root'"=="" local root "`c(pwd)'/.."
global ROOT "`root'"
global B   : environment BUILD_DIR
if "$B" == "" global B "$ROOT/output/build"      // output of 1_build.do
global M   : environment FIT_DIR
if "$M" == "" global M "$ROOT/output/fit"        // output of 2_fit.do (estimates and model inputs)
global V   "$M/variants"
global RAW "$ROOT/data/raw"
global OUT "$ROOT/output/analyze"
global TMP "$OUT/tmp"
cap mkdir "$ROOT/output"
cap mkdir "$OUT"
cap mkdir "$TMP"
global FIVE "England Spain Italy Germany France"

* does block X run?   ($RUN_x = 1/0)
foreach b in desc rate observed geometry outcome clogit pool moderators figures simplectx {
    global RUN_`b' = ("`only'"=="") | (strpos(" `only' "," `b' ")>0)
}

*==============================================================================
* 0. MATA LIBRARY
*==============================================================================
mata:
mata clear
mata set matastrict off

// group sums of the columns of a within the levels of key
real matrix grpsum(real matrix a, real colvector key)
{
    real colvector p, ks
    real matrix as, info, out
    real scalar i
    p = order(key, 1); ks = key[p]; as = a[p, .]
    info = panelsetup(ks, 1)
    out = J(rows(info), cols(a), .)
    for (i = 1; i <= rows(info); i++) out[i, .] = colsum(as[|info[i,1],1 \ info[i,2],.|])
    return(out)
}

// DerSimonian-Laird random-effects pooling with club-clustered variance
// (equation (3) of the paper).  v = (b, se, team, match).
// returns (theta, se_naive, se_team, se_twoway, tau, k, G_teams)
real rowvector re_pool(real matrix v)
{
    real colvector b, s, w, wr, tm, mt, a
    real scalar bf, Q, k, c, tau2, th, S, vt, vm, vi, Gt, Gm, Gi, v2
    real matrix g
    b = v[,1]; s = v[,2]; tm = v[,3]; mt = v[,4]
    w = 1:/s:^2; bf = sum(w:*b)/sum(w)
    Q = sum(w:*(b:-bf):^2); k = rows(b); c = sum(w) - sum(w:^2)/sum(w)
    tau2 = max((0, (Q-(k-1))/c))
    wr = 1:/(s:^2 :+ tau2); S = sum(wr); th = sum(wr:*b)/S
    a = wr:*(b:-th)
    g = grpsum(a, tm);             Gt = rows(g); vt = sum(g:^2)/S^2*(Gt>1 ? Gt/(Gt-1) : 1)
    g = grpsum(a, mt);             Gm = rows(g); vm = sum(g:^2)/S^2*(Gm>1 ? Gm/(Gm-1) : 1)
    g = grpsum(a, tm:*1e8 :+ mt);  Gi = rows(g); vi = sum(g:^2)/S^2*(Gi>1 ? Gi/(Gi-1) : 1)
    v2 = vt + vm - vi
    return((th, sqrt(1/S), sqrt(max((vt,0))), sqrt(max((v2,0))), sqrt(tau2), k, Gt))
}

// collector for pooled estimates
void pr_init()
{
    external string matrix PRS
    external real matrix PRM
    PRS = J(0,2,""); PRM = J(0,7,.)
}
// pool the terms of the dataset in memory (vars term b se team match)
void pool_spec(string scalar spec, string scalar termlist)
{
    external string matrix PRS
    external real matrix PRM
    string colvector T, terms
    real matrix D
    real scalar i
    T = st_sdata(., "term"); D = st_data(., ("b","se","team","match"))
    if (termlist == "") terms = uniqrows(T)
    else terms = tokens(termlist)'
    for (i = 1; i <= rows(terms); i++) {
        if (sum(T:==terms[i]) < 2) continue
        PRS = PRS \ (spec, terms[i])
        PRM = PRM \ re_pool(select(D, T:==terms[i]))
    }
}
// write the collector into the (empty) Stata dataset
void pr_dump()
{
    external string matrix PRS
    external real matrix PRM
    real scalar i
    string vector nm
    nm = ("pth","se_naive","se_team","se_twoway","tau","k","teams")
    st_addobs(rows(PRS))
    (void) st_addvar(("str40","str20"), ("spec","term"))
    (void) st_addvar("double", nm)
    st_sstore(., "spec", PRS[,1]); st_sstore(., "term", PRS[,2])
    for (i = 1; i <= 7; i++) st_store(., nm[i], PRM[,i])
}

// DL meta-regression with club-clustered variance:
//   - tau2: iterated DerSimonian-Laird estimate of the residual between-estimate variance, 50 passes
//     (stops before updating tau2 if |new-tau2|<1e-9; for several terms it does not converge but
//     cycles between two values);
//   - NOTE: after the loop the coefficients beta are those of the LAST PASS
//     (weights with the tau2 of the previous pass), while the weights w, the bread A and the clustered
//     variance use the final tau2.  Table S5 has this feature.
// Returns the last m coefficients and SEs.
void mreg(real colvector y, real colvector s, real matrix X, real colvector cl, real scalar m,
          real colvector bm, real colvector sem, real scalar tau2)
{
    real colvector w, res, beta
    real matrix A, XtW, U, meat, V
    real scalar n, k, it, Q, df, c, newt, G
    n = rows(X); k = cols(X); tau2 = 0
    for (it = 1; it <= 50; it++) {
        w = 1:/(s:^2 :+ tau2); XtW = (X:*w)'
        A = XtW*X; beta = lusolve(A, XtW*y); res = y - X*beta
        Q = sum(w:*res:^2); df = n - k
        c = sum(w) - trace(lusolve(A, cross(X, w:^2, X)))
        newt = max((0, (Q-df)/c))
        if (abs(newt - tau2) < 1e-9) break
        tau2 = newt
    }
    w = 1:/(s:^2 :+ tau2); XtW = (X:*w)'
    A = luinv(XtW*X); res = y - X*beta
    U = grpsum(X:*(w:*res), cl); G = rows(U)
    meat = cross(U, U)
    V = A*meat*A*(G/(G-1))
    bm = beta[|k-m+1 \ k|]; sem = sqrt(diagonal(V))[|k-m+1 \ k|]
}

// goal difference (own goals minus opponent goals) at time t, counting goals with time <= t.
// G = (match, team, t) sorted by match;  Q = (match, team, t)
real colvector goal_diff(real matrix G, real matrix Q)
{
    real colvector out
    real matrix info
    real scalar i, j, r1, r2
    transmorphic A
    info = panelsetup(G[,1], 1)
    A = asarray_create("real", 1)
    for (i = 1; i <= rows(info); i++) asarray(A, G[info[i,1],1], i)
    out = J(rows(Q), 1, 0)
    for (i = 1; i <= rows(Q); i++) {
        if (!asarray_contains(A, Q[i,1])) continue
        j = asarray(A, Q[i,1]); r1 = info[j,1]; r2 = info[j,2]
        out[i] = sum((G[|r1,3 \ r2,3|]:<=Q[i,3]) :* ((G[|r1,2 \ r2,2|]:==Q[i,2]) :- (G[|r1,2 \ r2,2|]:!=Q[i,2])))
    }
    return(out)
}
end

*==============================================================================
* 1. COMMON COVARIATES
*    team strength: points per game (ppg) and log squad value (lval), standardised
*    within league with the population SD (divisor n).
*==============================================================================
import delimited using "$B/team_cov.csv", clear varnames(1) case(preserve) asdouble stringcols(2)
gen byte lgn = (league=="England") + 2*(league=="Spain") + 3*(league=="Italy") + 4*(league=="Germany") + 5*(league=="France")
foreach v in ppg lval {
    bys lgn: egen double __m = mean(`v')
    gen double __d = (`v'-__m)^2
    bys lgn: egen double __s = total(__d)
    bys lgn: gen double __n = _N
    gen double z_`v' = (`v'-__m)/sqrt(__s/__n)
    drop __m __d __s __n
}
keep team league lgn ppg lval z_ppg z_lval
save "$TMP/team_cov.dta", replace

* match_cov: home indicator, z of own team (ppg, lval) and of the opponent
import delimited using "$B/match_cov.csv", clear varnames(1) case(preserve) asdouble
gen byte home = (side=="home")
merge m:1 team using "$TMP/team_cov.dta", keepusing(lgn z_ppg z_lval) nogenerate keep(master match)
rename (z_ppg z_lval) (z_team z_val)
rename team team_self
rename opp team
merge m:1 team using "$TMP/team_cov.dta", keepusing(z_ppg) nogenerate keep(master match)
rename team opp_team
rename team_self team
rename z_ppg z_opp
keep match team opp_team home lgn z_team z_val z_opp
sort match team
save "$TMP/match_cov.dta", replace

* goals (league match team_scoring t)
import delimited using "$B/goals.csv", clear varnames(1) case(preserve) asdouble
rename team_scoring team
keep if inlist(league,"England","Spain","Italy","Germany","France")
keep match team t
sort match t team
save "$TMP/goals.dta", replace

*==============================================================================
* 1b. SHARED INPUT FILES (pass records of the five leagues, player-on-pitch records)
*==============================================================================
* passes with coordinates are used by rate / observed / geometry / outcome;  Table 1 uses all rows
if $RUN_rate | $RUN_observed | $RUN_geometry | $RUN_outcome | $RUN_desc {
    import delimited using "$B/passes_all.csv", clear varnames(1) case(preserve) asdouble stringcols(1 6 10 15)
    keep if inlist(league,"England","Spain","Italy","Germany","France")
    drop period accurate subevent next_event
    compress
    save "$TMP/passes.dta", replace
}
* act_<league>.csv: players on the pitch (t_join, t_leave) and role
clear
tempfile actall
local i 0
foreach lg of global FIVE {
    import delimited using "$M/act_`lg'.csv", clear varnames(1) case(preserve) asdouble stringcols(4 12)
    keep match team player role t_join t_leave
    gen byte lgn = `i'+1
    local ++i
    if `i'>1 append using `actall'
    save `actall', replace
}
sort match team
save "$TMP/act_all.dta", replace
*==============================================================================
* 2. POOLING OF THE TEAM-MATCH ESTIMATES
*    Random-effects (DerSimonian-Laird) pooling, club-clustered variance eq. (3).
*    Feeds: Table 5 (J, JP), Table 7 (J JP R RP V U), Table 8 (variants prior_base/prior_simple),
*    Supplementary Tables S1 (X XP B BP), S2 (England recipient rules), S10 (JS).
*==============================================================================
* load est files (list of quoted paths): keep the successful fits used by the pooling scripts
program define load_est
    syntax anything(name=files)
    tempfile acc
    local i 0
    foreach f of local files {
        quietly import delimited using `"`f'"', clear varnames(1) case(preserve) asdouble stringcols(6)
        quietly keep if rc==0 & !missing(b) & !missing(se)
        quietly keep if se>0 & se<5 & abs(b)<10
        keep match team term b se
        local ++i
        quietly gen byte file = `i'
        if `i'>1 quietly append using `acc'
        quietly save `acc', replace
    }
    use `acc', clear
end

if $RUN_pool {
    di as txt _n "=== block pool: " c(current_time)
    mata: pr_init()
    foreach tag in J JP U V R RP X XP B BP JS {
        local fl ""
        foreach lg of global FIVE {
            local fl `"`fl' "$M/est_`lg'_`tag'.csv""'
        }
        load_est `fl'
        mata: pool_spec("`tag'", "")
    }
    * Table 8: five leagues, no tenure terms, earlier-match positions, all passes vs simple passes
    local t9 "inertia recip indeg trans same diff sim alter egoalterint"
    foreach v in base simple {
        local fl ""
        foreach lg of global FIVE {
            local fl `"`fl' "$V/prior_`v'_`lg'/est_`lg'_VP.csv""'
        }
        load_est `fl'
        mata: pool_spec("T8_`v'", "`t9'")
    }
    * Supplementary Table S2: England, recipient rules (baseline = next event)
    foreach v in base loose t5 acc {
        if "`v'"=="base" load_est "$V/prior_base_England/est_England_VP.csv"
        else             load_est "$V/prior_`v'/est_England_VP.csv"
        mata: pool_spec("S2_`v'", "")
    }
    clear
    mata: pr_dump()
    gen double z_team   = pth/se_team
    gen double z_twoway = pth/se_twoway
    gen double lo_team  = pth - 1.96*se_team
    gen double hi_team  = pth + 1.96*se_team
    gen double lo_twoway = pth - 1.96*se_twoway
    gen double hi_twoway = pth + 1.96*se_twoway
    save "$TMP/pooled.dta", replace

    * ---- paper tables (one CSV each) ----
    foreach x in "Table5_pooled_choice_effects|J JP" "Table8_robustness_specs|J JP R RP V U" ///
                 "TableS1_tenure_specs|X XP B BP" "TableS10_strict_fallback|JS" ///
                 "Table9_recipient_pass_type|T8_base T8_simple" "TableS2_recipient_rules|S2_base S2_loose S2_t5 S2_acc" {
        local fn = substr("`x'",1,strpos("`x'","|")-1)
        local sp = substr("`x'",strpos("`x'","|")+1,.)
        use "$TMP/pooled.dta", clear
        gen byte keep_ = 0
        foreach s of local sp {
            replace keep_ = 1 if spec=="`s'"
        }
        keep if keep_
        drop keep_
        order spec term pth se_team z_team lo_team hi_team se_twoway z_twoway lo_twoway hi_twoway se_naive tau k teams
        export delimited using "$OUT/`fn'.csv", replace
    }
}
* NOTE: section 4 (pass frequency) is placed before section 3 on purpose: the Benjamini-Hochberg count in section 3 uses its results.
*==============================================================================
* 4. PASS FREQUENCY: PIECEWISE POISSON RATE MODELS WITH PLAYER-TIME EXPOSURE
*    Interval = one team-match with constant set of players on the pitch, constant score state and
*    constant clock band.  Poisson with exposure, club-clustered errors.
*    Feeds: Table 3 (venue / strength), Table 4 (score and clock), Supp. Tables S3, S4.
*==============================================================================
mata:
// build the intervals.  ACT=(match,team,t_join,t_leave) sorted by match,team;
// GL=(match,team,t) sorted by match;  PS=(match,team,t,has_receiver) sorted by match,team,t
real matrix make_intervals(real matrix ACT, real matrix GL, real matrix PS)
{
    real scalar LATE, i, j, a, b, n, tend, diff, g1, g2, p1, p2, c
    real matrix actinfo, glinfo, psinfo, out, A, glm, psm
    real colvector cuts, gt, gteam, pt, phr, r
    real scalar m, t, nall, nrec, band
    transmorphic GM, PM
    LATE = 3600
    actinfo = panelsetup(ACT[,1]:*1e7 :+ ACT[,2], 1); glinfo = panelsetup(GL[,1], 1); psinfo = panelsetup(PS[,1]:*1e7 :+ PS[,2], 1)
    GM = asarray_create("real", 1); PM = asarray_create("real", 1)
    for (i = 1; i <= rows(glinfo); i++) asarray(GM, GL[glinfo[i,1],1], i)
    for (i = 1; i <= rows(psinfo); i++) asarray(PM, PS[psinfo[i,1],1]*1e7 + PS[psinfo[i,1],2], i)
    out = J(0, 10, .)
    for (i = 1; i <= rows(actinfo); i++) {
        A = ACT[|actinfo[i,1],1 \ actinfo[i,2],4|]
        m = A[1,1]; t = A[1,2]
        if (!asarray_contains(PM, m*1e7 + t)) continue
        j = asarray(PM, m*1e7 + t); p1 = psinfo[j,1]; p2 = psinfo[j,2]
        pt = PS[|p1,3 \ p2,3|]; phr = PS[|p1,4 \ p2,4|]
        tend = max(A[,4])
        if (asarray_contains(GM, m)) {
            j = asarray(GM, m); g1 = glinfo[j,1]; g2 = glinfo[j,2]
            gt = GL[|g1,3 \ g2,3|]; gteam = GL[|g1,2 \ g2,2|]
            r = (gt:>0) :& (gt:<tend); gt = select(gt, r); gteam = select(gteam, r)
        }
        else {
            gt = J(0,1,.); gteam = J(0,1,.)
        }
        cuts = (0 \ tend \ gt \ select(A[,3], (A[,3]:>0):&(A[,3]:<tend)) \ select(A[,4], (A[,4]:>0):&(A[,4]:<tend)))
        if (LATE < tend) cuts = cuts \ LATE
        cuts = uniqrows(cuts)
        for (c = 1; c < rows(cuts); c++) {
            a = cuts[c]; b = cuts[c+1]
            if (b - a < 1e-6) continue
            n = sum((A[,3]:<=a+1e-6) :& (A[,4]:>=b-1e-6))
            if (n == 0) continue
            diff = sum((gt:<=a+1e-6) :& (gteam:==t)) - sum((gt:<=a+1e-6) :& (gteam:!=t))
            r = (pt:>=a) :& (pt:<b)
            nall = sum(r); nrec = sum(r:*phr)
            band = min((floor(a/900), 6))
            out = out \ (m, t, (diff>0), (diff<0), (a >= LATE-1e-6), band, n*(b-a), nall, nrec, 0)
        }
    }
    return(out)
}
end

* the cluster-robust variance is multiplied by the small-sample factor (N-1)/(N-k) on top of G/(G-1),
* k = number of columns of the design matrix (without omitted base levels)
program define pysf, rclass
    local k 0
    local cn : colnames e(b)
    foreach c of local cn {
        if !regexm("`c'","(^|[0-9])b\.|o\.") local ++k
    }
    return scalar sf = sqrt((e(N)-1)/(e(N)-`k'))
    return scalar k = `k'
end

* post helper: coefficients of the last poisson fit, scaled by the small-sample factor
program define post_coefs
    * args: handle outcome model sf  then pairs  var=label ...
    gettoken ph 0 : 0
    gettoken oc 0 : 0
    gettoken md 0 : 0
    gettoken sf 0 : 0
    foreach pr of local 0 {
        local v = substr("`pr'",1,strpos("`pr'","=")-1)
        local l = substr("`pr'",strpos("`pr'","=")+1,.)
        capture local bb = _b[`v']
        if _rc continue
        post `ph' ("`oc'") ("`md'") ("`l'") (_b[`v']) (_se[`v']*`sf')
    }
end
program define post_lin
    * args: handle outcome model sf label expr
    args ph oc md sf label expr
    quietly lincom `expr'
    post `ph' ("`oc'") ("`md'") ("`label'") (r(estimate)) (r(se)*`sf')
end

if $RUN_rate {
    di as txt _n "=== block rate: " c(current_time)
    * ---- intervals ----
    use "$TMP/passes.dta", clear
    keep if !missing(x0)
    gen byte hasrec = !missing(receiver)
    keep match team t hasrec
    sort match team t
    tempfile ps
    save `ps'
    use "$TMP/goals.dta", clear
    sort match
    tempfile gl
    save `gl'
    use "$TMP/act_all.dta", clear
    sort match team
    mata: ACT = st_data(., ("match","team","t_join","t_leave"))
    use `gl', clear
    mata: GL = st_data(., ("match","team","t"))
    use `ps', clear
    mata: PS = st_data(., ("match","team","t","hasrec"))
    clear
    mata: IV = make_intervals(ACT, GL, PS)
    mata: st_addobs(rows(IV)); (void) st_addvar("double", ("match","team","lead","trail","late","band","expo","n_all","n_rec","dummy")); st_store(., ("match","team","lead","trail","late","band","expo","n_all","n_rec","dummy"), IV)
    drop dummy
    mata: mata drop ACT GL PS IV
    * attach covariates (rows whose team-match is not in match_cov are dropped)
    merge m:1 match team using "$TMP/match_cov.dta", keep(match) nogenerate
    gen double logexpo = ln(expo)
    gen double ll = lead*late
    gen double tl = trail*late
    rename (team z_opp z_team z_val) (club opp teamz val)
    rename teamz team
    rename lgn lg
    save "$TMP/rate_intervals.dta", replace
    quietly summarize n_all
    local att = r(sum)
    quietly summarize n_rec
    local rec = r(sum)
    local nint = _N
    preserve
    keep match club
    duplicates drop
    local ntm = _N
    restore
    di as txt "intervals `nint'  team-matches `ntm'  attempts `att'  with receiver `rec'"
    preserve
    clear
    set obs 1
    gen n_intervals = `nint'
    gen n_team_matches = `ntm'
    gen attempts = `att'
    gen receiver_passes = `rec'
    export delimited using "$OUT/rate_poisson_counts.csv", replace
    restore

    * ---- models ----
    tempfile rr
    tempname ph
    postfile `ph' str24 outcome str16 model str24 term double(b se) using `rr', replace
    foreach oc in n_all n_rec {
        local lab = cond("`oc'"=="n_all","all attempts","passes with receiver")
        * A: venue, club FE
        quietly poisson `oc' home opp i.club, exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "`lab'" A `sf' home=home opp=opp
        * A2: strength, league FE
        quietly poisson `oc' home opp team i.lg, exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "`lab'" A2 `sf' home=home opp=opp team=team
        * A3: squad value instead of team strength;  A4 both
        quietly poisson `oc' home opp val i.lg, exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "`lab'" A3 `sf' home=home opp=opp val=val
        quietly poisson `oc' home opp team val i.lg, exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "`lab'" A4 `sf' home=home opp=opp team=team val=val
        * B: score state x late, clock bands, club FE
        quietly poisson `oc' home opp lead trail ll tl ib0.band i.club, exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "`lab'" B `sf' home=home opp=opp lead=lead trail=trail ll=lead:late tl=trail:late
        post_lin `ph' "`lab'" B `sf' lead_late "lead + ll"
        post_lin `ph' "`lab'" B `sf' trail_late "trail + tl"
        post_lin `ph' "`lab'" B `sf' trail_vs_lead "trail - lead"
        post_lin `ph' "`lab'" B `sf' gap_growth "tl - ll"
        post_lin `ph' "`lab'" B `sf' trail_vs_lead_late "trail + tl - lead - ll"
        forvalues j = 1/6 {
            post_coefs `ph' "`lab'" B `sf' `j'.band=band`j'
        }
        * B without score state
        quietly poisson `oc' home opp ib0.band i.club, exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "`lab'" B_noscore `sf' home=home opp=opp
        forvalues j = 1/6 {
            post_coefs `ph' "`lab'" B_noscore `sf' `j'.band=band`j'
        }
    }
    * by league (all attempts), venue model A
    forvalues l = 1/5 {
        local lgname : word `l' of $FIVE
        quietly poisson n_all home opp i.club if lg==`l', exposure(expo) vce(cluster club)
        pysf
        local sf = r(sf)
        post_coefs `ph' "all attempts" "by_`lgname'" `sf' home=home opp=opp
    }
    postclose `ph'
    use `rr', clear
    gen double rr = exp(b)
    gen double lo = exp(b - 1.96*se)
    gen double hi = exp(b + 1.96*se)
    gen double z = b/se
    save "$TMP/rate_res.dta", replace
    export delimited using "$OUT/rate_poisson_all.csv", replace
    foreach x in "Table3_frequency_venue_strength|A A2" "Table4_frequency_score_clock|B B_noscore" ///
                 "TableS3_frequency_squad_value|A3 A4" "TableS4_frequency_by_league|by_England by_Spain by_Italy by_Germany by_France" {
        local fn = substr("`x'",1,strpos("`x'","|")-1)
        local sp = substr("`x'",strpos("`x'","|")+1,.)
        use "$TMP/rate_res.dta", clear
        gen byte keep_ = 0
        foreach s of local sp {
            replace keep_ = 1 if model=="`s'"
        }
        keep if keep_
        drop keep_
        export delimited using "$OUT/`fn'.csv", replace
    }
}
*==============================================================================
* 3. SECOND-STAGE (META-)REGRESSIONS ON CONTEXT
*    Team-match estimates regressed on team strength (ppg), squad value, home, opponent strength,
*    DL random-effects weights, league fixed effects, club-clustered variance.
*    Feeds: Supplementary Table S5 (tag JP), Benjamini-Hochberg count, club-FE meta-regression.
*==============================================================================
mata:
// one second-stage regression for one term.  data in memory: term b se team lgn z_team z_val home z_opp
// cols: which of (z_team, z_val, home, z_opp) enter; fe: "league" or "club".
// returns (n, tau, b_1..b_m, se_1..se_m)
real rowvector mod_run(string scalar term, real rowvector cc, string scalar fe)
{
    string colvector T
    real matrix D, M, X
    real colvector y, s, tm, lg, bm, sem, teams
    real scalar tau2, m
    T = st_sdata(., "term")
    D = select(st_data(., ("b","se","team","lgn","z_team","z_val","home","z_opp")), T:==term)
    y = D[,1]; s = D[,2]; tm = D[,3]; lg = D[,4]
    M = D[|1,5 \ .,8|]; M = M[, cc]; m = cols(cc)
    if (fe == "league") X = (J(rows(y),1,1), (lg:==2), (lg:==3), (lg:==4), (lg:==5), M)
    else {
        teams = uniqrows(tm)
        X = (J(rows(y),1,1), ((tm*J(1,rows(teams)-1,1)) :== teams[|2,1 \ rows(teams),1|]'), M)
    }
    mreg(y, s, X, tm, m, bm, sem, tau2)
    return((rows(y), sqrt(tau2), bm', sem'))
}
end

* team-match estimates of a tag with the covariates attached
program define load_mod
    args tag
    local fl ""
    foreach lg of global FIVE {
        local fl `"`fl' "$M/est_`lg'_`tag'.csv""'
    }
    load_est `fl'
    rename file lgn
    quietly merge m:1 match team using "$TMP/match_cov.dta", keep(match) nogenerate
    quietly drop if missing(z_team) | missing(z_opp)
    * (rows whose team is not in team_cov or whose opponent has no strength are not used)
end

if $RUN_moderators {
    di as txt _n "=== block moderators: " c(current_time)
    tempfile modres
    tempname ph
    postfile `ph' str6 tag str8 model str12 term str10 moderator double(b se n tau) using `modres', replace
    local terms "inertia recip indeg trans same diff sim alter egoalterint tie_ten1 tie_sm270 tie_sm900"
    foreach tag in JP J {
        load_mod `tag'
        foreach spec in "strength|1 3 4|strength home opponent" "value|2 3 4|value home opponent" ///
                        "both|1 2 3 4|strength value home opponent" {
            tokenize "`spec'", parse("|")
            local model "`1'"
            local cl "`3'"
            local nm "`5'"
            foreach t of local terms {
                local cols = subinstr("`cl'"," ",",",.)
                mata: st_matrix("MR", mod_run("`t'", (`cols'), "league"))
                local j 0
                foreach mname of local nm {
                    local ++j
                    local m = colsof(MR)/2 - 1
                    post `ph' ("`tag'") ("`model'") ("`t'") ("`mname'") (MR[1,2+`j']) (MR[1,2+`m'+`j']) (MR[1,1]) (MR[1,2])
                }
            }
        }
    }
    postclose `ph'
    * ---- club fixed effects (home advantage within club), tag JP
    tempfile fe
    tempname pf
    postfile `pf' str12 term double(n tau b_home se_home z_home b_opp se_opp z_opp) using `fe', replace
    load_mod JP
    foreach t in inertia recip indeg trans same diff sim alter egoalterint tie_ten1 tie_sm270 {
        mata: st_matrix("MR", mod_run("`t'", (3,4), "club"))
        post `pf' ("`t'") (MR[1,1]) (MR[1,2]) (MR[1,3]) (MR[1,5]) (MR[1,3]/MR[1,5]) (MR[1,4]) (MR[1,6]) (MR[1,4]/MR[1,6])
    }
    postclose `pf'
    use `fe', clear
    export delimited using "$OUT/TableS5b_metareg_club_fe.csv", replace

    use `modres', clear
    gen double z = b/se
    save "$TMP/modres.dta", replace
    export delimited using "$OUT/TableS5_moderators.csv", replace

    * ---- Benjamini-Hochberg over all tests (strength model, tag JP, plus 3 frequency tests)
    keep if tag=="JP" & model=="strength"
    gen str30 test = term + " | " + moderator
    keep test z
    capture confirm file "$TMP/rate_res.dta"
    if _rc==0 {
        preserve
        use "$TMP/rate_res.dta", clear
        keep if outcome=="all attempts" & model=="A2" & inlist(term,"team","home","opp")
        gen str30 test = "frequency | " + cond(term=="team","strength",cond(term=="opp","opponent","home"))
        keep test z
        tempfile fz
        save `fz'
        restore
        append using `fz'
    }
    else di as err "rate_res.dta missing: run block rate first (BH count without frequency tests)"
    gen double p = 2*normal(-abs(z))
    gsort p
    gen int rank = _n
    local m = _N
    gen byte bh_ok = p <= 0.05*rank/`m'
    quietly summarize rank if bh_ok
    local thr = r(max)
    gen byte bh_signif = rank <= `thr'
    quietly count if p <= 0.05/`m'
    local bonf = r(N)
    di as txt "Benjamini-Hochberg (q=0.05): tests `m'  significant `thr'  Bonferroni `bonf'"
    gen int n_tests = `m'
    gen int n_bh_significant = `thr'
    gen int n_bonferroni = `bonf'
    export delimited using "$OUT/TableS5c_benjamini_hochberg.csv", replace
}
*==============================================================================
* 2b. DESCRIPTIVE STATISTICS
*     Feeds: Table 1 (sample by league) and Table 2 (descriptive statistics).
*==============================================================================
* posts N, mean, population sd, min, max of a variable
program define stat_row
    args ph label var
    quietly summarize `var'
    local n = r(N)
    local sdp = sqrt(r(Var)*(r(N)-1)/r(N))
    if r(N)==1 local sdp = 0
    post `ph' ("`label'") (`n') (r(mean)) (`sdp') (r(min)) (r(max))
end

if $RUN_desc {
    di as txt _n "=== block desc: " c(current_time)
    * FIFA rating of the player (ov18, else ov19, else ov17)
    import delimited using "$B/fifa_overall.csv", clear varnames(1) case(preserve) asdouble
    gen double ovr = cond(!missing(ov18), ov18, cond(!missing(ov19), ov19, ov17))
    keep if !missing(ovr)
    rename wy_id player
    keep player ovr
    duplicates drop player, force
    tempfile ov
    save `ov'
    import delimited using "$B/nodes.csv", clear varnames(1) case(preserve) asdouble stringcols(1 4 6 9 14)
    merge m:1 player using `ov', keep(master match) nogenerate
    gen byte rated = !missing(ovr)
    gen double tin = cond(missing(t_in), 0, t_in)
    gen double mins = max(0, (t_out - tin)/60)
    save "$TMP/nodes.dta", replace

    * ---- Table 1 ----
    use "$TMP/passes.dta", clear
    gen byte hasrec = !missing(receiver)
    gen double one = 1
    collapse (count) passes=one (sum) rec=hasrec, by(league)
    tempfile t1a
    save `t1a'
    use "$TMP/passes.dta", clear
    keep league match team
    duplicates drop
    bys league match: gen byte m1 = _n==1
    bys league team: gen byte c1 = _n==1
    collapse (sum) matches=m1 clubs=c1, by(league)
    tempfile t1b
    save `t1b'
    use "$TMP/nodes.dta", clear
    keep league match team player rated
    gen double one = 1
    collapse (count) nodes=one (sum) nrated=rated, by(league)
    tempfile t1c
    save `t1c'
    use "$TMP/nodes.dta", clear
    keep league match team
    duplicates drop
    gen double one = 1
    collapse (count) tmatches=one, by(league)
    tempfile t1d
    save `t1d'
    use "$TMP/nodes.dta", clear
    keep league player
    duplicates drop
    gen double one = 1
    collapse (count) players=one, by(league)
    tempfile t1e
    save `t1e'
    import delimited using "$B/goals.csv", clear varnames(1) case(preserve) asdouble stringcols(1)
    keep if inlist(league,"England","Spain","Italy","Germany","France")
    gen double one = 1
    collapse (count) goals=one, by(league)
    merge 1:1 league using `t1a', nogenerate
    merge 1:1 league using `t1b', nogenerate
    merge 1:1 league using `t1c', nogenerate
    merge 1:1 league using `t1d', nogenerate
    merge 1:1 league using `t1e', nogenerate
    gen byte ord = (league=="England") + 2*(league=="Spain") + 3*(league=="Italy") + 4*(league=="Germany") + 5*(league=="France")
    sort ord
    * all-league unique players
    preserve
    use "$TMP/nodes.dta", clear
    keep player
    duplicates drop
    local nplayers_all = _N
    restore
    set obs 6
    replace league = "Total" in 6
    foreach v in goals passes rec matches clubs nodes nrated tmatches players {
        quietly summarize `v' if ord<6
        replace `v' = r(sum) in 6
    }
    replace players = `nplayers_all' in 6
    gen double recipient_pct = 100*rec/passes
    gen double goals_per_match = goals/matches
    gen double rated_pct = 100*nrated/nodes
    keep league matches clubs tmatches passes recipient_pct players goals_per_match rated_pct
    rename tmatches team_matches
    export delimited using "$OUT/Table1_sample_by_league.csv", replace

    * ---- Table 2 ----
    tempfile t2
    tempname ph
    postfile `ph' str60 variable double(n mean sd min max) using `t2', replace
    use "$TMP/passes.dta", clear
    gen byte hasrec = !missing(receiver)
    gen double one = 1
    collapse (count) passes=one (sum) rec=hasrec, by(league match team)
    stat_row `ph' "Passes per team-match" passes
    stat_row `ph' "Passes with a recipient" rec
    use "$TMP/nodes.dta", clear
    gen double one = 1
    preserve
    collapse (count) players=one, by(league match team)
    stat_row `ph' "Players appearing in the team-match" players
    restore
    stat_row `ph' "Minutes on the pitch" mins
    stat_row `ph' "Started the match (share)" started
    stat_row `ph' "Age (years)" age
    stat_row `ph' "FIFA overall rating" ovr
    import delimited using "$B/team_cov.csv", clear varnames(1) case(preserve) asdouble
    stat_row `ph' "Points per game" ppg
    stat_row `ph' "Squad value (log euro)" lval
    * pairs of teammates
    import delimited using "$B/dyads.csv", clear varnames(1) case(preserve) asdouble stringcols(6)
    preserve
    keep if shared_src=="tm_exact" & !missing(shared_tenure_days)
    stat_row `ph' "Shared tenure, exact dates, days" shared_tenure_days
    restore
    preserve
    keep if !missing(shared_tenure_days)
    gen double ten1 = shared_tenure_days>=365
    stat_row `ph' "Tenure of 1 year or more, share" ten1
    restore
    * shared minutes before the match (cumulative overlap in earlier matches of the same club)
    keep match team i j
    gen double pa = min(i,j)
    gen double pb = max(i,j)
    tempfile dy
    save `dy'
    use "$TMP/nodes.dta", clear
    keep league match team date player tin t_out
    rename (player tin t_out) (pa tin_a tout_a)
    tempfile na
    save `na'
    rename (pa tin_a tout_a) (pb tin_b tout_b)
    joinby league match team using `na'
    keep if pa < pb
    gen double ov = max(0, (min(tout_a, tout_b) - max(tin_a, tin_b))/60)
    sort team pa pb date match
    by team pa pb: gen double cum = sum(ov) - ov
    keep match team pa pb cum
    merge 1:m match team pa pb using `dy', keep(using match) nogenerate
    replace cum = 0 if missing(cum)
    stat_row `ph' "Shared minutes before the match" cum
    gen double zero = cum==0
    stat_row `ph' "Shared minutes = 0 (share)" zero
    postclose `ph'
    use `t2', clear
    export delimited using "$OUT/Table2_descriptives.csv", replace
}
*==============================================================================
* 5. OBSERVED PASSING OUTCOMES IN FOOTBALL UNITS
*    Per team-match: mean pass length, mean forward displacement, share of passes with an inferred
*    receiver, share of same-role passes, share of return passes in chains; OLS with club fixed
*    effects (home, opponent strength) and with league fixed effects (+ team strength), club-clustered.
*    Feeds: Table 6, Figure 2 (panels B to F).
*==============================================================================
if $RUN_observed {
    di as txt _n "=== block observed: " c(current_time)
    use "$TMP/passes.dta", clear
    keep if !missing(x0) & !missing(x1)
    * role of passer and receiver in this team-match (act_<league>.csv)
    tempfile roles
    preserve
    use "$TMP/act_all.dta", clear
    keep match team player role
    save `roles'
    restore
    rename passer player
    merge m:1 match team player using `roles', keepusing(role) keep(master match) nogenerate
    rename (player role) (passer role_p)
    rename receiver player
    merge m:1 match team player using `roles', keepusing(role) keep(master match) nogenerate
    rename (player role) (receiver role_r)
    gen double length = sqrt((x1-x0)^2 + (y1-y0)^2)
    gen double fwd = x1 - x0
    sort match team seq
    by match team: gen double prev_passer   = passer[_n-1]
    by match team: gen double prev_receiver = receiver[_n-1]
    by match team: gen double prev_cont     = next_team_same[_n-1]
    gen byte has = !missing(receiver)
    gen double ident = has
    gen double same = (role_p==role_r & role_p!="") if has
    gen byte chain = (prev_receiver==passer & prev_cont==1 & !missing(prev_receiver))
    gen double ret = (has & receiver==prev_passer) if chain
    collapse (mean) length fwd same ret ident (count) n=seq, by(match team)
    keep if n>=30
    foreach v in same ret ident {
        replace `v' = 100*`v'
    }
    merge 1:1 match team using "$TMP/match_cov.dta", keep(master match) nogenerate keepusing(home z_opp lgn)
    replace home = 0 if missing(home)
    merge m:1 team using "$TMP/team_cov.dta", keep(master match) nogenerate keepusing(z_ppg)
    rename (z_opp z_ppg) (opp team_z)
    drop if missing(opp) | missing(same)
    local ntm = _N
    * baseline (mean over team-matches)
    tempfile obs
    tempname ph
    postfile `ph' str8 outcome str12 term double(b se) using `obs', replace
    foreach y in length fwd same ret ident {
        quietly summarize `y'
        post `ph' ("`y'") ("baseline") (r(mean)) (.)
        quietly regress `y' home opp i.team, vce(cluster team)
        post `ph' ("`y'") ("home") (_b[home]) (_se[home])
        post `ph' ("`y'") ("opp") (_b[opp]) (_se[opp])
        quietly regress `y' home opp team_z i.lgn, vce(cluster team)
        post `ph' ("`y'") ("team") (_b[team_z]) (_se[team_z])
    }
    postclose `ph'
    use `obs', clear
    gen double z = b/se
    gen double lo = b - 1.96*se
    gen double hi = b + 1.96*se
    gen int n_team_matches = `ntm'
    save "$TMP/obs_res.dta", replace
    export delimited using "$OUT/Table6_observed_shares.csv", replace
}
*==============================================================================
* 6. DIRECT PASS GEOMETRY BY VENUE, OPPONENT, SCORE STATE AND CLOCK
*    Pass-level OLS (pass length, forward displacement, 0-100 pitch scale) with club fixed effects
*    (within transformation), club-clustered SEs with the factor G/(G-1) only.
*    Feeds: Supplementary Table S6.
*==============================================================================
mata:
// YX = (y1, y2, X(12 columns)), cl = club id (any real).  Within-club OLS, clustered variance G/(G-1).
// returns for each y a row: b(12), se(12), b_lead_late, se_lead_late, b_trail_late, se_trail_late
real matrix geom_fit(real matrix YX, real colvector cl)
{
    real colvector p, ks, e
    real matrix info, W, Xw, XtXi, S, V, out, sc
    real scalar i, g, G, K, iy
    real colvector yw, b
    p = order(cl, 1); W = YX[p, .]; ks = cl[p]
    info = panelsetup(ks, 1); G = rows(info); K = cols(W)
    // within transformation of all columns
    for (g = 1; g <= G; g++) {
        W[|info[g,1],1 \ info[g,2],K|] = W[|info[g,1],1 \ info[g,2],K|] :- mean(W[|info[g,1],1 \ info[g,2],K|])
    }
    Xw = W[|1,3 \ .,K|]
    XtXi = invsym(cross(Xw, Xw))
    out = J(2, 28, .)
    for (iy = 1; iy <= 2; iy++) {
        yw = W[,iy]
        b = XtXi*cross(Xw, yw); e = yw - Xw*b
        S = J(12, 12, 0)
        for (g = 1; g <= G; g++) {
            sc = cross(Xw[|info[g,1],1 \ info[g,2],12|], e[|info[g,1] \ info[g,2]|])
            S = S + sc*sc'
        }
        V = (G/(G-1))*XtXi*S*XtXi
        out[iy, 1..12] = b'
        out[iy, 13..24] = sqrt(diagonal(V))'
        out[iy, 25] = b[3] + b[5]; out[iy, 26] = sqrt(V[3,3] + V[5,5] + 2*V[3,5])
        out[iy, 27] = b[4] + b[6]; out[iy, 28] = sqrt(V[4,4] + V[6,6] + 2*V[4,6])
    }
    return(out)
}
end

if $RUN_geometry {
    di as txt _n "=== block geometry: " c(current_time)
    use "$TMP/passes.dta", clear
    keep if !missing(x0, y0, x1, y1)
    merge m:1 match team using "$TMP/match_cov.dta", keep(match) nogenerate keepusing(home z_opp)
    drop if missing(z_opp)
    keep match team t x0 y0 x1 y1 home z_opp
    gen double length = sqrt((x1-x0)^2 + (y1-y0)^2)
    gen double fwd = x1 - x0
    preserve
    use "$TMP/goals.dta", clear
    sort match t team
    mata: GL = st_data(., ("match","team","t"))
    restore
    gen double diff = .
    mata: st_store(., "diff", goal_diff(GL, st_data(., ("match","team","t"))))
    mata: mata drop GL
    gen double lead  = diff>0
    gen double trail = diff<0
    gen double late  = t>=3600
    gen double ll = lead*late
    gen double tl = trail*late
    gen byte band = min(floor(t/900),6)
    forvalues j = 1/6 {
        gen double b`j' = band==`j'
    }
    local npass = _N
    quietly summarize length
    local m_length = r(mean)
    quietly summarize fwd
    local m_fwd = r(mean)
    mata: YX = st_data(., ("length","fwd","home","z_opp","lead","trail","ll","tl","b1","b2","b3","b4","b5","b6"))
    mata: CL = st_data(., "team")
    mata: st_matrix("GR", geom_fit(YX, CL))
    mata: mata drop YX CL
    clear
    local names "home opp lead trail ll tl b1 b2 b3 b4 b5 b6"
    set obs 28
    gen str12 outcome = ""
    gen str12 term = ""
    gen double b = .
    gen double se = .
    local r 0
    forvalues iy = 1/2 {
        local yn : word `iy' of length fwd
        local j 0
        foreach nm of local names {
            local ++j
            local ++r
            replace outcome = "`yn'" in `r'
            replace term = "`nm'" in `r'
            replace b = GR[`iy',`j'] in `r'
            replace se = GR[`iy',`=12+`j''] in `r'
        }
        local ++r
        replace outcome = "`yn'" in `r'
        replace term = "lead_late" in `r'
        replace b = GR[`iy',25] in `r'
        replace se = GR[`iy',26] in `r'
        local ++r
        replace outcome = "`yn'" in `r'
        replace term = "trail_late" in `r'
        replace b = GR[`iy',27] in `r'
        replace se = GR[`iy',28] in `r'
    }
    keep if outcome!=""
    gen double z = b/se
    gen double lo = b - 1.96*se
    gen double hi = b + 1.96*se
    gen long n_passes = `npass'
    gen double mean_y = cond(outcome=="length", `m_length', `m_fwd')
    save "$TMP/geom_res.dta", replace
    export delimited using "$OUT/TableS6_pass_geometry.csv", replace
}
*==============================================================================
* 7. TEMPO AND GOAL DIFFERENCE
*    Tempo = log(pass attempts / (players on the pitch x seconds)) in a window, standardised.
*    Outcome: final goal difference (official score from raw/matches/matches_<league>.json);
*    club fixed effects, home, opponent strength; club-clustered SEs.
*    Feeds: Supplementary Table S7.
*==============================================================================
mata:
// total player-time and number of passes in window [a,b) of one team-match
void win_np(real matrix A, real colvector pt, real scalar a, real scalar b, real scalar n, real scalar e)
{
    real colvector ts
    real scalar i, u, v
    ts = (a \ b \ select(A[,3], (A[,3]:>a):&(A[,3]:<b)) \ select(A[,4], (A[,4]:>a):&(A[,4]:<b)))
    ts = uniqrows(ts)
    e = 0
    for (i = 1; i < rows(ts); i++) {
        u = ts[i]; v = ts[i+1]
        e = e + (v-u)*sum((A[,3]:<=u+1e-6) :& (A[,4]:>=v-1e-6))
    }
    n = sum((pt:>=a) :& (pt:<b))
}
// ACT=(match,team,t_join,t_leave) sorted by match,team; GL=(match,team,t) sorted by match;
// PS=(match,team,t) sorted by match,team.   Returns (match,team,n_f,e_f,n_l,e_l,n_w,e_w,goals_before_LM)
real matrix make_windows(real matrix ACT, real matrix GL, real matrix PS)
{
    real scalar LM, i, j, tend, bf, g1, g2, p1, p2, m, t, n, e
    real matrix actinfo, glinfo, psinfo, out, A
    real colvector gt, pt, pos
    transmorphic GM, PM
    LM = 1200
    actinfo = panelsetup(ACT[,1]:*1e7 :+ ACT[,2], 1); glinfo = panelsetup(GL[,1], 1); psinfo = panelsetup(PS[,1]:*1e7 :+ PS[,2], 1)
    GM = asarray_create("real", 1); PM = asarray_create("real", 1)
    for (i = 1; i <= rows(glinfo); i++) asarray(GM, GL[glinfo[i,1],1], i)
    for (i = 1; i <= rows(psinfo); i++) asarray(PM, PS[psinfo[i,1],1]*1e7 + PS[psinfo[i,1],2], i)
    out = J(0, 9, .)
    for (i = 1; i <= rows(actinfo); i++) {
        A = ACT[|actinfo[i,1],1 \ actinfo[i,2],4|]
        m = A[1,1]; t = A[1,2]
        if (!asarray_contains(PM, m*1e7 + t)) continue
        j = asarray(PM, m*1e7 + t); p1 = psinfo[j,1]; p2 = psinfo[j,2]
        pt = PS[|p1,3 \ p2,3|]
        tend = max(A[,4])
        if (asarray_contains(GM, m)) {
            j = asarray(GM, m); gt = GL[|glinfo[j,1],3 \ glinfo[j,2],3|]
        }
        else {
            gt = J(0,1,.)
        }
        pos = select(gt, gt:>0)
        bf = min((pos \ tend))
        out = out \ (m, t, ., ., ., ., ., ., sum(gt:<LM))
        win_np(A, pt, 0, bf, n, e);   out[rows(out),3] = n; out[rows(out),4] = e
        win_np(A, pt, 0, LM, n, e);   out[rows(out),5] = n; out[rows(out),6] = e
        win_np(A, pt, 0, tend, n, e); out[rows(out),7] = n; out[rows(out),8] = e
    }
    return(out)
}
end

* final scores from the Wyscout match files: one line per "side", "wyId" marker (filefilter), then read
program define read_scores
    args lg
    tempfile f1 f2
    quietly filefilter "$RAW/matches/matches_`lg'.json" `f1', from("\Qside\Q: \Q") to("\n@S ") replace
    quietly filefilter `f1' `f2', from("\QwyId\Q: ") to("\n@W ") replace
    tempname fh
    file open `fh' using `f2', read text
    tempname ph
    tempfile out
    postfile `ph' double(match team score) str5 side using `out', replace
    local nk 0
    file read `fh' line
    while r(eof)==0 {
        if substr(`"`macval(line)'"',1,3)=="@S " {
            local ln = substr(`"`macval(line)'"',1,80)
            if regexm(`"`ln'"', `"^@S (home|away)", "teamId": ([0-9]+), "score": ([0-9]+)"') {
                local nk = `nk' + 1
                local s`nk'_side = regexs(1)
                local s`nk'_team = regexs(2)
                local s`nk'_score = regexs(3)
            }
        }
        else if substr(`"`macval(line)'"',1,3)=="@W " {
            local ln = substr(`"`macval(line)'"',1,30)
            if regexm(`"`ln'"', "^@W ([0-9]+)") {
                local mid = regexs(1)
                forvalues q = 1/`nk' {
                    post `ph' (`mid') (`s`q'_team') (`s`q'_score') ("`s`q'_side'")
                }
                local nk 0
            }
        }
        file read `fh' line
    }
    file close `fh'
    postclose `ph'
    use `out', clear
end

if $RUN_outcome {
    di as txt _n "=== block outcome: " c(current_time)
    * ---- official scores, goal difference, opponent ----
    clear
    tempfile sc
    local i 0
    foreach lg of global FIVE {
        read_scores `lg'
        local ++i
        if `i'>1 append using `sc'
        save `sc', replace
    }
    use `sc', clear
    bys match: gen double score_o = score[3-_n]
    bys match: gen double opp_team = team[3-_n]
    gen double gd = score - score_o
    gen byte home = side=="home"
    keep match team gd home opp_team
    sort match team
    save "$TMP/scores.dta", replace
    di as txt "scores: " _N " team-matches"

    * ---- windows ----
    use "$TMP/passes.dta", clear
    keep if !missing(x0)
    keep match team t
    sort match team t
    mata: PS = st_data(., ("match","team","t"))
    use "$TMP/goals.dta", clear
    sort match t team
    mata: GL = st_data(., ("match","team","t"))
    use "$TMP/act_all.dta", clear
    sort match team
    mata: ACT = st_data(., ("match","team","t_join","t_leave"))
    clear
    mata: WN = make_windows(ACT, GL, PS)
    mata: st_addobs(rows(WN)); (void) st_addvar("double", ("match","team","n_f","e_f","n_l","e_l","n_w","e_w","g_lm")); st_store(., ("match","team","n_f","e_f","n_l","e_l","n_w","e_w","g_lm"), WN)
    mata: mata drop PS GL ACT WN
    merge 1:1 match team using "$TMP/scores.dta", keep(match) nogenerate
    merge m:1 team using "$TMP/team_cov.dta", keep(master match) nogenerate keepusing(z_ppg)
    rename z_ppg zself_unused
    rename opp_team oteam
    rename team own_team
    rename oteam team
    merge m:1 team using "$TMP/team_cov.dta", keep(master match) nogenerate keepusing(ppg)
    rename ppg ppg_o
    rename team opp_team
    rename own_team team
    drop zself_unused
    save "$TMP/outcome_windows.dta", replace

    tempfile res
    tempname ph
    postfile `ph' str10 kind str24 spec str14 term double(n b se) using `res', replace
    foreach kind in first landmark whole {
        use "$TMP/outcome_windows.dta", clear
        if "`kind'"=="first" {
            keep if n_f>=30 & e_f>0
            gen double tempo = ln(n_f/e_f)
        }
        else if "`kind'"=="landmark" {
            keep if g_lm==0 & n_l>0
            gen double tempo = ln(n_l/e_l)
        }
        else {
            gen double tempo = ln(n_w/e_w)
        }
        keep match team opp_team gd home ppg_o tempo
        * opponent's tempo in the same match (must also be in this window set)
        preserve
        keep match team tempo
        rename (team tempo) (opp_team tempo_o)
        tempfile oth
        save `oth'
        restore
        merge m:1 match opp_team using `oth', keep(match) nogenerate
        drop if missing(ppg_o)
        local n = _N
        quietly summarize tempo
        local sdb = sqrt(r(Var)*(r(N)-1)/r(N))
        local mb = r(mean)
        gen double b_s = (tempo - `mb')/`sdb'
        quietly summarize tempo_o
        gen double bo_s = (tempo_o - r(mean))/`sdb'
        quietly summarize ppg_o
        gen double opp_s = (ppg_o - r(mean))/sqrt(r(Var)*(r(N)-1)/r(N))
        quietly regress gd i.team home b_s opp_s, vce(cluster team)
        post `ph' ("`kind'") ("own tempo") ("home") (`n') (_b[home]) (_se[home])
        post `ph' ("`kind'") ("own tempo") ("own tempo") (`n') (_b[b_s]) (_se[b_s])
        post `ph' ("`kind'") ("own tempo") ("opp strength") (`n') (_b[opp_s]) (_se[opp_s])
        quietly regress gd i.team home b_s bo_s opp_s, vce(cluster team)
        post `ph' ("`kind'") ("own + opponent tempo") ("home") (`n') (_b[home]) (_se[home])
        post `ph' ("`kind'") ("own + opponent tempo") ("own tempo") (`n') (_b[b_s]) (_se[b_s])
        post `ph' ("`kind'") ("own + opponent tempo") ("opp tempo") (`n') (_b[bo_s]) (_se[bo_s])
        post `ph' ("`kind'") ("own + opponent tempo") ("opp strength") (`n') (_b[opp_s]) (_se[opp_s])
    }
    postclose `ph'
    use `res', clear
    gen double z = b/se
    save "$TMP/outcome_res.dta", replace
    export delimited using "$OUT/TableS7_tempo_goal_difference.csv", replace
}
*==============================================================================
* 8. ONE-STAGE CONDITIONAL LOGIT WITH CONTEXT INTERACTIONS
*    Each pass with an inferred receiver is one choice among the teammates on the pitch.
*    Statistics: depth difference, repeat, return;  interactions with home, opponent strength,
*    team strength, leading, trailing, late, leading x late, trailing x late (main effects of context
*    drop out of the conditional logit).  Newton-Raphson on the conditional likelihood, club-clustered
*    sandwich variance C/(C-1) (identical to clogit, vce(cluster club)); predicted shares under
*    counterfactual contexts by the delta method.
*    Feeds: Supplementary Tables S8 (coefficients) and S9 (predicted shares and contrasts).
*    (Implemented in Mata, because the 1.3 million
*     choice sets with up to 16 alternatives would need ~14 million long-format rows for clogit.)
*==============================================================================
mata:
// design matrix of one candidate slot: (base stats, base x ctx) ; base n x 3, ctx n x 8
real matrix xdes(real matrix b, real matrix cx)
{
    real matrix X
    real scalar m
    X = J(rows(b), 27, .)
    X[|1,1 \ .,3|] = b
    for (m = 1; m <= 8; m++) X[|1,3*m+1 \ .,3*m+3|] = b :* cx[,m]
    return(X)
}

// build the choice data from act (positions, on-pitch spells), events and goals
// ACT = (match,team,player,avgx,t_join,t_leave) sorted by match,team;  EV = (match,team,t,sender,receiver,home,z_opp,z_team)
// sorted by match,team,t,sender(str),receiver(str);  GL = (match,team,t) sorted by match
void clg_build(real matrix ACT, real matrix EV, real matrix GL, real scalar KM,
               real matrix BA, real matrix CX, real colvector CH, real colvector NK, real colvector CL, real rowvector info)
{
    real matrix evinfo, actinfo, glinfo, seen, bm
    real colvector pid, pax, ptj, ptl, uid, uidx, axu, gt, gteam, ci, pos, sel
    real scalar e, r, r1, r2, q, a1, a2, ai, key, nu, u, us, ur, nev, nd, N, nc, i, uc, m, tt, s, rc, ld, tr, late, diff, hm, zo, zt, Nmax, ok_s, j, tm
    transmorphic AM, GM
    evinfo = panelsetup(EV[,1]:*1e7 :+ EV[,2], 1); actinfo = panelsetup(ACT[,1]:*1e7 :+ ACT[,2], 1); glinfo = panelsetup(GL[,1], 1)
    AM = asarray_create("real", 1); GM = asarray_create("real", 1)
    for (i = 1; i <= rows(actinfo); i++) asarray(AM, ACT[actinfo[i,1],1]*1e7 + ACT[actinfo[i,1],2], i)
    for (i = 1; i <= rows(glinfo); i++) asarray(GM, GL[glinfo[i,1],1], i)
    Nmax = rows(EV)
    BA = J(Nmax, 3*KM, 0); CX = J(Nmax, 8, .); CH = J(Nmax, 1, .); NK = J(Nmax, 1, .); CL = J(Nmax, 1, .)
    nev = 0; nd = 0; N = 0
    for (e = 1; e <= rows(evinfo); e++) {
        r1 = evinfo[e,1]; r2 = evinfo[e,2]
        m = EV[r1,1]; tm = EV[r1,2]
        key = m*1e7 + tm
        if (!asarray_contains(AM, key)) continue
        ai = asarray(AM, key); a1 = actinfo[ai,1]; a2 = actinfo[ai,2]
        pid = ACT[|a1,3 \ a2,3|]; pax = ACT[|a1,4 \ a2,4|]; ptj = ACT[|a1,5 \ a2,5|]; ptl = ACT[|a1,6 \ a2,6|]
        uid = uniqrows(pid); nu = rows(uid)
        uidx = J(rows(pid), 1, .); axu = J(nu, 1, .)
        for (q = 1; q <= rows(pid); q++) {
            u = sum(uid:<pid[q]) + 1; uidx[q] = u; axu[u] = pax[q]
        }
        seen = J(nu, nu, 0)
        if (asarray_contains(GM, m)) {
            j = asarray(GM, m); gt = GL[|glinfo[j,1],3 \ glinfo[j,2],3|]; gteam = GL[|glinfo[j,1],2 \ glinfo[j,2],2|]
        }
        else {
            gt = J(0,1,.); gteam = J(0,1,.)
        }
        hm = EV[r1,6]; zo = EV[r1,7]; zt = EV[r1,8]
        for (r = r1; r <= r2; r++) {
            nev++
            tt = EV[r,3]; s = EV[r,4]; rc = EV[r,5]
            u = sum(uid:<s) + 1
            ok_s = (u <= nu) ? (uid[u] == s) : 0
            us = u
            ci = selectindex((ptj:<=tt) :& (tt:<ptl) :& (pid:!=s)); nc = rows(ci)
            pos = (nc > 0) ? selectindex(pid[ci]:==rc) : J(0,1,.)
            if (!ok_s | nc < 1 | nc > KM | rows(pos) == 0) {
                nd++
                if (ok_s) {
                    ur = sum(uid:<rc) + 1
                    if (ur <= nu) {
                        if (uid[ur] == rc) seen[us, ur] = 1
                    }
                }
                continue
            }
            diff = sum((gt:<=tt) :& (gteam:==tm)) - sum((gt:<=tt) :& (gteam:!=tm))
            late = (tt >= 3600); ld = (diff > 0); tr = (diff < 0)
            N++
            for (i = 1; i <= nc; i++) {
                uc = uidx[ci[i]]
                BA[N, 3*(i-1)+1] = abs(axu[us] - axu[uc])
                BA[N, 3*(i-1)+2] = seen[us, uc]
                BA[N, 3*(i-1)+3] = seen[uc, us]
            }
            CH[N] = pos[1]; NK[N] = nc; CL[N] = tm
            CX[N,.] = (hm, zo, zt, ld, tr, late, ld*late, tr*late)
            ur = sum(uid:<rc) + 1
            seen[us, ur] = 1
        }
    }
    BA = BA[|1,1 \ N,.|]; CX = CX[|1,1 \ N,.|]; CH = CH[|1 \ N|]; NK = NK[|1 \ N|]; CL = CL[|1 \ N|]
    // sort by club; club index 1..C
    sel = order(CL, 1)
    BA = BA[sel,.]; CX = CX[sel,.]; CH = CH[sel]; NK = NK[sel]; CL = CL[sel]
    info = (nev, nd, N, rows(uniqrows(CL)))
}

// choice probabilities of a chunk
real matrix clg_probs(real matrix BAc, real matrix CXc, real colvector NKc, real colvector beta, real scalar KM)
{
    real matrix E
    real scalar j
    E = J(rows(BAc), KM, 0)
    for (j = 1; j <= KM; j++) E[,j] = exp(xdes(BAc[|1,3*j-2 \ .,3*j|], CXc)*beta) :* (NKc :>= j)
    return(E :/ rowsum(E))
}

// log likelihood, gradient, information (Hessian of -ll), club scores
void clg_stats(real matrix BA, real matrix CX, real colvector CH, real colvector NK, real colvector CLI, real scalar C,
               real colvector beta, real scalar hess, real scalar score,
               real scalar ll, real colvector g, real matrix H, real matrix S)
{
    real scalar N, KM, CHN, a, b, j, n, F, i
    real matrix BAc, CXc, P, xb, xc, Xj, sc, info
    real colvector CHc, NKc, CLc
    N = rows(BA); KM = cols(BA)/3; CHN = 16384; F = 27
    ll = 0; g = J(F,1,0); H = J(F,F,0); S = J(C,F,0)
    for (a = 1; a <= N; a = a + CHN) {
        b = min((a+CHN-1, N)); n = b-a+1
        BAc = BA[|a,1 \ b,.|]; CXc = CX[|a,1 \ b,.|]; CHc = CH[|a \ b|]; NKc = NK[|a \ b|]
        P = clg_probs(BAc, CXc, NKc, beta, KM)
        xb = J(n,F,0); xc = J(n,F,0)
        for (j = 1; j <= KM; j++) {
            Xj = xdes(BAc[|1,3*j-2 \ .,3*j|], CXc)
            xb = xb + P[,j] :* Xj
            xc = xc + (CHc:==j) :* Xj
            if (hess) H = H + cross(Xj, P[,j], Xj)
        }
        ll = ll + sum(ln(rowsum(P :* ((CHc*J(1,KM,1)) :== (J(n,1,1)*(1..KM))))))
        sc = xc - xb
        g = g + colsum(sc)'
        if (hess) H = H - cross(xb, xb)
        if (score) {
            CLc = CLI[|a \ b|]
            info = panelsetup(CLc, 1)
            for (i = 1; i <= rows(info); i++) S[CLc[info[i,1]],.] = S[CLc[info[i,1]],.] + colsum(sc[|info[i,1],1 \ info[i,2],.|])
        }
    }
}

// counterfactual prediction: mean predicted share to repeat / return partner and mean depth difference
// spec: rows (ctxcol, value)
void clg_predict(real matrix BA, real matrix CX, real colvector NK, real colvector beta, real matrix spec,
                 real rowvector res, real matrix gr)
{
    real scalar N, KM, CHN, a, b, j, n, F, k, c, mm
    real matrix BAc, CXc, P, xb, Xj, w
    real colvector NKc, sj
    N = rows(BA); KM = cols(BA)/3; CHN = 16384; F = 27
    res = J(1,3,0); gr = J(3,F,0)
    for (a = 1; a <= N; a = a + CHN) {
        b = min((a+CHN-1, N)); n = b-a+1
        BAc = BA[|a,1 \ b,.|]; CXc = CX[|a,1 \ b,.|]; NKc = NK[|a \ b|]
        for (k = 1; k <= rows(spec); k++) CXc[,spec[k,1]] = J(n,1,spec[k,2])
        if (anyof(spec[,1], 4) | anyof(spec[,1], 5) | anyof(spec[,1], 6)) {
            CXc[,7] = CXc[,4] :* CXc[,6]; CXc[,8] = CXc[,5] :* CXc[,6]
        }
        P = clg_probs(BAc, CXc, NKc, beta, KM)
        xb = J(n,F,0)
        for (j = 1; j <= KM; j++) xb = xb + P[,j] :* xdes(BAc[|1,3*j-2 \ .,3*j|], CXc)
        // statistic columns: 2 = repeat, 3 = return, 1 = depth
        for (c = 1; c <= 3; c++) {
            mm = (c==1 ? 2 : (c==2 ? 3 : 1))
            w = J(n,1,0)
            for (j = 1; j <= KM; j++) {
                sj = BAc[,3*(j-1)+mm]
                w = w + P[,j] :* sj
                Xj = xdes(BAc[|1,3*j-2 \ .,3*j|], CXc)
                gr[c,.] = gr[c,.] + colsum((P[,j] :* sj) :* Xj)
            }
            res[c] = res[c] + sum(w)
            gr[c,.] = gr[c,.] - colsum(w :* xb)
        }
    }
    res = res/N; gr = gr/N
}

// full estimation + predictions.  Results to Stata matrices.
void clg_run(real scalar KM)
{
    real matrix BA, CX, H, S, A, V, gr, GRL, spec, LEV, GRs
    real colvector CH, NK, CL, CLI, beta, g, step, cu
    real rowvector info, res
    real scalar it, ll, ll0, C, F, N, s, i, j, k, a, b, un
    real matrix CON
    external real matrix ACT_, EV_, GL_
    clg_build(ACT_, EV_, GL_, KM, BA, CX, CH, NK, CL, info)
    N = info[3]; C = info[4]; F = 27
    cu = uniqrows(CL); CLI = J(N,1,.)
    for (i = 1; i <= C; i++) CLI[selectindex(CL:==cu[i])] = J(sum(CL:==cu[i]),1,i)
    printf("events %g used %g dropped %g clubs %g\n", info[1], N, info[2], C)
    beta = J(F,1,0)
    for (it = 1; it <= 12; it++) {
        clg_stats(BA, CX, CH, NK, CLI, C, beta, 1, 0, ll, g, H, S)
        step = lusolve(H + 1e-9*I(F), g); beta = beta + step
        printf("iter %g  loglik %12.2f  max|step| %9.2e\n", it, ll, max(abs(step)))
        if (max(abs(step)) < 1e-6) break
    }
    clg_stats(BA, CX, CH, NK, CLI, C, beta, 1, 1, ll, g, H, S)
    A = invsym(H)
    V = A*cross(S, S)*A*(C/(C-1))
    clg_stats(BA, CX, CH, NK, CLI, C, J(F,1,0), 0, 0, ll0, g, H, S)
    st_matrix("CL_COEF", (beta, sqrt(diagonal(V))))
    st_matrix("CL_INFO", (info[1], N, info[2], C, ll, ll0))
    // scenarios (rows: ctx column, value);  ctx columns: 1 home 2 opp 3 team 4 lead 5 trail 6 late
    LEV = J(12, 3, .); GRs = J(36, F, .)
    for (s = 1; s <= 12; s++) {
        if (s==1)  spec = (1,0)
        if (s==2)  spec = (1,1)
        if (s==3)  spec = (2,0)
        if (s==4)  spec = (2,1)
        if (s==5)  spec = (3,0)
        if (s==6)  spec = (3,1)
        if (s==7)  spec = (4,0 \ 5,0 \ 6,0)
        if (s==8)  spec = (4,1 \ 5,0 \ 6,0)
        if (s==9)  spec = (4,0 \ 5,1 \ 6,0)
        if (s==10) spec = (4,0 \ 5,0 \ 6,1)
        if (s==11) spec = (4,1 \ 5,0 \ 6,1)
        if (s==12) spec = (4,0 \ 5,1 \ 6,1)
        clg_predict(BA, CX, NK, beta, spec, res, gr)
        LEV[s,.] = res
        GRs[|3*(s-1)+1,1 \ 3*s,F|] = gr
    }
    // contrasts (a - b) over scenario indices: 1 away 2 home 3 opp0 4 opp+1 5 team0 6 team+1 7 level early 8 leading early
    // 9 trailing early 10 level late 11 leading late 12 trailing late
    CON = (2,1 \ 4,3 \ 6,5 \ 8,7 \ 9,7 \ 11,10 \ 12,10 \ 10,7 \ 11,8 \ 12,9)
    un = 1   // units: shares in percentage points (x100); depth per 10 pitch points is already in units of 10 points -> x10
    GRL = J(10, 6, .)
    for (k = 1; k <= rows(CON); k++) {
        a = CON[k,1]; b = CON[k,2]
        for (j = 1; j <= 3; j++) {
            gr = GRs[3*(a-1)+j,.] - GRs[3*(b-1)+j,.]
            GRL[k, 2*j-1] = (LEV[a,j] - LEV[b,j])
            GRL[k, 2*j]   = sqrt(gr*V*gr')
        }
    }
    st_matrix("CL_LEV", LEV); st_matrix("CL_CON", GRL)
}
end

program define run_clogit
    * args: suffix ("" or "strict")  -> reads act_<league>_prior<suffix>.csv
    args suf
    local sfx = cond("`suf'"=="","","_strict")
    local lgs : environment PACK_CL_LEAGUES
    if "`lgs'"=="" local lgs "$FIVE"
    tempfile actf evf
    local i 0
    foreach lg of local lgs {
        quietly import delimited using "$M/act_`lg'_prior`suf'.csv", clear varnames(1) case(preserve) asdouble
        keep match team player avgx t_join t_leave
        local ++i
        if `i'>1 quietly append using `actf'
        quietly save `actf', replace
    }
    sort match team, stable
    mata: ACT_ = st_data(., ("match","team","player","avgx","t_join","t_leave"))
    keep match team
    duplicates drop
    tempfile actkeys
    save `actkeys'
    local i 0
    foreach lg of local lgs {
        quietly import delimited using "$M/ev_`lg'.csv", clear varnames(1) case(preserve) asdouble stringcols(3 4)
        local ++i
        if `i'>1 quietly append using `evf'
        quietly save `evf', replace
    }
    * (team-matches with covariates and with positions only; events sorted by t, sender, receiver)
    merge m:1 match team using `actkeys', keep(match) nogenerate
    merge m:1 match team using "$TMP/match_cov.dta", keep(match) nogenerate keepusing(home z_opp z_team)
    gen double snd = real(sender)
    gen double rcv = real(receiver)
    sort match team t sender receiver
    mata: EV_ = st_data(., ("match","team","t","snd","rcv","home","z_opp","z_team"))
    use "$TMP/goals.dta", clear
    sort match t team
    mata: GL_ = st_data(., ("match","team","t"))
    mata: clg_run(16)
    mata: mata drop ACT_ EV_ GL_
    * ---- coefficient table ----
    clear
    local ctx "home opp team lead trail late lead_late trail_late"
    local names "depth repeat return"
    foreach c of local ctx {
        foreach b in depth repeat return {
            local names "`names' `b'_x_`c'"
        }
    }
    set obs 27
    gen str24 name = ""
    gen double b = .
    gen double se = .
    local j 0
    foreach n of local names {
        local ++j
        replace name = "`n'" in `j'
        replace b  = CL_COEF[`j',1] in `j'
        replace se = CL_COEF[`j',2] in `j'
    }
    gen double z = b/se
    gen double events  = CL_INFO[1,1]
    gen double used    = CL_INFO[1,2]
    gen double dropped = CL_INFO[1,3]
    gen double clubs   = CL_INFO[1,4]
    gen double loglik  = CL_INFO[1,5]
    gen double loglik_null = CL_INFO[1,6]
    export delimited using "$OUT/TableS8_conditional_logit_coef`sfx'.csv", replace
    save "$TMP/clogit_coef`sfx'.dta", replace
    * ---- predicted levels and contrasts ----
    clear
    set obs 22
    gen str40 row = ""
    gen str12 kind = ""
    gen double repeat_pp = .
    gen double repeat_se = .
    gen double return_pp = .
    gen double return_se = .
    gen double depth_pts = .
    gen double depth_se = .
    local sc "away home opp0 opp+1 team0 team+1 level_early leading_early trailing_early level_late leading_late trailing_late"
    local j 0
    foreach s of local sc {
        local ++j
        replace row = "`s'" in `j'
        replace kind = "level" in `j'
        replace repeat_pp = 100*CL_LEV[`j',1] in `j'
        replace return_pp = 100*CL_LEV[`j',2] in `j'
        replace depth_pts = 10*CL_LEV[`j',3] in `j'
    }
    local cn `" "home vs away" "opponent +1 SD vs mean" "team strength +1 SD vs mean" "leading vs level, early" "trailing vs level, early" "leading vs level, late" "trailing vs level, late" "late vs early, level" "late vs early, leading" "late vs early, trailing" "'
    local k 0
    foreach c of local cn {
        local ++k
        local r = 12 + `k'
        replace row = "`c'" in `r'
        replace kind = "contrast" in `r'
        replace repeat_pp = 100*CL_CON[`k',1] in `r'
        replace repeat_se = 100*CL_CON[`k',2] in `r'
        replace return_pp = 100*CL_CON[`k',3] in `r'
        replace return_se = 100*CL_CON[`k',4] in `r'
        replace depth_pts = 10*CL_CON[`k',5] in `r'
        replace depth_se  = 10*CL_CON[`k',6] in `r'
    }
    export delimited using "$OUT/TableS9_conditional_logit_predicted`sfx'.csv", replace
    save "$TMP/clogit_pred`sfx'.dta", replace
end

if $RUN_clogit {
    di as txt _n "=== block clogit: " c(current_time)
    run_clogit
    run_clogit strict
}
*==============================================================================
* 9. FIGURES
*    Figure 1: variation of the team-match estimates around the pooled effect (tag JP).
*    Figure 2: frequency and observed outcomes in football units (home, opponent, team strength).
*    Needs $TMP/pooled.dta, rate_res.dta, obs_res.dta (blocks pool, rate, observed).
*==============================================================================
mata:
// percentile with linear interpolation of a column vector
real scalar nppct(real colvector x, real scalar p)
{
    real colvector s
    real scalar pos, lo, hi
    s = sort(x, 1); pos = (rows(s)-1)*p/100; lo = floor(pos); hi = ceil(pos)
    return(s[lo+1] + (s[hi+1]-s[lo+1])*(pos-lo))
}
end

* TIF: Stata for Mac has no built-in TIF translator in this installation (Graph2tif missing);
* try the native export first and otherwise convert the PNG with the macOS tool sips (LZW, dpi set).
program define export_tif
    args name width dpi
    capture graph export "$OUT/`name'.tif", name(`name') width(`width') replace
    if _rc {
        quietly shell sips -s format tiff -s formatOptions lzw -s dpiWidth `dpi' -s dpiHeight `dpi' "$OUT/`name'.png" --out "$OUT/`name'.tif" >/dev/null 2>&1
        capture confirm file "$OUT/`name'.tif"
        if _rc di as err "TIF export failed for `name' (no native translator and no sips)"
    }
end

* tick labels with a leading zero (0.5 instead of .5): ticks from _natscale within [lo,hi]
capture program drop mklab
program define mklab, rclass
    args lo hi
    quietly _natscale `lo' `hi' 5
    local mn = r(min)
    local dl = r(delta)
    local s ""
    local v = `mn'
    while `v' <= `hi' + 1e-9 {
        if `v' >= `lo' - 1e-9 {
            local t = strofreal(`v', "%9.0g")
            if substr("`t'",1,1)=="." local t "0`t'"
            if substr("`t'",1,2)=="-." local t "-0`=substr("`t'",2,.)'"
            local s `"`s' `v' "`t'""'
        }
        local v = round(`v' + `dl', `dl'/1000)
    }
    return local lab `"`s'"'
end

if $RUN_figures {
    di as txt _n "=== block figures: " c(current_time)
    set graphics off
    graph set eps fontface Helvetica
    set scheme s1mono      // minimal black-and-white scheme; all colours below are black or grey
    * ---------------- Figure 1 ----------------
    local fl ""
    foreach lg of global FIVE {
        local fl `"`fl' "$M/est_`lg'_JP.csv""'
    }
    load_est `fl'
    tempfile est
    save `est'
    use "$TMP/pooled.dta", clear
    keep if spec=="JP"
    tempfile pj
    save `pj'
    local panels `" "diff|Depth difference" "sim|Lateral similarity" "inertia|Repeat passing" "recip|Return passing" "alter|Receiver rating" "tie_ten1|Shared tenure, 1+ year" "'

    local letters "A B C D E F"
    local gl ""
    local p 0
    foreach pn of local panels {
        local ++p
        local t = substr("`pn'",1,strpos("`pn'","|")-1)
        local ttl = substr("`pn'",strpos("`pn'","|")+1,.)
        local L : word `p' of `letters'
        use `pj', clear
        keep if term=="`t'"
        local th = pth[1]
        local tau = tau[1]
        local kk = k[1]
        use `est', clear
        keep if term=="`t'"
        mata: x = st_data(., "b")
        mata: st_numscalar("p1", nppct(x,1)); st_numscalar("p99", nppct(x,99))
        mata: mata drop x
        scalar lo = min(p1, `th'-2.2*`tau', 0)
        scalar hi = max(p99, `th'+2.2*`tau', 0)
        scalar pad = .04*(hi-lo)
        scalar lo = lo - pad
        scalar hi = hi + pad
        local lo = lo
        local hi = hi
        local w = (hi-lo)/40
        local a1 = `th' - 1.96*`tau'
        local a2 = `th' + 1.96*`tau'
        local thf : display %5.3f `th'
        local tauf : display %5.3f `tau'
        local yt ""
        if inlist("`L'","A","D") local yt "Team-matches"
        mklab `lo' `hi'
        local xlab `"`r(lab)'"'
        twoway (histogram b if inrange(b,`lo',`hi'), frequency start(`lo') width(`w') fcolor(white) lcolor(black) lwidth(vthin)) ///
               , xline(`a1' `a2', lcolor(black) lpattern(dot) lwidth(medthin)) ///
                 xline(0, lcolor(gs10) lwidth(thin) lpattern(solid)) xline(`th', lcolor(black) lwidth(medthick) lpattern(solid)) ///
                 title("{bf:`L'  `ttl'}", position(11) size(medsmall) color(black)) ///
                 xtitle("pooled `thf', {&tau} = `tauf'", size(small) color(black)) ///
                 ytitle("`yt'", size(small) color(black)) ylabel(none) yscale(lstyle(none)) ytick(none) ///
                 xlabel(`xlab', labsize(small) labcolor(black) nogrid) xscale(range(`lo' `hi') lcolor(black)) ///
                 legend(off) graphregion(color(white) margin(small)) plotregion(color(white) style(none)) ///
                 name(f1_`p', replace) nodraw
        local gl "`gl' f1_`p'"
    }
    graph combine `gl', cols(3) xsize(6.5) ysize(4.4) graphregion(color(white) margin(small)) imargin(small) name(Figure1, replace)
    graph export "$OUT/Figure1.png", name(Figure1) width(3900) replace
    graph export "$OUT/Figure1.eps", name(Figure1) replace
    export_tif Figure1 3900 600

    * ---------------- Figure 2 ----------------
    * One panel: pass frequency and five observed outcomes, as change relative to the outcome's mean (%),
    * for home vs away, opponent strength +1 SD and team strength +1 SD.
    use "$TMP/rate_res.dta", clear
    keep if outcome=="all attempts" & ((model=="A" & inlist(term,"home","opp")) | (model=="A2" & term=="team"))
    gen double est = 100*(rr-1)
    gen double elo = 100*(lo-1)
    gen double ehi = 100*(hi-1)
    gen str8 var = "freq"
    keep var term est elo ehi
    tempfile f2
    save `f2'
    use "$TMP/obs_res.dta", clear
    rename (outcome b) (var est)
    bysort var: egen double __base = max(cond(term=="baseline", est, .))
    keep if term!="baseline"
    gen double elo = est - 1.96*se
    gen double ehi = est + 1.96*se
    foreach v in est elo ehi {
        replace `v' = 100*`v'/__base
    }
    keep var term est elo ehi
    append using `f2'
    gen double row = .
    replace row = 6 if var=="freq"
    replace row = 5 if var=="ident"
    replace row = 4 if var=="same"
    replace row = 3 if var=="ret"
    replace row = 2 if var=="length"
    replace row = 1 if var=="fwd"
    gen byte k = cond(term=="home",1,cond(term=="opp",2,3))
    gen double y = row + cond(k==1, .17, cond(k==2, 0, -.17))
    save "$TMP/fig2_data.dta", replace
    local c1 "0 114 178"
    local c2 "213 94 0"
    local c3 "0 158 115"
    twoway (rcap elo ehi y if k==1, horizontal lcolor("`c1'") lwidth(medthick) msize(0)) ///
           (rcap elo ehi y if k==2, horizontal lcolor("`c2'") lwidth(medthick) msize(0)) ///
           (rcap elo ehi y if k==3, horizontal lcolor("`c3'") lwidth(medthick) msize(0)) ///
           (scatter y est if k==1, mcolor("`c1'") msymbol(O) msize(medium)) ///
           (scatter y est if k==2, mcolor("`c2'") msymbol(D) msize(medium)) ///
           (scatter y est if k==3, mcolor("`c3'") msymbol(S) msize(medium)) ///
           , xline(0, lcolor(gs8) lwidth(thin)) ///
             ylabel(6 "Pass frequency" 5 "Receiver inferred" 4 "Same-role share" 3 "Return share" 2 "Pass length" 1 "Forward displacement", angle(0) labsize(small) noticks nogrid) ///
             yscale(range(.35 6.65)) ///
             xlabel(-20(5)15, labsize(small) grid glcolor(gs14) glpattern(solid)) ///
             xtitle("Change relative to the outcome's mean, %", size(small)) ytitle("") ///
             legend(order(4 "Home vs away" 5 "Opponent +1 SD" 6 "Team strength +1 SD") rows(1) position(6) size(small) colgap(*3) symxsize(*.6) region(lstyle(none))) ///
             graphregion(color(white)) plotregion(color(white)) xsize(8) ysize(6.4) name(Figure2, replace)
    graph export "$OUT/Figure2.png", name(Figure2) width(1950) replace
    graph export "$OUT/Figure2.eps", name(Figure2) replace
    export_tif Figure2 1950 300
}


*==============================================================================
* 4. CONTEXT CONTRASTS FOR SIMPLE PASSES ONLY  (Supplementary Table S11)
*    Same second-stage regressions as Table S5 (strength model, league fixed effects, club-clustered
*    variance), on team-match estimates from the model without tenure terms, earlier-match positions,
*    fitted to all passes with a receiver (VP_base) and to simple passes only (VP_simple).
*==============================================================================
if $RUN_simplectx {
    di as txt _n "=== block simplectx: " c(current_time)
    tempfile sres
    tempname ps
    postfile `ps' str8 variant str12 term str10 moderator double(b se n tau) using `sres', replace
    foreach vr in base simple {
        local fl ""
        foreach lg of global FIVE {
            local fl `"`fl' "$V/prior_`vr'_`lg'/est_`lg'_VP.csv""'
        }
        load_est `fl'
        rename file lgn
        quietly merge m:1 match team using "$TMP/match_cov.dta", keep(match) nogenerate
        quietly drop if missing(z_team) | missing(z_opp)
        foreach t in inertia recip indeg trans same diff sim alter {
            mata: st_matrix("MR", mod_run("`t'", (1,3,4), "league"))
            local j 0
            foreach mname in strength home opponent {
                local ++j
                post `ps' ("`vr'") ("`t'") ("`mname'") (MR[1,2+`j']) (MR[1,5+`j']) (MR[1,1]) (MR[1,2])
            }
        }
    }
    postclose `ps'
    use `sres', clear
    gen double z = b/se
    export delimited using "$OUT/TableS11_simple_pass_moderators.csv", replace
    list, sepby(term) noobs
}

*==============================================================================
* END
*==============================================================================
timer off 99
quietly timer list 99
di as txt _n "3_analyze.do finished; run time " %6.1f r(t99) " seconds (" %5.1f r(t99)/60 " minutes)"
