*==============================================================================
* 2_fit.do   Stage 2 of the replication pack: model inputs and fits
*
* WHAT IT DOES
*   Stage "inputs": reads the build files (the csv files of 1_build.do in output/build/, data/raw/players.json and the match json files in data/raw/matches/)
*       and writes the model input files to OUT (default PACK/output/fit):
*       ev_<lg>.csv (pass events), act_<lg>.csv (player-match records), dy_<lg>.csv (dyad
*       covariates incl. shared minutes and exact-tenure bands), act_<lg>_prior(strict).csv
*       (earlier-match position profiles).
*   Stage "fits": one dynamic network actor model (nwdynam) per team-match, for every
*       requested specification ("job") and league; writes est_<lg>_<tag>.csv to OUT
*       (variant jobs: OUT/variants/<variant>_<lg>/est_<lg>_V.csv, see below).
*
* HOW TO RUN (from any directory)
*   global PACK "/path/to/replication_pack"          // or environment variable PACK_ROOT
*   global NWDIR "/path/to/folder_with_nwcommands_and_nwdynam"   // or environment NWDIR;
*                                                    // if empty, installed packages are used
*   do "$PACK/stata/2_fit.do"
*   Batch:  stata-mp -b do 2_fit.do      (set the globals in a small wrapper do-file or
*           export PACK_ROOT and NWDIR in the shell first)
*
* SWITCHES (set before the do-file is run; defaults in brackets)
*   global STAGES        "inputs fits"   [inputs fits]   either word or both
*   global LEAGUES       "England"       [all five]  leagues to be FITTED
*   global INPUT_LEAGUES                 [all five]  leagues for the input files (the earlier-match
*                                         position files need all five; leave the default)
*   global TAGS          "J JP"          [paper set, below]  jobs to be fitted
*   global MAXN          5               [99999]     only the first MAXN team-matches (test run)
*   global PARALLEL      1               [1]  fit the leagues in separate Stata batch processes
*                                         (only if more than one league is fitted; 0 = loop)
*   global OUT           "..."           [PACK/output/fit]
*   global STATA_EXE     "..."           [stata binary of the running flavor]
*
*   Jobs (TAGS): J U V R X B (match positions), JP XP BP RP (earlier-match positions),
*     JS (strict fallback, Table S10), V_simple V_base V_loose V_t5 V_acc (tag V),
*     VP_simple VP_base VP_loose VP_t5 VP_acc (tag VP; Table 9, Table S2).
*     Default set = all of these.  loose, t5, acc exist for England only.
*
* nwcommands/nwdynam: the paper's estimates used nwdynam 1.9.1 with nwcommands 3.26.0 (0_install.do
* installs nwdynam 1.9.1 by default). NWDIR is optional: if set, it must contain the folders nwcommands/
* and nwdynam/ and is used instead of the installed packages. Newer releases (nwdynam 1.11.1) gave
* the same estimates in the checks made, apart from the unstable team-matches described in the README.
* NWDIR/versions are recorded in OUT/versions.txt.  Parallel runs need a Unix-like Stata (Mac/Linux);
* on Windows set PARALLEL 0.  Worker logs: OUT/_work/<league>/2_fit.log.
*
* Run time (Mac, 18 cores): inputs about 50 s; one league and one specification about 1 minute.
*
* No random numbers are used anywhere. Version 17 syntax.
*==============================================================================
version 17
set more off
set linesize 120

*------------------------------------------------------------------------------
* 0. Settings
*------------------------------------------------------------------------------
local isChild : env FIT_CHILD
if "`isChild'" == "1" {                      // worker process started by run_parallel
    global PACK      : env PACK_ROOT
    global NWDIR     : env NWDIR
    global OUT       : env FIT_OUT
    global LEAGUES   : env FIT_LG
    global TAGS      : env FIT_TAGS
    global MAXN      : env FIT_MAXN
    global STAGES    "fits"
    global PARALLEL  0
}
if "$PACK" == "" global PACK : env PACK_ROOT
if "$PACK" == "" {
    capture confirm file "`c(pwd)'/stata/2_fit.do"
    if !_rc global PACK "`c(pwd)'"
    else {
        capture confirm file "`c(pwd)'/../stata/2_fit.do"
        if !_rc global PACK "`c(pwd)'/.."
    }
}
if "$PACK" == "" {
    display as error "set global PACK (the replication pack folder with stata/ and data/) or the environment variable PACK_ROOT"
    exit 198
}
if "$NWDIR" == "" global NWDIR : env NWDIR
if "$OUT" == ""      global OUT "$PACK/output/fit"
if "$BUILD" == ""    global BUILD : env BUILD_DIR
if "$BUILD" == ""    global BUILD "$PACK/output/build"   // output of 1_build.do (set BUILD to use other build files)
if "$STAGES" == ""   global STAGES "inputs fits"
if "$LEAGUES" == ""  global LEAGUES "England Spain Italy Germany France"
if "$INPUT_LEAGUES" == "" global INPUT_LEAGUES "England Spain Italy Germany France"
if "$MAXN" == ""     global MAXN 99999
if "$PARALLEL" == "" global PARALLEL 1
if "$TAGS" == "" global TAGS "J U V R X B JP XP BP RP JS V_simple V_base V_loose V_t5 V_acc VP_simple VP_base VP_loose VP_t5 VP_acc"
capture mkdir "$PACK/output"
capture mkdir "$OUT"
capture mkdir "$OUT/_work"
capture mkdir "$OUT/variants"
global W "$OUT/_work"

*------------------------------------------------------------------------------
* 1. Mata helper (used by the fits): dense matrix of a dyad covariate in the node order _lab
*------------------------------------------------------------------------------
mata:
void mkmat(string scalar nm, string scalar vvar)
{
    external string colvector _lab
    string colvector I, J
    real colvector v, ia, ib
    real matrix M
    real scalar r, n
    n = rows(_lab)
    I = st_sdata(., "i"); J = st_sdata(., "j"); v = st_data(., vvar)
    M = J(n, n, 0)
    for (r = 1; r <= rows(I); r++) {
        ia = selectindex(_lab :== I[r]); ib = selectindex(_lab :== J[r])
        if (rows(ia)==1 & rows(ib)==1) M[ia, ib] = v[r]
    }
    st_matrix(nm, M)
}
end

*------------------------------------------------------------------------------
* 2. INPUT FILES
*------------------------------------------------------------------------------

* 2a. Static tables used by all leagues: nodes, dyads, FIFA ratings with role-wise z scores,
*     match meta data (date, teams) and player roles from the Wyscout json files.
capture program drop fx_static
program define fx_static
    local P "$PACK"
    * --- nodes (all columns kept as text so that pass-through values are written back unchanged)
    import delimited using "$BUILD/nodes.csv", clear stringcols(_all) varnames(1)
    gen double tin  = cond(t_in=="", 0, real(t_in))
    gen double tout = real(t_out)
    gen double prn  = real(prior_rank)
    save "$W/nodes.dta", replace
    import delimited using "$BUILD/dyads.csv", clear stringcols(_all) varnames(1)
    save "$W/dyads.dta", replace

    * --- players.json: role code2 per wyId (one record per line after the filter)
    tempfile pj mj
    filefilter "`P'/data/raw/players.json" `pj', from("}, {\QpassportArea\Q") to("}\n{\QpassportArea\Q") replace
    infix strL line 1-100000 using `pj', clear
    gen str12 p = ""
    replace p = ustrregexs(1) if ustrregexm(line, `""wyId": ([0-9]+)"')
    gen str8 role = ""
    replace role = ustrregexs(1) if ustrregexm(line, `""role": \{"code2": "([A-Za-z0-9]*)""')
    keep p role
    drop if p==""
    duplicates drop p, force
    tempfile roles
    save `roles'

    * --- FIFA overall rating (ov18, else ov19, else ov17), z scores overall and within role
    import delimited using "$BUILD/fifa_overall.csv", clear stringcols(_all) varnames(1)
    gen str10 ovs = cond(ov18!="", ov18, cond(ov19!="", ov19, ov17))
    keep if ovs != ""
    gen double ov = real(ovs)
    rename wy_id p
    merge 1:1 p using `roles', keep(master match)
    gen byte injson = (_merge==3)
    drop _merge
    summarize ov
    scalar GM = r(mean)
    scalar GS = r(sd)*sqrt((r(N)-1)/r(N))
    bysort role: gen long nr = _N
    bysort role: egen double rmean = mean(ov)
    bysort role: egen double rsd0 = sd(ov)
    gen double rsd = rsd0*sqrt((nr-1)/nr)
    gen double ovz = real(string((ov-GM)/GS, "%21.4f"))
    gen double ovr = 0
    replace ovr = real(string((ov-rmean)/rsd, "%21.4f")) if injson & nr>5
    gen byte rated = 1
    keep p ov ovz ovr rated
    save "$W/fifa.dta", replace

    * --- matches: id, date, the two teams (five leagues)
    tempfile mm
    local first 1
    foreach lg in England Spain Italy Germany France {
        filefilter "`P'/data/raw/matches/matches_`lg'.json" `mj', from("}, {\Qstatus\Q") to("}\n{\Qstatus\Q") replace
        infix strL line 1-100000 using `mj', clear
        gen str8 league = "`lg'"
        gen str10 match = ""
        replace match = ustrregexs(1) if ustrregexm(line, `""wyId": ([0-9]+)"')
        gen str20 dateutc = ""
        replace dateutc = ustrregexs(1) if ustrregexm(line, `""dateutc": "([^"]+)""')
        gen str10 t1 = ""
        replace t1 = ustrregexs(1) if ustrregexm(line, `""teamId": ([0-9]+)"')
        gen strL l2 = ustrregexrf(line, `""teamId": [0-9]+"', "")
        gen str10 t2 = ""
        replace t2 = ustrregexs(1) if ustrregexm(l2, `""teamId": ([0-9]+)"')
        keep league match dateutc t1 t2
        drop if match==""
        if `first' == 0 append using `mm'
        save `mm', replace
        local first 0
    }
    gen double mnum = real(match)
    sort dateutc mnum
    gen long mrank = _n                       // order of matches: date, then match id
    save "$W/matches.dta", replace
    keep league match t1 t2
    gen long id = _n
    reshape long t, i(id) j(k)
    rename t team
    drop id k
    save "$W/matchteams.dta", replace
end

* 2b. Model inputs of one league from a pass file: ev_, act_, dy_ (with shared minutes if sm is given)
capture program drop fx_export
program define fx_export
    syntax , LG(string) PASSFILE(string) OUTDIR(string) [SM]
    tempfile nl nkeys pall posf evq pl evkeys actf dyf

    * nodes of the league; mean and (population) sd of the prior rank over all node rows
    use "$W/nodes.dta", clear
    keep if league == "`lg'"
    save `nl'
    quietly summarize prn
    local RM = r(mean)
    local RS = r(sd)*sqrt((r(N)-1)/r(N))
    keep match team
    duplicates drop
    save `nkeys'

    * passes of the league
    import delimited using "`passfile'", clear varnames(1) stringcols(1 2 3 7 8) asdouble
    keep if league == "`lg'"
    keep match team t passer receiver x0 y0
    recast double t
    save `pall'

    * mean pass origin per passer and match (/10), written as avgx/avgy
    keep if !missing(x0) & !missing(y0)
    recast double x0 y0
    collapse (mean) ax=x0 ay=y0, by(match team passer)
    replace ax = ax/10
    replace ay = ay/10
    rename passer p
    save `posf'

    * events: passes with a known receiver, at least 30 per team-match, team-match in nodes
    use `pall', clear
    keep if receiver != ""
    rename passer sender
    keep match team sender receiver t
    order match team sender receiver t
    bysort match team: keep if _N >= 30
    merge m:1 match team using `nkeys', keep(match) nogen
    sort match team t sender receiver
    save `evq'
    export delimited using "`outdir'/ev_`lg'.csv", replace
    keep match team
    duplicates drop
    save `evkeys'

    * players seen in the events, first and last event time, first event of the team-match
    use `evq', clear
    keep match team sender t
    rename sender p
    tempfile pa
    save `pa'
    use `evq', clear
    keep match team receiver t
    rename receiver p
    append using `pa'
    collapse (min) first=t (max) last=t, by(match team p)
    bysort match team: egen double e0 = min(first)
    save `pl'

    * act_<lg>.csv
    use `nl', clear
    rename player p
    merge m:1 match team p using `pl', keep(match) nogen
    merge m:1 match team p using `posf', keep(master match) nogen
    merge m:1 p using "$W/fifa.dta", keep(master match) nogen
    replace ovz   = 0 if missing(ovz)
    replace ovr   = 0 if missing(ovr)
    replace rated = 0 if missing(rated)
    gen double tj = min(tin, first)
    gen double tl = max(tout, last + 0.001)
    gen byte present0 = (tj <= e0 + 1e-9)
    gen double t_join  = real(string(tj, "%21.3f"))
    gen double t_leave = real(string(tl, "%21.3f"))
    gen double avgx = ax
    gen double avgy = ay
    gen double rz = 0
    replace rz = real(string((prn - `RM')/`RS', "%21.4f")) if !missing(prn)
    gen byte rolec = .
    replace rolec = 0 if role=="GK"
    replace rolec = 1 if role=="DF"
    replace rolec = 2 if role=="MD"
    replace rolec = 3 if role=="FW"
    * star: the two best FIFA-rated players of the team-match (ties: larger player id string)
    gsort match team -rated -ov -p
    by match team: gen byte star = (_n <= 2 & rated == 1)
    rename p player
    keep match team player rolec age height prior_rank present0 t_join t_leave tenure_days tenure_src started avgx avgy rz ovz ovr rated star
    rename rolec role
    order match team player role age height prior_rank present0 t_join t_leave tenure_days tenure_src started avgx avgy rz ovz ovr rated star
    sort match team player
    export delimited using "`outdir'/act_`lg'.csv", replace

    * dy_<lg>.csv: both directions of every pair of players seen in the events
    tempfile pli plj
    use `pl', clear
    keep match team p
    rename p i
    save `pli'
    rename i j
    save `plj'
    use "$W/dyads.dta", clear
    merge m:1 match team using `evkeys', keep(match) nogen
    merge m:1 match team i using `pli', keep(match) nogen
    merge m:1 match team j using `plj', keep(match) nogen
    gen byte miss = (shared_tenure_days == "")
    gen double shn = cond(miss, 0, real(shared_tenure_days))
    gen byte ex = (shared_src == "tm_exact")
    gen byte ten1  = (shn >= 365  & !miss)
    gen byte ten3  = (shn >= 1095 & !miss)
    gen byte tenmiss = miss
    gen byte ten15 = (shn >= 548  & !miss)
    gen byte ten2  = (shn >= 730  & !miss)
    gen byte ten1x = (ex & shn >= 365)
    gen byte ten3x = (ex & shn >= 1095)
    gen byte tenmissx = !ex
    gen byte exb1 = (ex & shn >= 365 & shn < 1095)
    gen byte exb3 = (ex & shn >= 1095 & shn < 1825)
    gen byte exb5 = (ex & shn >= 1825)
    rename coapp_prior coapp
    rename same_role samerole
    rename shared_src src
    keep match team i j ten1 ten3 tenmiss ten15 ten2 ten1x ten3x tenmissx coapp samerole src exb1 exb3 exb5
    tempfile fwd
    save `fwd'
    rename (i j) (j i)
    append using `fwd'
    gen str20 pa = cond(i < j, i, j)
    gen str20 pb = cond(i < j, j, i)

    * shared minutes before the match: cumulative overlap of earlier matches of the team
    if "`sm'" != "" {
        tempfile pairs
        use `nl', clear
        keep team match date player tin tout
        rename (player tin tout) (pa tina touta)
        tempfile side1 side2
        save `side1'
        rename (pa tina touta) (pb tinb toutb)
        save `side2'
        use `side1', clear
        joinby match team using `side2'
        keep if pa < pb
        gen double ov = min(touta, toutb) - max(tina, tinb)
        gen double o60 = cond(ov > 0, ov/60, 0)
        gen double mnum = real(match)
        sort team pa pb date mnum
        by team pa pb: gen double cumb = 0
        by team pa pb: replace cumb = cumb[_n-1] + o60[_n-1] if _n > 1
        keep match team pa pb cumb
        duplicates drop match team pa pb, force
        save `pairs'
        use `fwd', clear
        rename (i j) (j i)
        append using `fwd'
        gen str20 pa = cond(i < j, i, j)
        gen str20 pb = cond(i < j, j, i)
        tempfile dyfull
        save `dyfull'
        merge m:1 match team pa pb using `pairs', keep(master match) nogen
        replace cumb = 0 if missing(cumb)
        sort match team cumb
        by match team: gen double q1 = cumb[floor(_N/3) + 1]
        by match team: gen double q2 = cumb[floor(2*_N/3) + 1]
        gen byte sm0   = (cumb > 0)
        gen byte sm270 = (cumb >= 270)
        gen byte sm900 = (cumb >= 900)
        gen byte smhi  = (cumb > q2 & cumb > 0)
        gen byte smlo  = (cumb > 0 & cumb <= q1)
        drop cumb q1 q2 pa pb
        order match team i j ten1 ten3 tenmiss ten15 ten2 ten1x ten3x tenmissx coapp samerole src sm0 sm270 sm900 smhi smlo exb1 exb3 exb5
    }
    else {
        drop pa pb
        order match team i j ten1 ten3 tenmiss ten15 ten2 ten1x ten3x tenmissx coapp samerole src exb1 exb3 exb5
    }
    sort match team i j
    export delimited using "`outdir'/dy_`lg'.csv", replace
end

* 2c. Earlier-match position profiles; strict = role fallback from earlier matches only
capture program drop fx_prior
program define fx_prior
    syntax , [STRICT]
    local MINP 10
    tempfile pm roles rb
    * pass origin sums per player and match, all five leagues
    import delimited using "$BUILD/passes_all.csv", clear varnames(1) stringcols(1 2 3 7 8) asdouble
    keep if inlist(league, "England", "Spain", "Italy", "Germany", "France")
    keep if !missing(x0) & !missing(y0)
    keep match passer x0 y0
    recast double x0 y0
    collapse (sum) sx=x0 sy=y0 (count) n=x0, by(passer match)
    rename passer player
    merge m:1 match using "$W/matches.dta", keep(master match) nogen keepusing(mrank)
    replace mrank = 1e9 if missing(mrank)
    save `pm'
    * role of each player (act files of all leagues, later league overrides)
    local first 1
    foreach lg in England Spain Italy Germany France {
        import delimited using "$OUT/act_`lg'.csv", clear stringcols(_all) varnames(1)
        keep player role
        if `first' == 0 append using `roles'
        save `roles', replace
        local first 0
    }
    duplicates drop player, force
    * (all players in all act files; duplicates across leagues keep one row)
    save `roles', replace
    * role-level sums per match (all passers; players without act row have role "")
    use `pm', clear
    merge m:1 player using `roles', keep(master match) nogen
    replace role = "" if missing(role)
    * overall role means (fallback) over all passes
    preserve
    collapse (sum) sx sy n, by(role)
    gen double rmx = sx/n/10
    gen double rmy = sy/n/10
    keep role rmx rmy
    tempfile rmean
    save `rmean'
    restore
    * role sums per match rank, cumulative before each match rank (strict fallback)
    collapse (sum) sx sy n, by(role mrank)
    tempfile rbm
    save `rbm'
    local nm = 0
    use "$W/matches.dta", clear
    keep mrank
    tempfile ranks
    save `ranks'
    tempfile rbc
    local firstr 1
    foreach r in x 0 1 2 3 {           // role codes of act_*.csv; x = no role
        use `ranks', clear
        gen str1 role = cond("`r'" == "x", "", "`r'")
        merge 1:1 role mrank using `rbm', keep(master match) nogen
        foreach v in sx sy n {
            replace `v' = 0 if missing(`v')
        }
        sort mrank
        foreach v in sx sy n {
            gen double c`v' = sum(`v') - `v'
        }
        keep role mrank csx csy cn
        if `firstr' == 0 append using `rbc'
        save `rbc', replace
        local firstr 0
    }

    foreach lg in England Spain Italy Germany France {
        import delimited using "$OUT/act_`lg'.csv", clear stringcols(_all) varnames(1)
        gen long ord = _n
        merge m:1 match using "$W/matches.dta", keep(master match) nogen keepusing(mrank)
        tempfile act
        save `act'
        * sums over the player's passes in earlier matches
        keep ord player mrank
        rename mrank arank
        joinby player using `pm'
        keep if mrank < arank
        collapse (sum) psx=sx psy=sy pn=n, by(ord)
        tempfile prev
        save `prev'
        use `act', clear
        merge 1:1 ord using `prev', keep(master match) nogen
        foreach v in psx psy pn {
            replace `v' = 0 if missing(`v')
        }
        merge m:1 role using `rmean', keep(master match) nogen
        replace rmx = 5 if missing(rmx)
        replace rmy = 5 if missing(rmy)
        if "`strict'" != "" {
            merge m:1 role mrank using `rbc', keep(master match) nogen
            replace rmx = cond(cn > 0 & !missing(cn), csx/cn/10, 5)
            replace rmy = cond(cn > 0 & !missing(cn), csy/cn/10, 5)
        }
        gen double nx = real(string(cond(pn >= `MINP', psx/pn/10, rmx), "%21.3f"))
        gen double ny = real(string(cond(pn >= `MINP', psy/pn/10, rmy), "%21.3f"))
        drop avgx avgy
        rename (nx ny) (avgx avgy)
        sort ord
        keep match team player role age height prior_rank present0 t_join t_leave tenure_days tenure_src started avgx avgy rz ovz ovr rated star
        order match team player role age height prior_rank present0 t_join t_leave tenure_days tenure_src started avgx avgy rz ovz ovr rated star
        export delimited using "$OUT/act_`lg'_prior`=cond("`strict'"!="","strict","")'.csv", replace
    }
end

*------------------------------------------------------------------------------
* 3. FITS
*------------------------------------------------------------------------------

* 3a. nwcommands/nwdynam on the adopath (the folder NWDIR holds nwcommands/ and nwdynam/)
capture program drop nw_setup
program define nw_setup
    local nwdir "$NWDIR"
    if "`nwdir'" == "" {
        capture which nwdynam
        if _rc {
            display as error "nwdynam not installed: set NWDIR to the folder with nwcommands/ and nwdynam/"
            exit 601
        }
        exit
    }
    local pid = subinstr("$LEAGUES", " ", "_", .)    // own PLUS/PERSONAL folders per worker
    capture mkdir "$W/plus_`pid'"
    capture mkdir "$W/pers_`pid'"
    sysdir set PLUS "$W/plus_`pid'"
    sysdir set PERSONAL "$W/pers_`pid'"
    foreach d in nwcommands nwcommands/lib nwcommands/lib/plugins/unix nwdynam nwdynam/lib nwdynam/lib/plugins/unix {
        adopath ++ "`nwdir'/`d'"
    }
end

* 3b. One league, one specification: one dynamic network actor model per team-match
capture program drop fit_one
program define fit_one
    syntax , LG(string) TAG(string) EVF(string) ACTF(string) DYF(string) OUTF(string) ///
        [NETS(string) SPEC(string) MAXN(integer 99999)]
    local tlast = 0
    tempfile evall actall dyall list res
    import delimited using "`evf'", clear
    tostring sender receiver, replace
    save `evall'
    bysort match team: keep if _n == 1
    keep match team
    gen k = _n
    keep if k <= `maxn'
    save `list'
    import delimited using "`actf'", clear stringcols(3)
    rename player _nwnode
    save `actall'
    import delimited using "`dyf'", clear stringcols(3 4)
    save `dyall'
    tempname P
    postfile `P' long match long team int n_ev int n_act int rc str40 term double b double se double ll using `res', replace
    use `list', clear
    local K = _N
    forvalues k = 1/`K' {
        use `list', clear
        local m = match[`k']
        local t = team[`k']
        quietly nwclear
        tempfile act_k dy_k
        quietly {
            use `actall' if match==`m' & team==`t', clear
            save `act_k'
            use `dyall' if match==`m' & team==`t', clear
            save `dy_k'
            use `evall' if match==`m' & team==`t', clear
        }
        local nev = _N
        quietly summ t
        local tlast = r(max)
        quietly nwset sender receiver, eventtime(t) name(pass)
        quietly nwload pass, xvars
        local nact = _N
        quietly {
            gen long ord = _n
            merge 1:1 _nwnode using `act_k', keep(master match) nogen
            sort ord
            summ avgx
            replace avgx = r(mean) if missing(avgx)
            summ avgy
            replace avgy = r(mean) if missing(avgy)
            replace rz = 0 if missing(rz)
            replace ovz = 0 if missing(ovz)
            replace ovr = 0 if missing(ovr)
            gen double ovz2 = ovz
            gen byte star2 = star
        }
        mata: _lab = st_sdata(., "_nwnode")
        capture noisily {
            local nt ""
            foreach nm in `nets' {
                preserve
                quietly use `dy_k', clear
                mata: mkmat("M_`nm'", "`nm'")
                restore
                quietly nwset, mat(M_`nm') directed name(`nm') nooutput
                mata: st_local("sm", strofreal(sum(st_matrix("M_`nm'"))))
                local full = `nact'*(`nact'-1)
                if `sm' > 0 & `sm' < `full' local nt "`nt' tie(`nm')"
            }
            local nn ""
            if "`nt'" != "" local nn "netterms(`nt')"
            * joins and leaves of the players
            tempfile cc
            local ccopt
            preserve
            quietly {
                use `act_k', clear
                gen str40 id2 = _nwnode
                gen byte lv = t_leave < `tlast'
                tempfile jn
                keep _nwnode id2 present0 lv t_join t_leave
                save `jn'
                keep if present0==0
                gen double time = t_join
                gen byte joined = 1
                keep _nwnode id2 joined time
                save `cc'
                use `jn', clear
                keep if lv==1
                gen double time = t_leave
                gen byte joined = 0
                keep _nwnode id2 joined time
                append using `cc'
                rename _nwnode id1
            }
            local ncomp = _N
            if `ncomp' > 0 {
                quietly nwset id1 id2 joined, eventtime(time) selfloop name(membership)
                local ccopt "compchange(membership)"
            }
            restore
            quietly nwdynam pass, `spec' `nn' present(present0) `ccopt'
        }
        local rc = _rc
        if `rc'==0 {
            matrix b = e(b)
            matrix V = e(V)
            local cn : colnames b
            local i = 0
            foreach c of local cn {
                local ++i
                post `P' (`m') (`t') (`nev') (`nact') (0) ("`c'") (b[1,`i']) (sqrt(V[`i',`i'])) (e(ll))
            }
        }
        else post `P' (`m') (`t') (`nev') (`nact') (`rc') ("") (.) (.) (.)
    }
    postclose `P'
    use `res', clear
    export delimited using "`outf'", replace
end

* 3c. Job table: specification of every tag -> locals returned with c_local
*     vr = passes variant ("" = main passes), act = position source ("", prior, priorstrict)
capture program drop job_def
program define job_def
    args job
    local H  "inertia recip indeg trans same(role) diff(avgx) sim(avgy) alter(ovz) egoalterint(ovz ovz2)"
    local B7 "inertia recip indeg trans same(role) diff(avgx) sim(avgy)"
    local tag "`job'"
    local nets ""
    local spec "`H'"
    local act ""
    local vr ""
    if "`job'" == "J"  local nets "ten1 sm270 sm900"
    if "`job'" == "U"  local spec "inertia recip indeg trans alter(ovz) egoalterint(ovz ovz2)"
    if "`job'" == "V"  local spec "inertia recip indeg trans same(role) alter(ovz) egoalterint(ovz ovz2)"
    if "`job'" == "R"  local spec "`B7' alter(ovr)"
    if "`job'" == "X"  local nets "ten1x ten3x tenmissx"
    if "`job'" == "B"  local nets "exb1 exb3 exb5"
    if "`job'" == "JP" {
        local nets "ten1 sm270 sm900"
        local act prior
    }
    if "`job'" == "XP" {
        local nets "ten1x ten3x tenmissx"
        local act prior
    }
    if "`job'" == "BP" {
        local nets "exb1 exb3 exb5"
        local act prior
    }
    if "`job'" == "RP" {
        local spec "`B7' alter(ovr)"
        local act prior
    }
    if "`job'" == "JS" {
        local nets "ten1 sm270 sm900"
        local act priorstrict
    }
    if substr("`job'", 1, 2) == "V_" {
        local tag V
        local vr = substr("`job'", 3, .)
    }
    if substr("`job'", 1, 3) == "VP_" {
        local tag VP
        local vr = substr("`job'", 4, .)
        local act prior
    }
    c_local tag "`tag'"
    c_local nets `"`nets'"'
    c_local spec `"`spec'"'
    c_local act "`act'"
    c_local vr "`vr'"
end

* 3d. Run all requested jobs for one league
capture program drop run_league
program define run_league
    args lg
    foreach job of global TAGS {
        job_def `job'
        * variants exist for England only for loose, t5, acc
        if inlist("`vr'", "loose", "t5", "acc") & "`lg'" != "England" continue
        local pre ""
        if "`tag'" == "VP" local pre "prior_"
        if "`vr'" == "" {
            local indir "$OUT"
            local outdir "$OUT"
        }
        else {
            local vname "`vr'"
            if inlist("`vr'", "simple", "base") local vname "`vr'_`lg'"
            local indir "$OUT/variants/`vname'"
            local outdir "$OUT/variants/`pre'`vname'"
            capture mkdir "$OUT/variants/`vname'"
            capture mkdir "`outdir'"
            * inputs of the variant (recipient rule or pass type) from its pass file; made once per league
            if "${built_`vr'_`lg'}" == "" {
                local pf "$BUILD/passes_`vr'.csv"
                if "`vr'" == "base" local pf "$BUILD/passes_all.csv"
                display as text "[`lg'] building inputs of variant `vr' from `pf'"
                fx_export, lg(`lg') passfile("`pf'") outdir("`indir'")
                global built_`vr'_`lg' 1
            }
        }
        local actf "`indir'/act_`lg'.csv"
        if "`act'" == "prior"       local actf "$OUT/act_`lg'_prior.csv"
        if "`act'" == "priorstrict" local actf "$OUT/act_`lg'_priorstrict.csv"
        display as text "[`lg'] job `job' (tag `tag'): " c(current_time)
        fit_one, lg(`lg') tag(`tag') evf("`indir'/ev_`lg'.csv") actf("`actf'") dyf("`indir'/dy_`lg'.csv") ///
            outf("`outdir'/est_`lg'_`tag'.csv") nets("`nets'") spec("`spec'") maxn($MAXN)
        display as text "[`lg'] job `job' done: " c(current_time)
    }
end

* 3e. One Stata batch process per league (shell); each worker runs this same do-file
capture program drop run_parallel
program define run_parallel
    local dofile "$DOFILE"
    local exe "$STATA_EXE"
    if "`exe'" == "" {
        local fl = cond(c(MP), "stata-mp", cond(c(SE), "stata-se", "stata"))
        if strpos(c(machine_type), "Mac") {
            local app = cond(c(MP), "StataMP", cond(c(SE), "StataSE", "StataBE"))
            local exe "`c(sysdir_stata)'`app'.app/Contents/MacOS/`fl'"
        }
        else local exe "`c(sysdir_stata)'`fl'"
        capture confirm file "`exe'"
        if _rc {
            display as error "Stata executable not found (`exe'): set global STATA_EXE or use PARALLEL 0"
            exit 601
        }
    }
    local cmd ""
    foreach lg of global LEAGUES {
        capture mkdir "$W/`lg'"
        local cmd `"`cmd' (cd "$W/`lg'" && env FIT_CHILD=1 FIT_LG="`lg'" FIT_TAGS="$TAGS" FIT_MAXN="$MAXN" FIT_OUT="$OUT" BUILD_DIR="$BUILD" PACK_ROOT="$PACK" NWDIR="$NWDIR" "`exe'" -b do '"`dofile'"' </dev/null >/dev/null 2>&1) &"'
    }
    shell `cmd' wait
end

*------------------------------------------------------------------------------
* 4. Driver
*------------------------------------------------------------------------------
* the do-file's own path for the worker processes
if "$DOFILE" == "" global DOFILE "$PACK/stata/2_fit.do"

local needstatic = 0
if strpos("$STAGES", "inputs") local needstatic = 1
foreach j of global TAGS {
    if substr("`j'", 1, 2) == "V_" | substr("`j'", 1, 3) == "VP_" local needstatic = 1
}
if "`isChild'" == "1" local needstatic = 0

if `needstatic' {
    display as text "== static tables: " c(current_time)
    fx_static
}

if strpos("$STAGES", "inputs") {
    foreach lg of global INPUT_LEAGUES {
        display as text "== inputs `lg': " c(current_time)
        fx_export, lg(`lg') passfile("$BUILD/passes_all.csv") outdir("$OUT") sm
    }
    display as text "== earlier-match positions: " c(current_time)
    fx_prior
    fx_prior, strict
}

if strpos("$STAGES", "fits") {
    nw_setup
    if "`isChild'" != "1" {
        * record the versions of the packages that produce the estimates
        capture {
            file open vf using "$OUT/versions.txt", write text replace
            file write vf "Stata `c(stata_version)' `c(flavor)', date `c(current_date)'" _n
            foreach d in nwcommands nwdynam {
                local vfn "$NWDIR/`d'/VERSION"
                capture confirm file "`vfn'"
                if !_rc {
                    file open vv using "`vfn'", read text
                    file read vv vl
                    file close vv
                    file write vf "`d' `vl'" _n
                }
            }
            file close vf
        }
    }
    local nl : word count $LEAGUES
    if "`isChild'" != "1" & $PARALLEL == 1 & `nl' > 1 {
        display as text "== fits, one process per league: " c(current_time)
        run_parallel
        display as text "== workers finished; their logs are in $W/<league>/2_fit.log: " c(current_time)
    }
    else {
        foreach lg of global LEAGUES {
            run_league `lg'
        }
    }
}
display as text "== 2_fit.do finished: " c(current_time)
