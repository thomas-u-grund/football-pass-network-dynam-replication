*==============================================================================
* 1_build.do   --   DATA BUILD STAGE of the replication pack
*
* Paper: "Fitting dynamic network actor models to thousands of team-matches: pass frequency and
*         receiver choice across match contexts in European football"
*
* Reads   data/raw/  Wyscout 2017/18 events + matches (7 competitions), players.json, playerank.json
*         data/tm/   Transfermarkt players.csv.gz, transfers.csv.gz
*         data/kg/   sptd/clubs.csv, sptd/transfers.csv, fifa/players_17|18|19.csv
* Writes  output/build/
*   passes_all.csv  passes_simple.csv  passes_acc.csv  passes_loose.csv  passes_t5.csv   passes with inferred receivers
*   player_map.csv                                                                      Wyscout player -> Transfermarkt player
*   nodes.csv  dyads.csv                                                                player-match records; pairs with shared tenure
*   fifa_overall.csv                                                                    FIFA ratings
*   goals.csv                                                                           goals and score state
*   team_cov.csv  match_cov.csv                                                         club covariates, home/away
* Which paper results use which file (see README, Section 5): passes_* feed every pass-network and
* frequency analysis (Tables 1-9, S2-S11, Figure 2); nodes/dyads feed the tenure networks, positions and
* descriptives (Tables 1, 2, 5, 8, S1); fifa_overall the rating covariate (Tables 5, 8); goals the
* goal-state files (Tables 3, 4, S3, S4, S11); team_cov / match_cov the club covariates and home/away of the
* meta-regressions (Tables 3-5, S8).
*
* Usage:   run from the stata/ folder (or set the environment variable PACK_ROOT to the pack folder):
*          stata-mp -b do 1_build.do                        (about 2 minutes)
*          set  global KEEPTMP 1  to keep the intermediate .dta files in output/build/tmp
*
* Requirements: Stata 17 or later (developed and checked with StataNow 19.5); no user-written
* packages. The only external command is `gunzip` (macOS / Linux) to unpack the two .csv.gz files
* into output/build/tmp; on Windows unpack players.csv.gz and transfers.csv.gz there by hand.
*
* Notes on the implementation:
*   - json: the 190 MB event files are cut into one line per event and parsed line by line in Mata
*     (parse_events); the small files are flattened by a generic Mata json reader (jflat).
*   - Fuzzy name matching uses the Ratcliff-Obershelp similarity (sm_matches/sm_ratio; longest
*     block = earliest in a, then in b; no junk characters).
*   - Rounding is exact decimal rounding of the binary value, ties to even (pyround/pround).
*   - Transfer lists are sorted by the string order of date, from-club, to-club, source.
*   - Ties in the club vote (nodes) are broken by the order of the loop; the run prints the number of
*     teams with a tie (0 in these data).
*==============================================================================
*==============================================================================
* 0. SETUP
*==============================================================================
version 17.0
clear all
set more off
set varabbrev off
set linesize 120
timer clear
timer on 99
display as text "1_build.do started: `c(current_date)' `c(current_time)'"

* Root folder of the replication pack: environment variable PACK_ROOT, else default
global root : environment PACK_ROOT
if `"$root"' == "" global root "`c(pwd)'/.."          // run from the stata/ folder, or set PACK_ROOT
global raw  `"$root/data/raw"'
global tm   `"$root/data/tm"'
global kg   `"$root/data/kg"'
global out  `"$root/output/build"'
global tmp  `"$out/tmp"'
capture mkdir `"$root/output"'
capture mkdir `"$out"'
capture mkdir `"$tmp"'
* leagues (the first five are the club leagues)
global LEAGUES England Spain Italy Germany France European_Championship World_Cup
global CLUBLG  England Spain Italy Germany France

* check that the downloaded input files are in place (see README, Section 3)
local need `"`"$raw/players.json"' `"$raw/playerank.json"' `"$tm/players.csv.gz"' `"$tm/transfers.csv.gz"' `"$kg/sptd/clubs.csv"' `"$kg/sptd/transfers.csv"'"'
foreach y in 17 18 19 {
    local need `"`need' `"$kg/fifa/players_`y'.csv"'"'
}
foreach lg of global LEAGUES {
    local need `"`need' `"$raw/events/events_`lg'.json"' `"$raw/matches/matches_`lg'.json"'"'
}
foreach f of local need {
    capture confirm file `"`f'"'
    if _rc {
        display as error `"input file not found: `f'  (see README, Section 3: Data to obtain)"'
        exit 601
    }
}

*==============================================================================
* MATA LIBRARY (all helper logic of this do-file lives here)
*==============================================================================
mata:
mata set matastrict off

// ---- number / string extraction from one JSON text line (flat key search) ----
real scalar jn(string scalar l, string scalar key)
{
    real scalar p, e, q
    string scalar t
    p = strpos(l, key)
    if (p==0) return(.)
    t = substr(l, p+strlen(key), 32)
    e = strpos(t, ",")
    q = strpos(t, "}")
    if (e==0 | (q>0 & q<e)) e = q
    return(strtoreal(substr(t, 1, e-1)))
}
string scalar js(string scalar l, string scalar key)
{
    real scalar p, e
    string scalar t
    p = strpos(l, key)
    if (p==0) return("")
    t = substr(l, p+strlen(key), 80)
    e = strpos(t, `"""')
    return(substr(t, 1, e-1))
}

// ---- Wyscout event file -> Stata dataset in memory (one row per event, file order in idx).
//      The 190 MB single-line json is cut into one line per event (subinstr on the whole text),
//      written to a temporary text file and then read line by line, so memory use stays small.
void parse_events(string scalar fn, string scalar tmpf)
{
    real scalar fh, n, i, p, q
    string scalar s, s1, l, tg, ps
    real matrix X
    string colvector EN, SN, PN
    fh = fopen(fn, "r")
    s = fread(fh, 2000000000)
    fclose(fh)
    s1 = subinstr(s, `"}, {"eventId": "', `"}"'+char(10)+`"{"eventId": "')
    n = strlen(s)-strlen(s1)+1
    s = ""
    unlink(tmpf); fh = fopen(tmpf, "w"); fwrite(fh, s1); fclose(fh)
    s1 = ""
    X = J(n, 13, .); EN = SN = PN = J(n,1,"")
    fh = fopen(tmpf, "r")
    for (i=1; i<=n; i++) {
        l = fget(fh)
        X[i,1] = jn(l, `""matchId": "')
        X[i,2] = jn(l, `""eventSec": "')
        X[i,3] = jn(l, `""teamId": "')
        X[i,4] = jn(l, `""playerId": "')
        EN[i] = js(l, `""eventName": ""')
        SN[i] = js(l, `""subEventName": ""')
        PN[i] = js(l, `""matchPeriod": ""')
        p = strpos(l, `""tags": ["')
        tg = substr(l, p, strpos(substr(l,p,.), "]"))
        X[i,5] = strpos(tg, `"{"id": 1801}"')>0
        X[i,6] = strpos(tg, `"{"id": 101}"')>0
        X[i,7] = strpos(tg, `"{"id": 102}"')>0
        X[i,8] = (strpos(tg, `"{"id": 1701}"')>0 | strpos(tg, `"{"id": 1703}"')>0)
        p = strpos(l, `""positions": ["')
        ps = substr(l, p, strpos(substr(l,p,.), "]"))
        if (strpos(ps, "{")>0) {
            X[i,9]  = jn(ps, `""x": "')
            X[i,10] = jn(ps, `""y": "')
            q = strpos(ps, `"}, {"')
            if (q>0) {
                ps = substr(ps, q+3, .)
                X[i,11] = jn(ps, `""x": "')
                X[i,12] = jn(ps, `""y": "')
            }
        }
        X[i,13] = i
    }
    fclose(fh)
    unlink(tmpf)
    stata("clear")
    st_addobs(n)
    st_store(., st_addvar("double", ("match","sec","team","player","f1801","f101","f102","fred","x0","y0","x1","y1","idx")), X)
    st_sstore(., st_addvar("str30","ename"), EN)
    st_sstore(., st_addvar("str40","sname"), SN)
    st_sstore(., st_addvar("str2","period"), PN)
}

// ---- generic JSON flattener: writes a dataset (rec, path, val) with one row per leaf value.
//      rec = index of the element of the top-level array; path = keys joined by "."
//      (array positions are 0,1,2...). Used for players, matches and playerank files (small files).
string scalar jdecode(string scalar s)
{
    string scalar bs
    bs = char(92)
    if (strpos(s, bs)==0) return(s)
    s = subinstr(s, bs+bs, char(1))
    s = subinstr(s, bs+char(34), char(34))
    s = subinstr(s, bs+"/", "/")
    s = subinstr(s, bs+"t", char(9))
    s = subinstr(s, bs+"n", char(10))
    s = subinstr(s, char(1), bs)
    return(s)
}
void jflat(string scalar fn)
{
    real scalar fh, n, pos, c, depth, needkey, nv, cap, st0, j
    real rowvector b, kind, aidx
    string scalar s, v, pth
    string rowvector comp
    real colvector R
    string colvector P, V
    fh = fopen(fn, "r"); s = fread(fh, 2000000000); fclose(fh)
    b = ascii(s); n = cols(b)
    comp = J(1, 64, ""); kind = J(1,64,0); aidx = J(1,64,0)
    depth = 0; needkey = 0; pos = 1; nv = 0; cap = 500000
    R = J(cap,1,.); P = J(cap,1,""); V = J(cap,1,"")
    while (pos <= n) {
        c = b[pos]
        if (c==32 | c==10 | c==13 | c==9 | c==58) {
            pos = pos + 1
            continue
        }
        if (c==123) {
            depth = depth + 1; kind[depth] = 1; needkey = 1; pos = pos + 1
            continue
        }
        if (c==91) {
            depth = depth + 1; kind[depth] = 2; aidx[depth] = 0; comp[depth] = "0"; needkey = 0; pos = pos + 1
            continue
        }
        if (c==125 | c==93) {
            depth = depth - 1; needkey = 0; pos = pos + 1
            continue
        }
        if (c==44) {
            if (kind[depth]==1) {
                needkey = 1
            }
            else {
                aidx[depth] = aidx[depth] + 1
                comp[depth] = strofreal(aidx[depth], "%12.0f")
            }
            pos = pos + 1
            continue
        }
        if (c==34) {              // string
            st0 = pos + 1; pos = st0
            while (b[pos]!=34) {
                if (b[pos]==92) pos = pos + 1
                pos = pos + 1
            }
            v = jdecode(substr(s, st0, pos-st0)); pos = pos + 1
            if (kind[depth]==1 & needkey) {
                comp[depth] = v; needkey = 0
                continue
            }
        }
        else {                    // number / true / false / null
            st0 = pos
            while (pos<=n) {
                c = b[pos]
                if (c==44 | c==125 | c==93 | c==32 | c==10 | c==13) break
                pos = pos + 1
            }
            v = substr(s, st0, pos-st0)
        }
        // emit leaf
        nv = nv + 1
        if (nv > cap) {
            R = R \ J(cap,1,.); P = P \ J(cap,1,""); V = V \ J(cap,1,"")
            cap = 2*cap
        }
        pth = comp[2]
        for (j=3; j<=depth; j++) pth = pth + "." + comp[j]
        R[nv] = aidx[1]; P[nv] = pth; V[nv] = v
    }
    stata("clear")
    st_addobs(nv)
    st_store(., st_addvar("long","rec"), R[|1\nv|])
    st_sstore(., st_addvar("str244","path"), P[|1\nv|])
    st_sstore(., st_addvar("str244","val"), V[|1\nv|])
}
// ---- round(x, d): exact decimal rounding of the binary value (ties to even), using an
//      error-free product so that values like 3722.0915 round the same way everywhere.
real colvector pround(real colvector x, real scalar d)
{
    real colvector p, c, ahi, alo, e, base, r, sg
    real scalar m, bhi, blo, cc
    m = 10^d
    p = x*m
    c = 134217729*x; ahi = c-(c-x); alo = x-ahi
    cc = 134217729*m; bhi = cc-(cc-m); blo = m-bhi
    e = (((ahi*bhi):-p) :+ ahi*blo :+ alo*bhi) :+ alo*blo
    base = floor(p)
    sg = ((p-base):-0.5) :+ e
    r = base :+ (sg:>0) :+ (sg:==0):*mod(base,2)
    r = r/m
    r[selectindex(x:>=.)] = J(rows(selectindex(x:>=.)),1,.)
    return(r)
}


// ============================ name matching helpers ============================
// sorted unique tokens of a normalised name string (letters a-z and blanks only)
string rowvector toks(string scalar s)
{
    string rowvector t
    t = tokens(s)
    t = select(t, t:!="")
    if (cols(t)==0) return(J(1,0,""))
    return(uniqrows(t')')
}
real scalar nint(string rowvector a, string rowvector b)
{
    real scalar i, j, n
    n = 0
    for (i=1; i<=cols(a); i++) {
        for (j=1; j<=cols(b); j++) {
            if (a[i]==b[j]) { n = n + 1; break; }
        }
    }
    return(n)
}
// total size of the matching blocks (Ratcliff-Obershelp) of strings a and b
// (Ratcliff-Obershelp; longest match = earliest in a, then earliest in b; no junk, strings < 200 chars)
real scalar sm_matches(string scalar sa, string scalar sb)
{
    real rowvector a, b, prev, cur, S
    real scalar na, nb, top, alo, ahi, blo, bhi, i, j, k, besti, bestj, bestsize, total
    a = ascii(sa); b = ascii(sb); na = cols(a); nb = cols(b)
    if (na==0 | nb==0) return(0)
    S = J(1, 400, .); top = 4; S[1..4] = (1, na+1, 1, nb+1)    // stack of half-open ranges [alo,ahi) x [blo,bhi)
    total = 0
    while (top > 0) {
        bhi = S[top]; blo = S[top-1]; ahi = S[top-2]; alo = S[top-3]; top = top - 4
        besti = alo; bestj = blo; bestsize = 0
        prev = J(1, nb+2, 0)
        for (i=alo; i<ahi; i++) {
            cur = J(1, nb+2, 0)
            for (j=blo; j<bhi; j++) {
                if (a[i]==b[j]) {
                    k = prev[j] + 1          // prev[j] holds the length ending at j-1 (index shifted by one)
                    cur[j+1] = k
                    if (k > bestsize) { besti = i-k+1; bestj = j-k+1; bestsize = k; }
                }
            }
            prev = cur
        }
        if (bestsize > 0) {
            total = total + bestsize
            if (alo < besti & blo < bestj) {
                S[top+1] = alo; S[top+2] = besti; S[top+3] = blo; S[top+4] = bestj; top = top + 4
            }
            if (besti+bestsize < ahi & bestj+bestsize < bhi) {
                S[top+1] = besti+bestsize; S[top+2] = ahi; S[top+3] = bestj+bestsize; S[top+4] = bhi; top = top + 4
            }
        }
    }
    return(total)
}
real scalar sm_ratio(string scalar sa, string scalar sb)
{
    real scalar t
    t = strlen(sa) + strlen(sb)
    if (t==0) return(1)
    return(2.0*sm_matches(sa, sb)/t)
}

// ---- wyscout player -> Transfermarkt player. Result written to variables
//      of the dataset in memory: tm_id (blank = none) and how.
void match_players(real colvector wid, string colvector wfirst, string colvector wlast, string colvector wshort,
                   string colvector wdob, real colvector tid, string colvector tname, string colvector tdob)
{
    real scalar np, i, r, c, n, ov, lastov, ok, jb
    transmorphic scalar A
    real rowvector cand, ovs, rat, ord
    string rowvector full, last, tn, tt
    real colvector tmres
    string colvector how
    string scalar lvl, wn
    real rowvector idx, o2
    np = rows(wid)
    // index TM players by date of birth (file order kept) and tokenise their names
    A = asarray_create("string")
    for (r=1; r<=rows(tid); r++) {
        if (tdob[r]!="") {
            if (asarray_contains(A, tdob[r])) asarray(A, tdob[r], (asarray(A, tdob[r]), r))
            else asarray(A, tdob[r], r)
        }
    }
    tmres = J(np,1,.); how = J(np,1,"")
    for (i=1; i<=np; i++) {
        full = toks(wfirst[i] + " " + wlast[i] + " " + wshort[i])
        last = toks(wlast[i])
        lvl = ""
        cand = J(1,0,.)
        if (wdob[i]!="" & asarray_contains(A, wdob[i])) cand = asarray(A, wdob[i])
        n = cols(cand)
        if (n > 0) {
            ovs = J(1, n, 0)
            for (c=1; c<=n; c++) ovs[c] = nint(full, toks(tname[cand[c]]))
            ord = order((-ovs', (1::n)), (1,2))'          // stable sort by descending overlap
            cand = cand[ord]; ovs = ovs[ord]
            ov = ovs[1]
            tn = toks(tname[cand[1]])
            lastov = nint(last, tn)
            ok = (n==1 | (n>1 ? ovs[2] < ov : 0))
            if (ov>=2 & ok) { tmres[i] = tid[cand[1]]; lvl = "dob+name(>=2 tokens)"; }
            else if (ov>=1 & lastov>=1 & ok) { tmres[i] = tid[cand[1]]; lvl = "dob+name(1 token, last)"; }
            else if (ov>=1 & ok) { tmres[i] = tid[cand[1]]; lvl = "dob+one token"; }
            else {
                wn = invtokens(full, "")
                rat = J(1, n, .)
                for (c=1; c<=n; c++) rat[c] = sm_ratio(wn, invtokens(toks(tname[cand[c]]), ""))
                o2 = order((-rat', (1::n)), (1,2))'
                rat = rat[o2]; cand = cand[o2]
                if (rat[1]>=0.8 & (n==1 | (n>1 ? rat[2] < rat[1]-0.1 : 0))) {
                    tmres[i] = tid[cand[1]]; lvl = "dob+fuzzy name"
                }
            }
        }
        how[i] = (lvl=="" ? "unmatched" : lvl)
    }
    st_store(., st_addvar("double","tm_id"), tmres)
    st_sstore(., st_addvar("str30","how"), how)
}

// ---- FIFA row -> Wyscout player: same birth date,
//      name-token overlap, unique best overlap >= 1. Adds variable wid to the dataset in memory.
void fifa_match(real colvector wid, string colvector wtok, string colvector wdob, string colvector ftok, string colvector fdob)
{
    real scalar i, n, c, r
    transmorphic scalar A
    real rowvector cand, ovs, ord
    real colvector res
    string rowvector ft
    A = asarray_create("string")
    for (r=1; r<=rows(wid); r++) {
        if (asarray_contains(A, wdob[r])) asarray(A, wdob[r], (asarray(A, wdob[r]), r))
        else asarray(A, wdob[r], r)
    }
    res = J(rows(ftok),1,.)
    for (i=1; i<=rows(ftok); i++) {
        if (!asarray_contains(A, fdob[i])) continue
        cand = asarray(A, fdob[i]); n = cols(cand)
        ft = toks(ftok[i])
        ovs = J(1,n,0)
        for (c=1; c<=n; c++) ovs[c] = nint(ft, toks(wtok[cand[c]]))
        ord = order((-ovs', (1::n)), (1,2))'
        cand = cand[ord]; ovs = ovs[ord]
        if (ovs[1]>=1 & (n==1 | (n>1 ? ovs[2] < ovs[1] : 0))) res[i] = wid[cand[1]]
    }
    st_store(., st_addvar("double","wid"), res)
}

// ============================ tenure at the club ============================
// merged transfer list of one player: exact Transfermarkt rows + SPTD year rows not already present
// (same to-club within one year), sorted by (date, from, to, source)
string matrix mk_merged(string matrix ex, string matrix sp)
{
    string matrix L
    real scalar i, k, dup
    if (rows(ex)>0) L = (ex, J(rows(ex),1,"tm_exact"))
    else L = J(0,4,"")
    for (i=1; i<=rows(sp); i++) {
        dup = 0
        for (k=1; k<=rows(L); k++) {
            if (sp[i,3]==L[k,3]) {
                if (abs(strtoreal(substr(sp[i,1],1,4)) - strtoreal(substr(L[k,1],1,4))) <= 1) {
                    dup = 1
                    break
                }
            }
        }
        if (!dup) L = L \ (sp[i,.], "sptd_year")
    }
    if (rows(L)>1) L = L[order(L,(1,2,3,4)),]
    return(L)
}
// club c belongs to the organisation of the senior club (same id, or same normalised name up to youth/B-team words)
real scalar same_org(string scalar c, string scalar senior, transmorphic CN)
{
    string scalar a, b
    if (c==senior) return(1)
    a = ""; b = ""
    if (asarray_contains(CN, c)) a = asarray(CN, c)
    if (asarray_contains(CN, senior)) b = asarray(CN, senior)
    if (a=="" | b=="") return(0)
    return(a==b | substr(a, 1, strlen(b)+1)==b+" ")
}
// arrival date into the organisation of `senior` that holds on `date`
//   returns J(1,0,"") = no history, "club_mismatch", or (arrival date, source)
string rowvector org_arrival(string matrix L, string scalar senior, string scalar date, real scalar wid,
                             transmorphic CN, transmorphic FJa, string matrix FJs)
{
    real scalar n, k, i, y, fr
    string scalar d, src, j
    n = 0
    for (k=1; k<=rows(L); k++) {
        if (L[k,1] <= date) n = k
        else break
    }
    if (n==0) return(J(1,0,""))
    if (!same_org(L[n,3], senior, CN)) return("club_mismatch")
    i = n
    while (i>1) {
        if (!same_org(L[i-1,3], senior, CN)) break
        i = i - 1
    }
    d = L[i,1]; src = L[i,4]
    if (src=="sptd_year" & asarray_contains(FJa, wid)) {
        fr = asarray(FJa, wid)
        for (y=1; y<=3; y++) {
            j = FJs[fr,y]
            if (j!="") {
                if (abs(strtoreal(substr(j,1,4)) - strtoreal(substr(d,1,4))) <= 1 & j <= date) {
                    d = j; src = "fifa_exact"
                    break
                }
            }
        }
    }
    return((d, src))
}
real scalar dnum(string scalar s)
{
    return(mdy(strtoreal(substr(s,6,2)), strtoreal(substr(s,9,2)), strtoreal(substr(s,1,4))))
}
// Reads from the dataset in memory (one row per player-team-match of the club leagues):
//   team, date (str10), tid (Transfermarkt id or .), player (Wyscout id)
// and from global Mata objects EXid/EXs (exact transfers), SPid/SPs (SPTD), CNid/CNok (club name keys),
// FJid/FJs (FIFA joined dates). Writes ten_days, ten_src, and prints the club map statistics.
void tenure_calc()
{
    external real colvector EXid, SPid, FJid
    external string matrix EXs, SPs, FJs
    external string colvector CNid, CNok
    real colvector team, tid, player, days, ei, si
    string colvector date, src
    real matrix EXi, SPi
    transmorphic EXa, SPa, CN, FJa, MG, VT, TMAP, loc
    string matrix L, e0, s0
    string rowvector a
    string scalar senior, key, club
    real scalar N, r, k, n, t, nties, nteams, bestc, bests, c, ds
    real rowvector v
    team = st_data(., "team"); tid = st_data(., "tid"); player = st_data(., "player")
    date = st_sdata(., "date"); N = rows(team)
    // lookups
    CN = asarray_create("string")
    for (k=1; k<=rows(CNid); k++) asarray(CN, CNid[k], CNok[k])
    FJa = asarray_create("real")
    for (k=1; k<=rows(FJid); k++) asarray(FJa, FJid[k], k)
    EXi = panelsetup(EXid, 1); EXa = asarray_create("real")
    for (k=1; k<=rows(EXi); k++) asarray(EXa, EXid[EXi[k,1]], k)
    SPi = panelsetup(SPid, 1); SPa = asarray_create("real")
    for (k=1; k<=rows(SPi); k++) asarray(SPa, SPid[SPi[k,1]], k)
    MG = asarray_create("real")
    // ---- 1. vote: senior club of a Wyscout team = most frequent last-known club of its players
    VT = asarray_create("string")
    k = 0
    for (r=1; r<=N; r++) {
        if (tid[r]>=.) continue
        if (!asarray_contains(MG, tid[r])) {
            e0 = J(0,3,""); s0 = J(0,3,"")
            if (asarray_contains(EXa, tid[r])) {
                n = asarray(EXa, tid[r])
                e0 = EXs[|EXi[n,1],1 \ EXi[n,2],3|]
            }
            if (asarray_contains(SPa, tid[r])) {
                n = asarray(SPa, tid[r])
                s0 = SPs[|SPi[n,1],1 \ SPi[n,2],3|]
            }
            asarray(MG, tid[r], mk_merged(e0, s0))
        }
        L = asarray(MG, tid[r])
        n = 0
        for (c=1; c<=rows(L); c++) {
            if (L[c,1] <= date[r]) n = c
            else break
        }
        if (n>0) {
            key = strofreal(team[r], "%12.0f") + "|" + L[n,3]
            if (asarray_contains(VT, key)) {
                v = asarray(VT, key)
                v[1] = v[1] + 1
                asarray(VT, key, v)
            }
            else {
                k = k + 1
                asarray(VT, key, (1, k))
            }
        }
    }
    TMAP = asarray_create("string")
    for (loc = asarray_first(VT); loc != NULL; loc = asarray_next(VT, loc)) {
        key = asarray_key(VT, loc)
        v = asarray(VT, key)
        t = strpos(key, "|")
        club = substr(key, 1, t-1)
        if (asarray_contains(TMAP, club)) {
            a = asarray(TMAP, club)           // (best count, best seq, club id, tie flag)
            bestc = strtoreal(a[1]); bests = strtoreal(a[2])
            if (v[1] > bestc | (v[1]==bestc & v[2] < bests)) {
                asarray(TMAP, club, (strofreal(v[1]), strofreal(v[2]), substr(key,t+1,.), (v[1]==bestc ? "1" : "0")))
            }
            else if (v[1]==bestc) {
                asarray(TMAP, club, (a[1], a[2], a[3], "1"))
            }
        }
        else asarray(TMAP, club, (strofreal(v[1]), strofreal(v[2]), substr(key,t+1,.), "0"))
    }
    nties = 0; nteams = 0
    for (loc = asarray_first(TMAP); loc != NULL; loc = asarray_next(TMAP, loc)) {
        a = asarray(TMAP, asarray_key(TMAP, loc)); nteams = nteams + 1; nties = nties + (a[4]=="1")
    }
    printf("club map: %g teams, %g with a tie for the top club\n", nteams, nties)
    // ---- 2. tenure of every player-team-match
    days = J(N,1,.); src = J(N,1,"")
    for (r=1; r<=N; r++) {
        senior = ""
        key = strofreal(team[r], "%12.0f")
        if (asarray_contains(TMAP, key)) {
            a = asarray(TMAP, key)
            senior = a[3]
        }
        if (tid[r]>=.) {
            src[r] = "none"
            continue
        }
        a = J(1,0,"")
        if (senior!="") {
            a = org_arrival(asarray(MG, tid[r]), senior, date[r], player[r], CN, FJa, FJs)
        }
        if (cols(a)==2) {
            ds = dnum(date[r]) - dnum(a[1])
            days[r] = (ds >= 0 ? ds : 0)
            src[r] = a[2]
        }
        else src[r] = (cols(a)==1 ? "club_mismatch" : "no_history")
    }
    st_store(., st_addvar("double","ten_days"), days)
    st_sstore(., st_addvar("str14","ten_src"), src)
}

end

* rounding wrapper: pyround newvar = oldvar, d(#)
capture program drop pyround
program define pyround
    syntax anything(equalok), d(integer)
    gettoken new rest : anything, parse("=")
    gettoken eq old : rest, parse("=")
    local new = trim("`new'")
    local old = trim("`old'")
    quietly gen double `new' = .
    mata: st_store(., "`new'", pround(st_data(., "`old'"), `d'))
end

* ---- name normalisation: unescape, NFKD, strip accents, lower case, only a-z
capture program drop normvar
program define normvar
    args var
    quietly {
        replace `var' = ustrunescape(`var') if strpos(`var', "\u")
        replace `var' = ustrnormalize(`var', "nfkd")
        replace `var' = ustrregexra(`var', "\p{M}", "")
        replace `var' = ustrlower(`var')
        replace `var' = subinstr(`var', "ø", "o", .)
        replace `var' = subinstr(`var', "đ", "d", .)
        replace `var' = subinstr(`var', "ł", "l", .)
        replace `var' = subinstr(`var', "ß", "ss", .)
        replace `var' = ustrregexra(`var', "[^a-z ]", " ")
        replace `var' = itrim(trim(`var'))
    }
end


*==============================================================================
* 1. EVENTS: parse the Wyscout event json of the seven competitions into one dataset
*
*==============================================================================
local i 0
tempfile evall
foreach lg of global LEAGUES {
    local ++i
    display as text "reading events_`lg'.json ..."
    mata: parse_events(`"$raw/events/events_`lg'.json"', `"$tmp/events_line.txt"')
    gen byte lgn = `i'
    gen str21 league = "`lg'"
    compress
    if `i' > 1 append using `evall'
    save `evall', replace
}
preserve
    clear
    set obs 7
    gen byte lgn = _n
    gen str21 league = ""
    local i 0
    foreach lg of global LEAGUES {
        local ++i
        replace league = "`lg'" in `i'
    }
    save `"$tmp/lgmap.dta"', replace
restore
use `evall', clear
* period index (1H=0, 2H=1, E1=2, E2=3, P=4, other=9)
gen byte per = cond(period=="1H",0,cond(period=="2H",1,cond(period=="E1",2,cond(period=="E2",3,cond(period=="P",4,9)))))
* order of first appearance of a match in its league file
bysort lgn match: egen long mfirst = min(idx)
* continuous match clock: each period starts after the previous one ended
preserve
    collapse (max) mx=sec, by(lgn match per)
    replace mx = max(mx, 0)
    sort lgn match per
    by lgn match: gen double off = 0 if _n==1
    by lgn match: replace off = off[_n-1] + (mx[_n-1] + 1) if _n>1
    by lgn match: gen double mend = off + (mx + 1) if _n==_N     // total length
    by lgn match: egen double mend2 = max(mend)
    drop mend
    rename mend2 mend
    save `"$tmp/off.dta"', replace
restore
merge m:1 lgn match per using `"$tmp/off.dta"', keepusing(off) nogenerate
sort lgn mfirst per sec idx
save `"$tmp/ev.dta"', replace

*==============================================================================
* 2. PASSES with inferred recipient  -> passes_all.csv
*    Rule: the recipient of a pass is the player of the very next event of the same match and
*    period, if that event belongs to the same team, to a different player, and playerId != 0.
*    Feeds all pass-network analyses (Tables 1-8 of the paper, via the export step).
*==============================================================================
use `"$tmp/ev.dta"', clear
* events are sorted by (league, first appearance of match, period, eventSec, file order)
by lgn mfirst: gen long k = _n - 1
* next event (n+1) and the one after (n+2) inside the same match and period
gen byte hn  = (_n < _N) & lgn[_n+1]==lgn & match[_n+1]==match & per[_n+1]==per
gen byte hn2 = (_n < _N-1) & lgn[_n+2]==lgn & match[_n+2]==match & per[_n+2]==per
gen byte same = team[_n+1]==team if hn
gen double rec = player[_n+1] if hn & same==1 & player[_n+1]!=player & player[_n+1]!=0
gen str30 next_event = ename[_n+1] if hn
* variant rules: 'loose' = skip one duel, 't5' = receiving event within 5 s
gen double rec_t5 = rec if hn & (sec[_n+1]-sec) <= 5
gen double rec_loose = rec
replace rec_loose = player[_n+2] if missing(rec) & hn & ename[_n+1]=="Duel" & hn2 ///
        & team[_n+2]==team & player[_n+2]!=player & player[_n+2]!=0
keep if ename=="Pass"
gen double t0 = off + sec
pyround t = t0, d(3)
drop t0
rename (team player period sname f1801 same k) (team_ passer period_ subevent accurate next_team_same seq)
rename team_ team
rename period_ period
rename rec receiver
keep league lgn match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same rec_t5 rec_loose
order league match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same
compress
save `"$tmp/passes_all.dta"', replace
export delimited league match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same ///
    using `"$out/passes_all.csv"', replace
* ---- passes_simple.csv: only 'Simple pass'
preserve
    keep if subevent=="Simple pass"
    export delimited league match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same ///
        using `"$out/passes_simple.csv"', replace
restore
* ---- passes_acc.csv: England, recipient missing for inaccurate passes
preserve
    keep if league=="England"
    replace receiver = . if accurate!=1
    export delimited league match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same ///
        using `"$out/passes_acc.csv"', replace
restore
* ---- passes_loose.csv / passes_t5.csv: England recipient-rule variants
preserve
    keep if league=="England"
    replace receiver = rec_loose
    export delimited league match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same ///
        using `"$out/passes_loose.csv"', replace
restore
preserve
    keep if league=="England"
    replace receiver = rec_t5
    export delimited league match team seq t period passer receiver accurate subevent x0 y0 x1 y1 next_event next_team_same ///
        using `"$out/passes_t5.csv"', replace
restore

*==============================================================================
* 3. WYSCOUT PLAYERS -> TRANSFERMARKT PLAYERS      -> player_map.csv
*    Same date of birth + overlap of name tokens; fuzzy fallback (see Mata block).
*==============================================================================
* ---- Wyscout players.json
mata: jflat(`"$raw/players.json"')
keep if inlist(path,"wyId","firstName","middleName","lastName","shortName") | inlist(path,"birthDate","height","role.code2","foot","birthArea.name")
replace path = subinstr(path, ".", "_", .)
reshape wide val, i(rec) j(path) string
rename val* *
rename (role_code2 birthArea_name) (role birtharea)
foreach v in firstName middleName lastName shortName birthDate foot role birtharea {
    capture confirm variable `v'
    if _rc gen str1 `v' = ""
}
destring wyId height, replace
sort rec
gen long wyorder = rec
* normalised name parts
gen strL nfirst0 = firstName
gen strL nfirst = firstName
gen strL nmid = middleName
gen strL nlast = lastName
gen strL nshort = shortName
foreach v in nfirst0 nfirst nmid nlast nshort {
    normvar `v'
}
replace nfirst = itrim(trim(nfirst + " " + nmid))
save `"$tmp/wy.dta"', replace

* ---- Transfermarkt players.csv.gz (gunzip needed once; temporary plain csv in output/build/tmp)
shell gunzip -c "$tm/players.csv.gz" > "$tmp/tm_players.csv"
import delimited using `"$tmp/tm_players.csv"', clear stringcols(_all) varnames(1) encoding(utf-8) bindquote(strict) maxquotedrows(unlimited)
keep player_id name date_of_birth
gen str10 dob = substr(date_of_birth, 1, 10)
destring player_id, replace
gen strL nname = name
normvar nname
save `"$tmp/tm_players.dta"', replace
erase `"$tmp/tm_players.csv"'

* ---- matching
use `"$tmp/tm_players.dta"', clear
mata: TMid = st_data(., "player_id"); TMname = st_sdata(., "nname"); TMdob = st_sdata(., "dob")
use `"$tmp/wy.dta"', clear
mata: match_players(st_data(., "wyId"), st_sdata(., "nfirst"), st_sdata(., "nlast"), st_sdata(., "nshort"), st_sdata(., "birthDate"), TMid, TMname, TMdob)
keep wyId tm_id how wyorder
rename wyId wy_id
gen double wyId = wy_id          // key for later merges
sort wyorder
save `"$tmp/player_map.dta"', replace
export delimited wy_id tm_id how using `"$out/player_map.csv"', replace
tabulate how

*==============================================================================
* 4. MATCH FILES, PLAYERANK, FIFA, TRANSFERS (inputs of the nodes / dyads / covariate steps)
*==============================================================================
* ---- matches (all seven competitions): match table, teams, line-ups, substitutions
local i 0
foreach lg of global LEAGUES {
    local ++i
    mata: jflat(`"$raw/matches/matches_`lg'.json"')
    gen str21 league = "`lg'"
    if `i' > 1 append using `"$tmp/matches_flat.dta"'
    save `"$tmp/matches_flat.dta"', replace
}
* match level
keep if inlist(path,"wyId","dateutc","seasonId") | regexm(path, "^teamsData\.[0-9]+\.(side|score)$") ///
    | regexm(path, "^teamsData\.[0-9]+\.formation\.(lineup\.[0-9]+\.playerId|substitutions\.[0-9]+\.(playerIn|playerOut|minute))$")
save `"$tmp/matches_flat.dta"', replace
keep if inlist(path,"wyId","dateutc","seasonId")
reshape wide val, i(league rec) j(path) string
gen double match = real(valwyId)
gen str10 date = substr(valdateutc, 1, 10)
keep league rec match date
save `"$tmp/mt.dta"', replace            // league rec match date
* teams: side and score
use `"$tmp/matches_flat.dta"', clear
keep if regexm(path, "^teamsData\.[0-9]+\.(side|score)$")
gen double team = real(ustrregexra(path, "^teamsData\.([0-9]+)\..*$", "$1"))
gen str5 fld = ustrregexra(path, "^teamsData\.[0-9]+\.", "")
keep league rec team fld val
reshape wide val, i(league rec team) j(fld) string
rename (valside valscore) (side score)
destring score, replace
merge m:1 league rec using `"$tmp/mt.dta"', assert(match) nogenerate
save `"$tmp/mteam.dta"', replace         // league rec team side score match date
* line-ups
use `"$tmp/matches_flat.dta"', clear
keep if regexm(path, "\.formation\.lineup\.[0-9]+\.playerId$")
gen double team = real(ustrregexra(path, "^teamsData\.([0-9]+)\..*$", "$1"))
gen double player = real(val)
keep league rec team player
merge m:1 league rec using `"$tmp/mt.dta"', assert(match) keepusing(match) nogenerate
save `"$tmp/lineup.dta"', replace
* substitutions (in order of appearance; later entries override earlier ones for the same player)
use `"$tmp/matches_flat.dta"', clear
keep if regexm(path, "\.formation\.substitutions\.[0-9]+\.(playerIn|playerOut|minute)$")
gen double team = real(ustrregexra(path, "^teamsData\.([0-9]+)\..*$", "$1"))
gen long sub = real(ustrregexra(path, "^.*substitutions\.([0-9]+)\..*$", "$1"))
gen str9 fld = ustrregexra(path, "^.*\.", "")
keep league rec team sub fld val
reshape wide val, i(league rec team sub) j(fld) string
rename (valplayerIn valplayerOut valminute) (pin pout minute)
destring pin pout minute, replace
merge m:1 league rec using `"$tmp/mt.dta"', assert(match) keepusing(match) nogenerate
save `"$tmp/subs.dta"', replace
erase `"$tmp/matches_flat.dta"'

* ---- playerank (only used for the 'prior_rank' node variable)
mata: jflat(`"$raw/playerank.json"')
keep if inlist(path,"matchId","playerId","playerankScore")
reshape wide val, i(rec) j(path) string
gen double match = real(valmatchId)
gen double player = real(valplayerId)
gen double score = real(valplayerankScore)
keep rec match player score
sort rec
save `"$tmp/playerank.dta"', replace

* ---- FIFA ratings 17/18/19: rows, matched to Wyscout players by date of birth and name tokens
local f 0
foreach y in 17 18 19 {
    import delimited using `"$kg/fifa/players_`y'.csv"', clear stringcols(_all) varnames(1) encoding(utf-8)
    keep short_name long_name dob overall potential value_eur wage_eur joined
    gen byte ed = `y'
    gen long rowno = _n
    if `f' > 0 append using `"$tmp/fifa_all.dta"'
    save `"$tmp/fifa_all.dta"', replace
    local ++f
}
gen strL ftok = short_name
gen strL ftok2 = long_name
normvar ftok
normvar ftok2
replace ftok = itrim(trim(ftok + " " + ftok2))
drop ftok2
save `"$tmp/fifa_all.dta"', replace
use `"$tmp/wy.dta"', clear
gen strL wtok = itrim(trim(nfirst0 + " " + nlast + " " + nshort))
mata: WYid = st_data(., "wyId"); WYtok = st_sdata(., "wtok"); WYdob = st_sdata(., "birthDate")
use `"$tmp/fifa_all.dta"', clear
mata: fifa_match(WYid, WYtok, WYdob, st_sdata(., "ftok"), st_sdata(., "dob"))
save `"$tmp/fifa_all.dta"', replace      // + wid = matched Wyscout id (missing if none)

*==============================================================================
* 5. FIFA OVERALL RATINGS  -> fifa_overall.csv
*    Rating/potential of every Wyscout player found in the FIFA 17/18/19 files (later rows overwrite).
*==============================================================================
use `"$tmp/fifa_all.dta"', clear
keep if !missing(wid)
destring overall potential, replace
sort wid ed rowno
by wid ed: keep if _n==_N
keep wid ed overall potential
reshape wide overall potential, i(wid) j(ed)
rename (wid overall17 overall18 overall19 potential18) (wy_id ov17 ov18 ov19 pot18)
keep wy_id ov17 ov18 ov19 pot18
export delimited wy_id ov17 ov18 ov19 pot18 using `"$out/fifa_overall.csv"', replace

*==============================================================================
* 6. GOALS (club leagues)  -> goals.csv
*    A goal = event with tag 101 of type Shot/Free Kick (scoring team) or tag 102 = own goal (other team scores).
*    Time on the continuous match clock (rounded to 0.01 s). Feeds the goal-state files of the models.
*==============================================================================
use `"$tmp/mteam.dta"', clear
keep league match team
sort league match team
by league match: gen byte j = _n
reshape wide team, i(league match) j(j)
rename (team1 team2) (tm1 tm2)
tempfile tpair
save `tpair'
use `"$tmp/ev.dta"', clear
keep if lgn <= 5 & (f102==1 | (f101==1 & inlist(ename,"Shot","Free Kick")))
merge m:1 league match using `tpair', keep(master match) nogenerate
gen double team_scoring = team
replace team_scoring = cond(team==tm1, tm2, tm1) if f102==1
gen double t0 = off + sec
pyround t = t0, d(2)
keep league lgn mfirst per sec idx match team_scoring t
sort lgn mfirst per sec idx
* number of derived goals vs official score
preserve
    gen byte one = 1
    collapse (sum) goals=one, by(league match team_scoring)
    rename team_scoring team
    merge 1:1 league match team using `"$tmp/mteam.dta"', keepusing(score) keep(using match) nogenerate
    replace goals = 0 if missing(goals)
    keep if inlist(league,"England","Spain","Italy","Germany","France")
    gen byte bad = goals != score
    bysort league match: egen byte mbad = max(bad)
    egen tagm = tag(league match)
    quietly count if tagm & mbad
    display as text "goals.csv: matches with goal mismatch (derived vs official score) = " r(N)
restore
export delimited league match team_scoring t using `"$out/goals.csv"', replace

*==============================================================================
* 7. TENURE INPUTS: transfer histories, club-name keys, FIFA joined dates
*==============================================================================
* matched Transfermarkt ids
use `"$tmp/player_map.dta"', clear
keep if !missing(tm_id)
keep tm_id
rename tm_id tid
duplicates drop
save `"$tmp/matched_tid.dta"', replace

* ---- exact Transfermarkt transfers (dcaribou), file order kept
shell gunzip -c "$tm/transfers.csv.gz" > "$tmp/tm_transfers.csv"
import delimited using `"$tmp/tm_transfers.csv"', clear stringcols(_all) varnames(1) encoding(utf-8) bindquote(strict) maxquotedrows(unlimited)
erase `"$tmp/tm_transfers.csv"'
gen long ord = _n
keep ord player_id transfer_date from_club_id to_club_id from_club_name to_club_name
save `"$tmp/tmtr_all.dta"', replace
* club names: SPTD clubs first, then Transfermarkt names (first occurrence; to-club before from-club in a row)
keep ord to_club_id to_club_name from_club_id from_club_name
gen long k1 = 2*ord
gen long k2 = 2*ord + 1
tempfile a1 a2
preserve
    keep k1 to_club_id to_club_name
    rename (k1 to_club_id to_club_name) (key id name)
    save `a1'
restore
keep k2 from_club_id from_club_name
rename (k2 from_club_id from_club_name) (key id name)
append using `a1'
sort key
by key: keep if _n==1
bysort id (key): keep if _n==1
gen byte prio = 2
save `"$tmp/cname_tm.dta"', replace
import delimited using `"$kg/sptd/clubs.csv"', clear stringcols(_all) varnames(1) encoding(utf-8) bindquote(strict) maxquotedrows(unlimited)
gen long key = _n
rename club_name name
keep id name key
bysort id (key): keep if _n==_N          // dict comprehension: last row wins
gen byte prio = 1
append using `"$tmp/cname_tm.dta"'
sort id prio key
by id: keep if _n==1
* organisation key of a club name (okey): normalised tokens without youth / B-team words
gen strL okey = name
normvar okey
replace okey = subinstr(okey, "yth ", "", .)
replace okey = subinstr(okey, " yth", "", .)
replace okey = ustrregexra(okey, "\b(u ?[0-9]{2}|ii|iii|b|yth|youth|juniors?|jugend|reserves?|primavera|jong|sub ?[0-9]+|academy|amateure|castilla|a|2|3|esp|olym)\b", " ")
replace okey = itrim(trim(okey))
keep id okey
save `"$tmp/cname.dta"', replace

* ---- exact transfers of matched players
use `"$tmp/tmtr_all.dta"', clear
destring player_id, gen(tid)
merge m:1 tid using `"$tmp/matched_tid.dta"', keep(match) nogenerate
rename (transfer_date from_club_id to_club_id) (d f t)
keep tid ord d f t
sort tid ord
save `"$tmp/ex.dta"', replace
* ---- SPTD transfers (year only -> 30 June of that year) of matched players
import delimited using `"$kg/sptd/transfers.csv"', clear stringcols(_all) varnames(1) encoding(utf-8) bindquote(strict) maxquotedrows(unlimited)
gen long ord = _n
destring player_id, gen(tid)
merge m:1 tid using `"$tmp/matched_tid.dta"', keep(match) nogenerate
gen str10 d = year + "-06-30"
rename (left_club_id joined_club_id) (f t)
keep tid ord d f t
sort tid ord
save `"$tmp/sp.dta"', replace

* ---- FIFA 'joined' dates by Wyscout id (editions 17, 18, 19; later rows overwrite)
use `"$tmp/fifa_all.dta"', clear
keep if !missing(wid) & joined!=""
sort wid ed rowno
by wid ed: keep if _n==_N
keep wid ed joined
reshape wide joined, i(wid) j(ed)
foreach y in 17 18 19 {
    capture confirm variable joined`y'
    if _rc gen str10 joined`y' = ""
}
sort wid
save `"$tmp/fj.dta"', replace

* load all into Mata
use `"$tmp/ex.dta"', clear
mata: EXid = st_data(., "tid"); EXs = st_sdata(., ("d","f","t"))
use `"$tmp/sp.dta"', clear
mata: SPid = st_data(., "tid"); SPs = st_sdata(., ("d","f","t"))
use `"$tmp/cname.dta"', clear
mata: CNid = st_sdata(., "id"); CNok = st_sdata(., "okey")
use `"$tmp/fj.dta"', clear
mata: FJid = st_data(., "wid"); FJs = st_sdata(., ("joined17","joined18","joined19"))

*==============================================================================
* 8. NODES (player-team-match records)  -> nodes.csv
*    One row per player who appears in the event stream of a club-league match: position/role, age,
*    minutes on the pitch (continuous match clock), tenure at the club with provenance, earlier
*    PlayeRank mean, number of earlier appearances for the team.
*==============================================================================
* match-level helper: start of the 2nd half and total length on the continuous clock
use `"$tmp/off.dta"', clear
gen double off2h = off if per==1
collapse (max) off2h mend, by(lgn match)
save `"$tmp/off2h.dta"', replace

* ---- players who appear in the events of a (match, team); club leagues only
use `"$tmp/ev.dta"', clear
keep if lgn<=5 & player!=0
keep lgn league match team player
duplicates drop
merge m:1 league match using `"$tmp/mt.dta"', keepusing(date) keep(master match)
assert _merge==3
drop _merge
rename player wyId
merge m:1 wyId using `"$tmp/wy.dta"', keepusing(role height foot birthDate birtharea) keep(master match) nogenerate
merge m:1 wyId using `"$tmp/player_map.dta"', keepusing(tm_id) keep(master match) nogenerate
rename (wyId tm_id) (player tid)
sort match team player
* ---- tenure at the club (Mata, see library)
mata: tenure_calc()
tab ten_src
* ---- age at match date
gen double bd = date(birthDate, "YMD")
gen double dd = date(date, "YMD")
gen double age0 = (dd - bd)/365.25 if !missing(bd)
pyround age = age0, d(2)
replace height = . if height==0
* ---- minutes on the pitch: line-up, substitutions, red cards (continuous clock)
merge m:1 lgn match using `"$tmp/off2h.dta"', keepusing(off2h mend) keep(master match) nogenerate
merge 1:1 league match team player using `"$tmp/lineup.dta"', keepusing(rec) keep(master match)
gen byte started = (_merge==3)
drop _merge rec
preserve                                   // substitutions: player in
    use `"$tmp/subs.dta"', clear
    merge m:1 league using `"$tmp/lgmap.dta"', nogenerate
    merge m:1 lgn match using `"$tmp/off2h.dta"', keepusing(off2h mend) keep(match) nogenerate
    gen double tsub = cond(minute<=45, min(minute*60, cond(!missing(off2h), off2h, mend)), ///
                           cond(!missing(off2h), off2h, 0) + (minute-45)*60)
    keep league match team sub pin pout tsub
    tempfile subt
    save `subt'
    rename pin player
    sort league match team player sub
    by league match team player: keep if _n==_N
    keep league match team player tsub
    rename tsub tin_sub
    save `"$tmp/tin_sub.dta"', replace
    use `subt', clear
    rename pout player
    sort league match team player sub
    by league match team player: keep if _n==_N
    keep league match team player tsub
    rename tsub tout_sub
    save `"$tmp/tout_sub.dta"', replace
    use `"$tmp/ev.dta"', clear              // red cards: Foul events with tag 1701/1703 (last one in file order)
    keep if lgn<=5 & fred==1 & ename=="Foul" & player!=0
    gen double tred = off + sec
    sort lgn match team player idx
    by lgn match team player: keep if _n==_N
    keep league match team player tred
    save `"$tmp/tred.dta"', replace
restore
foreach f in tin_sub tout_sub tred {
    merge 1:1 league match team player using `"$tmp/`f'.dta"', keep(master match) nogenerate
}
gen double t_in0 = tin_sub
replace t_in0 = 0 if started==1 & missing(t_in0)
gen double t_out0 = cond(!missing(tout_sub), tout_sub, mend)
replace t_out0 = min(t_out0, tred) if !missing(tred)
pyround t_in = t_in0, d(1)
pyround t_out = t_out0, d(1)
* ---- number of earlier appearances of the player for the team (earlier = strictly earlier date)
sort team player date match
by team player: gen long n_prior_apps = _n - 1
by team player: assert date != date[_n-1] if _n>1
* ---- PlayeRank: mean score and number of earlier matches (any competition in the data)
preserve
    use `"$tmp/playerank.dta"', clear
    merge m:1 match using `"$tmp/mt.dta"', keepusing(date) keep(match) nogenerate
    rename date pdate
    keep player pdate score
    save `"$tmp/pr_dated.dta"', replace
restore
preserve
    keep match team player date
    joinby player using `"$tmp/pr_dated.dta"', unmatched(none)
    keep if pdate < date
    gen long prior_n = 1
    collapse (mean) prior_rank=score (sum) prior_n, by(match team player)
    save `"$tmp/pr_agg.dta"', replace
restore
merge 1:1 match team player using `"$tmp/pr_agg.dta"', keep(master match) nogenerate
replace prior_n = 0 if missing(prior_n)
rename (ten_days ten_src) (tenure_days tenure_src)
save `"$tmp/nodes_all.dta"', replace
order league match team date player role age height foot t_in t_out started tenure_days tenure_src prior_rank prior_n n_prior_apps
export delimited league match team date player role age height foot t_in t_out started tenure_days tenure_src prior_rank prior_n n_prior_apps ///
    using `"$out/nodes.csv"', replace

*==============================================================================
* 9. DYADS (all pairs of players of a team in a match)  -> dyads.csv
*    Shared tenure = shorter of the two tenures (provenance of the weaker source), number of earlier
*    matches in which both appeared (coapp_prior), same role, same birth area.
*==============================================================================
use `"$tmp/nodes_all.dta"', clear
keep league match team date player role birtharea tenure_days tenure_src
save `"$tmp/n_pair.dta"', replace
rename (player role birtharea tenure_days tenure_src) (j role_j area_j days_j src_j)
tempfile nj
save `nj'
use `"$tmp/n_pair.dta"', clear
rename (player role birtharea tenure_days tenure_src) (i role_i area_i days_i src_i)
joinby league match team date using `nj'
keep if i < j
gen double shared_tenure_days = min(days_i, days_j) if !missing(days_i) & !missing(days_j)
gen str14 shared_src = "missing"
replace shared_src = cond(src_i=="tm_exact" & src_j=="tm_exact", "tm_exact", ///
        cond(src_i=="sptd_year" | src_j=="sptd_year", "sptd_year", "fifa_or_mixed")) if !missing(shared_tenure_days)
gen byte same_role = (role_i!="" & role_i==role_j)
gen byte same_birtharea = (area_i!="" & area_i==area_j)
* earlier matches (strictly earlier date) in which both appeared for the team
sort team i j date match
by team i j: gen long coapp_prior = _n - 1
by team i j: assert date != date[_n-1] if _n>1
keep match team i j shared_tenure_days shared_src coapp_prior same_role same_birtharea
export delimited match team i j shared_tenure_days shared_src coapp_prior same_role same_birtharea using `"$out/dyads.csv"', replace

*==============================================================================
* 10. CLUB COVARIATES  -> team_cov.csv, match_cov.csv
*     team_cov: points and goal difference per game (5 leagues), mean log market value / log wage / FIFA
*     overall of the 14 players with most minutes (>= 8 rated), used as club-level covariates in the
*     meta-regressions (Tables 3-5, S8). match_cov: home/away side and opponent of every team-match.
*==============================================================================
* FIFA covariates per Wyscout player: edition 18 preferred, then 17, then 19; rows with unreadable overall skipped
use `"$tmp/fifa_all.dta"', clear
keep if !missing(wid)
gen double ov = real(overall)
keep if !missing(ov)
gen double value = real(value_eur)
gen double wage = real(wage_eur)
replace value = 0 if missing(value)
replace wage = 0 if missing(wage)
gen long prio = cond(ed==18, 3, cond(ed==17, 2, 1))*1000000 + rowno
sort wid prio
by wid: keep if _n==_N
keep wid value wage ov
rename wid player
save `"$tmp/V.dta"', replace
* minutes per team and player
use `"$tmp/nodes_all.dta"', clear
gen double mins = (t_out - cond(missing(t_in), 0, t_in))/60
collapse (sum) mins (firstnm) league, by(team player)
gsort team -mins player
by team: keep if _n <= 14
merge m:1 player using `"$tmp/V.dta"', keep(master match) nogenerate
keep if !missing(value) & value > 0
gen double lv = ln(value)
gen double lw = ln(wage) if wage > 0
collapse (mean) lval=lv lwage=lw ovr=ov (count) nrated=lv (firstnm) league, by(team)
keep if nrated >= 8
tempfile cov
save `cov'
* points / goal difference per game from the match files
use `"$tmp/mteam.dta"', clear
keep if inlist(league,"England","Spain","Italy","Germany","France")
rename (team score side) (team_a score_a side_a)
tempfile ma
save `ma'
keep league rec team_a score_a
rename (team_a score_a) (team_b score_b)
joinby league rec using `ma'
keep if team_a != team_b
gen double pts = cond(score_a > score_b, 3, cond(score_a == score_b, 1, 0))
gen double gd = score_a - score_b
* match_cov.csv
preserve
    keep match team_a side_a team_b
    rename (team_a side_a team_b) (team side opp)
    export delimited match team side opp using `"$out/match_cov.csv"', replace
restore
collapse (sum) pts gd (count) gp=pts, by(team_a)
gen double ppg = pts/gp
gen double gdpg = gd/gp
rename team_a team
merge 1:1 team using `cov', keep(match) nogenerate
keep team league ppg gdpg lval lwage ovr nrated
order team league ppg gdpg lval lwage ovr nrated
export delimited team league ppg gdpg lval lwage ovr nrated using `"$out/team_cov.csv"', replace

* clean up
if `"$KEEPTMP"' != "1" {
    local fl : dir `"$tmp"' files "*"
    foreach f of local fl {
        capture erase `"$tmp/`f'"'
    }
    capture rmdir `"$tmp"'
}
timer off 99
quietly timer list 99
display as text "1_build.do finished; run time " %6.1f r(t99) " seconds"
