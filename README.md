# Replication materials: "Fitting dynamic network actor models to thousands of team-matches: pass frequency and receiver choice across match contexts in European football"

The whole analysis is run in Stata, with three do-files. This folder contains the code, the stored estimates and the tables and figures behind the manuscript. It does not contain third-party raw data (Section 3 says where to obtain them).

## 1. What is in this folder

```
README.md, LICENSE
stata/      0_install.do  1_build.do  2_fit.do  3_analyze.do
data/       empty: you put the downloaded third-party data here (Section 3)
output/     created by the do-files: build/, fit/, analyze/
stored/     the estimates and tables of the paper, for comparison with your own run
  estimates/   est_<league>_<tag>.csv (team-match estimates per specification; recipient-rule and pass-type variants in estimates/variants/)
  tables_figures/   one csv per table, Figures 1 and 2
```

- `stata/0_install.do`: installs the Stata packages `nwcommands` and `nwdynam`.
- `stata/1_build.do`: data build. Turns the downloaded data into the analysis files (passes with inferred receivers, player-match records, pairs with shared tenure, FIFA ratings, goals, club covariates).
- `stata/2_fit.do`: model inputs and fits. Fits one dynamic network actor model (`nwdynam`) per team-match for every specification in the paper, one estimate file `est_<league>_<tag>.csv` per league and specification.
- `stata/3_analyze.do`: pooling, regressions, tables and figures: Poisson models of pass frequency, observed shares, pass geometry, one-stage conditional logit, random-effects pooling with club-clustered variance, moderators, Figures 1 and 2.

## 2. Requirements

- Stata 19.5 MP (StataNow); the code uses version 17 syntax. The parallel fit in `2_fit.do` and the unpacking of two `.csv.gz` files use Unix commands (Mac or Linux); on Windows set `global PARALLEL 0` before running `2_fit.do` and unpack `data/tm/players.csv.gz` and `data/tm/transfers.csv.gz` by hand into `output/build/tmp`. The TIF export of the figures uses the macOS tool `sips`; elsewhere the PNG and EPS files are still written.
- The Stata packages `nwcommands` and `nwdynam`, installed by `stata/0_install.do` (`net install nwcommands, from(https://nwcommands.org/net/) replace`, `nwinstall, all`, then nwdynam v1.9.1 from its GitHub release; see https://nwcommands.org/nwdynam/install/). If you install this way, leave `NWDIR` unset. The estimates in the paper were made with nwcommands 3.26.0 and nwdynam 1.9.1; later versions (nwdynam 1.11.1, nwcommands 3.35.0) gave identical estimates for the test case. 
- About 3 GB of disk space for the data and the built files.

## 3. Data to obtain and where to put it

No third-party data are distributed with this pack. You need four sources. All are free, and two need a free login (Kaggle, IEEE DataPort). Put the files exactly as shown below, keeping the file names; `1_build.do` checks that all 23 files are in place and stops with a message naming any missing one. Please read the licence on each source page before use or redistribution.

Create these folders inside `data/` (the starting point is the folder that contains `stata/`):

```
data/raw/events/      data/raw/matches/      data/tm/      data/kg/sptd/      data/kg/fifa/
```

**Step 1. Wyscout match events (needed for everything).**
Source: Pappalardo et al. (2019), "A public data set of spatio-temporal match events in soccer competitions", Scientific Data 6, 236. The data are on Figshare: https://figshare.com/collections/Soccer_match_event_dataset/4415000 (free, no login).
1. Open the collection page and download the files `events.zip`, `matches.zip`, `players.json` and `playerank.json`. (The collection has further files, for example teams and competitions, that this analysis does not use.)
2. Unzip `events.zip`: it contains one file per competition, `events_England.json`, `events_Spain.json`, `events_Italy.json`, `events_Germany.json`, `events_France.json`, `events_European_Championship.json` and `events_World_Cup.json`. Put all seven in `data/raw/events/`.
3. Unzip `matches.zip` in the same way and put the seven `matches_<competition>.json` files in `data/raw/matches/`.
4. Put `players.json` and `playerank.json` directly in `data/raw/`.
The analysis uses only the five national leagues, but the build also reads the two tournament files, so keep all seven.

**Step 2. FIFA player ratings (receiver rating, Tables 5 and 8; squad value, Table S3).**
Source: Kaggle, "FIFA 20 complete player dataset" by user stefanoleone992 (scraped from sofifa.com; licence field CC0-1.0): https://www.kaggle.com/datasets/stefanoleone992/fifa-20-complete-player-dataset
1. Log in to Kaggle (free account) and download the dataset (a zip with one csv per FIFA edition).
2. From it, take `players_17.csv`, `players_18.csv` and `players_19.csv` and put them in `data/kg/fifa/`. The other editions are not used.

**Step 3. Transfermarkt (shared tenure, Table 8 tenure rows, Supplementary Table S1).**
Source: the open project `dcaribou/transfermarkt-datasets`: https://github.com/dcaribou/transfermarkt-datasets (licence listed as CC0-1.0). The repository README links to the download of the data files.
1. Download `players.csv.gz` and `transfers.csv.gz` (keep them compressed).
2. Put both in `data/tm/`. The build unpacks them itself.
The project updates its files over time. The copies used for the paper were downloaded on 18 November 2026; later snapshots can differ slightly, which can move the tenure results in the third decimal.

**Step 4. Sports Transfer and Player Dataset (shared tenure for players Transfermarkt misses).**
Source: IEEE Data Descriptions, doi 10.1109/IEEEDATA.2025.3569490 (the authors state CC BY 4.0). Open the DOI, go to the dataset page on IEEE DataPort (free login) and download `clubs.csv` and `transfers.csv`. Put both in `data/kg/sptd/`.

Final layout:

```
data/raw/events/events_<England|Spain|Italy|Germany|France|European_Championship|World_Cup>.json
data/raw/matches/matches_<the same seven>.json
data/raw/players.json    data/raw/playerank.json
data/tm/players.csv.gz   data/tm/transfers.csv.gz
data/kg/sptd/clubs.csv   data/kg/sptd/transfers.csv
data/kg/fifa/players_17.csv   players_18.csv   players_19.csv
```

Without the FIFA, Transfermarkt and SPTD files, the build stops. Everything about pass frequency, positions and earlier passes depends only on the Wyscout files; the rating and tenure results need the other three.

## 4. How to reproduce

Open a terminal in the `stata/` folder of this pack (the do-files find the pack from there; you can also set the environment variable `PACK_ROOT` to the pack folder) and run:

```
stata-mp -b do 0_install.do    # installs nwcommands and nwdynam (needs internet), a few minutes
stata-mp -b do 1_build.do      # about 2 to 3 minutes   -> output/build/
stata-mp -b do 2_fit.do        # about 30 minutes with the five leagues in parallel   -> output/fit/
stata-mp -b do 3_analyze.do    # about 10 minutes       -> output/analyze/
```

(Alternative to `0_install.do`: install the packages yourself and set `NWDIR` to the folder with `nwcommands/` and `nwdynam/`.) Each stage reads the output of the previous one. The header of each do-file lists switches: `3_analyze.do` takes block names as arguments to run only part of the analysis (`desc rate observed geometry outcome clogit pool moderators figures simplectx`); `2_fit.do` takes the globals `STAGES`, `LEAGUES`, `TAGS`, `MAXN` and `PARALLEL`. No random numbers are used. To check your run, compare `output/analyze/` with `stored/tables_figures/` and `output/fit/` with `stored/estimates/`.

## 5. Which output gives which table

`3_analyze.do` writes one csv per table to `output/analyze/`, named `Table<number>_...csv` and `TableS<number>_...csv`, and Figure1 and Figure2 as PNG, EPS and TIF.

| Output | Content |
|---|---|
| Tables 1, 2 | sample by league; descriptives |
| Tables 3, 4 and S3, S4 | pass frequency (Poisson with player-time exposure, clustered by club) |
| Table 5, Table 8, Table S1, S10 | pooled team-match effects by specification (tags below) |
| Table 6, Figure 2 | observed outcomes by context |
| Table 9, Table S2 | recipient rule and pass type variants |
| Figure 1 | variation of the team-match estimates (tag `JP`) |
| Table S5 (and S5b, S5c) | moderators of the receiver-choice effects; club fixed effects; Benjamini-Hochberg count |
| Table S6 | geometry of passes by clock and score |
| Table S7 | pass frequency and the goal difference |
| Tables 7, S8, S9 | one-stage conditional logit and predicted passing patterns (Table 7 is a subset of S9; takes most of the 10 minutes) |
| Table S11 | context contrasts for simple passes only (block `simplectx`; needs the variant estimates of `2_fit.do`, tags `VP_base` and `VP_simple`; runs in seconds) |

Specification tags: `J` main model with tenure and shared-minutes networks, positions from the match; `JP`, `XP`, `BP`, `RP` the same with earlier-match positions (`JP` is the main specification of the paper); `JS` strict role fallback; `U`, `V` models without position controls; `R` rating standardised within role; `X`, `B` exact-tenure models; `V_*` and `VP_*` recipient-rule and pass-type variants.

## 6. What was checked

- A complete run of the four do-files from the raw data in an empty folder (24 Nov 2026, Stata 19.5, nwdynam 1.11.1, nwcommands distribution of 5 Sep 2026; about 2, 36 and 8 minutes for the three stages) reproduced the stored estimate files and tables: all tables agree to rounding error, except that the exact-tenure specification `XP` differs in the fourth decimal (see the next point). An earlier complete run with nwdynam 1.9.1 and nwcommands 3.26.0 gave the same estimates.
- Unstable team-matches: in the main specification `J`, 31 of 3652 team-matches have coefficients that are not identified (they run off to about 1e7); they are dropped in pooling, and the pooled results are unchanged to four decimals. The exact-tenure specifications `X` and `XP` have 14 to 59 such team-matches per league. For these degenerate team-matches the fit is not reproducible bit for bit from run to run (a repeated fit of `XP` for Germany changed 7 of 612 team-matches), which moves the pooled `XP` estimates by less than 0.001 (Table S1).
- In the moderator regressions (Table S5), the iterated DerSimonian-Laird estimate of the between-match variance does not converge for some terms (it alternates between two values), and the coefficients use the weights of the previous pass. The code documents this at `mreg` in `3_analyze.do`.
- The one-stage conditional logit is written in Mata, because `clogit` would need about 14 million long-format rows.
- No random numbers are used in any step.

## 7. Terms of use

The code and the estimates in this folder accompany the manuscript. No third-party data are redistributed; Section 3 gives the link to each source. The code is released under the licence in `LICENSE`.
