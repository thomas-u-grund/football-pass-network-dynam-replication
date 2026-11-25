*==============================================================================
* 0_install.do   Install nwcommands (core) and nwdynam
*
* Paper: "Fitting dynamic network actor models to thousands of team-matches: pass frequency
*         and receiver choice across match contexts in European football"
*
* The estimates in the paper were made with nwcommands 3.26.0 and nwdynam 1.9.1. This file installs
* the current releases of both. Later releases (nwcommands 3.35.0, nwdynam 1.11.1) gave the same
* estimates in the checks made, apart from the unstable team-matches described in the README.
* If you use this file, leave NWDIR unset when you run 2_fit.do (the packages are then taken from the
* Stata ado path). Requires Stata 17 or later and an internet connection.
*
* Usage:  stata-mp -b do 0_install.do
*         To install into a separate folder instead of your PLUS folder, run first:
*             sysdir set PLUS "/path/to/folder"   and   adopath + "/path/to/folder"
*==============================================================================
version 17
set more off

* 1. core package (the models need it); route as on https://nwcommands.org/nwcommands/install/
net install nwcommands, from("https://nwcommands.org/net/") replace
nwinstall, all nogui

* 2. nwdynam, installed through the core (https://nwcommands.org/nwdynam/install/)
nwinstall, modelpkg(nwdynam)

* 3. check and record the versions
capture which nwdynam
if _rc {
    display as error "nwdynam was not installed; see https://www.nwcommands.org"
    exit 601
}
nwdynam_version
