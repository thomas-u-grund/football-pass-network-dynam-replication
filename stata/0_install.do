*==============================================================================
* 0_install.do   Install nwcommands (core) and nwdynam from GitHub
*
* Paper: "Fitting dynamic network actor models to thousands of team-matches: pass frequency
*         and receiver choice across match contexts in European football"
*
* The estimates in the paper were made with nwcommands 3.26.0 and nwdynam 1.9.1. This file installs
* the current release of nwcommands (nwcommands.org) and the tagged release RELEASE of nwdynam (global NWDYNAM_RELEASE,
* default v1.9.1; if that tag cannot be installed, the current release is installed instead). Later
* releases (nwcommands 3.35.0, nwdynam 1.11.1) gave identical estimates in the checks made.
* If you use this file, leave NWDIR unset when you run 2_fit.do (the packages are then taken from the
* Stata ado path). Requires Stata 17 or later and an internet connection.
*
* Usage:  stata-mp -b do 0_install.do
*         To install into a separate folder instead of your PLUS folder, run first:
*             sysdir set PLUS "/path/to/folder"   and   adopath + "/path/to/folder"
*==============================================================================
version 17
set more off
global NWDYNAM_RELEASE : environment NWDYNAM_RELEASE
if "$NWDYNAM_RELEASE" == "" global NWDYNAM_RELEASE "v1.9.1"

* 1. core package (the models need it); route as on https://nwcommands.org/nwcommands/install/
net install nwcommands, from("https://nwcommands.org/net/") replace
nwinstall, all nogui

* 2. nwdynam: tagged release, else the current release (https://nwcommands.org/nwdynam/install/)
capture noisily net install nwdynam, from("https://raw.githubusercontent.com/thomas-u-grund/nwdynam/$NWDYNAM_RELEASE/") replace
if _rc {
    display as text "release $NWDYNAM_RELEASE not available; installing the current release of nwdynam"
    net install nwdynam, from("https://raw.githubusercontent.com/thomas-u-grund/nwdynam/main/") replace
}

* 3. check and record the versions
capture which nwdynam
if _rc {
    display as error "nwdynam was not installed; see https://www.nwcommands.org"
    exit 601
}
nwdynam_version
