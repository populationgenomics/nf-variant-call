#!/usr/bin/env bash
# Print FREEMIX, the estimated contamination fraction, from a VerifyBamID2 .selfSM: the only parser of that file.
# The column is found by its header name, not its position. The value is printed as written, with no cap.
# Fails if the file is missing or empty, has no FREEMIX column, has other than one data row, or if the value
# isn't a number in [0, 1].
#   usage: freemix.sh <sg>.verify-bamid.selfSM
set -euo pipefail

selfsm=${1:?usage: freemix.sh <selfSM>}
if [[ ! -s ${selfsm} ]]; then
    echo "freemix.sh: ${selfsm} is missing or empty" >&2
    exit 1
fi

awk -F '\t' -v file="${selfsm}" '
    function fail(message) {
        print "freemix.sh: " file ": " message > "/dev/stderr"
        failed = 1
        exit 1
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "FREEMIX") {
                column = i
            }
        }
        if (!column) {
            fail("no FREEMIX column in the header")
        }
        next
    }
    NF == 0 { next }
    {
        rows++
        value = $column
    }
    END {
        if (failed) {
            exit 1
        }
        if (rows != 1) {
            fail("expected one data row, found " rows + 0)
        }
        if (value !~ /^([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/ || value + 0 > 1) {
            fail("FREEMIX is not a number in [0, 1]: \"" value "\"")
        }
        print value
    }
' "${selfsm}"
