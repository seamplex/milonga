#!/bin/bash
# run the 3D IAEA PWR benchmark over a structured mesh of
# characteristic length lc=10cm with the default scheme

. locateruntest.sh

$milongabin ${testdir}3dpwr-struct.mil
outcome=$?

callpandoc 3dpwr-struct
callgmsh 3dpwr-struct.pos

exit $outcome
