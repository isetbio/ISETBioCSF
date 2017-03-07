#!/bin/bash
logfile="runPoirsonAndWandell96Job3.log"
errorfile="runPoirsonAndWandell96Job3.errorlog"

if [ -f $logfile ] ; then
    rm $logfile
fi

if [ -f $errorfile ] ; then
    rm $errorfile
fi

# NOTE:
# in .bashrc add: alias matlab='/Applications/MATLAB_R2016a.app/bin/matlab -nodisplay'
# or whatever the path to the current MATLAB is

# Run the .bashrc
. ~/.bashrc

# Go !
matlab -nodisplay -noawt -nosplash -nodesktop -r 'try runPoirsonAndWandell96Job(3); catch; end; quit' 1> $logfile 2>$errorfile &
