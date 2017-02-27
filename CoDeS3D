#!/bin/sh

# Set up a CoDeS3D working environment.  This should either run the CoDeS3D command directly
# or else drop the user to a shell where all CoDeS3D commands are accessible from anywhere in the filesystem.
# Currently only compatible with BaSH.

CODES3DPATH="$(readlink -f $0)"
while [ -h "$CODES3DPATH" ]; do
    CODES3DDIR="$(cd -P "$(dirname "$CODES3DPATH" )" && pwd )"
    CODES3DPATH="$(readlink "$CODES3DPATH")"
    [[ $CODES3DPATH != /* ]] && SOURCE="$CURRDIR/$CODES3DPATH"
done
CODES3DDIR="$(cd -P "$(dirname "$CODES3DPATH" )" && pwd )"

if [ "$*" = "" ] ; then
    # This script may be run as a symlink - need relative path from target if so
    # In this instance, a CoDeS3D-enabled shell should be started
    if [ `basename "$SHELL"` = "bash" ] ; then
        #Set bash dot directory.
        echo "Running $SHELL --rcfile $CODES3DDIR/docs/.bashrc"
        CODES3DSHELL="$SHELL --rcfile $CODES3DDIR/docs/.bashrc"
    else
        if which bash > /dev/null ; then
            echo "Running bash --rcfile $CODES3DDIR/docs/.bashrc"
            CODES3DSHELL="bash --rcfile $CODES3DDIR/docs/.bashrc"
        else
            echo "To start a CoDeS3D shell, you must have the BaSH shell available"
            echo "and in your path."
            exit 1
        fi
    fi

echo """
Setting up CoDeS3D environment.

Type 'exit' or press Ctrl+D at any time to leave.
"""
    # Add codes3d scripts to $PYTHONPATH
    export PATH="$PATH:$CODES3DDIR/codes3d"
    # Run "CoDeS3D shell"
    eval exec $CODES3DSHELL

else
    # Otherwise, the user only wishes to run a single command
    export PATH="$PATH:$CODES3DDIR/codes3d"
    cmd=`basename "$1" .py`.py
    shift

    exec "$CODES3DDIR/codes3d/$cmd" "$@"
fi