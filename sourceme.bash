#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#
# Environment setup for package spaceranger.
# Source this file before running.
#
# Determine path to this script; resolve symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
#
# Source user's own environment first.
#
# Only source .bashrc if we're being sourced from shell10x script.
# Otherwise, we could end up in an infinite loop if user is
# sourcing this file from their .bashrc.
# Note: .bash_profile is for login shells only.
if [ ! -z $_RUN10X ] && [ -e ~/.bashrc ]; then
    source ~/.bashrc
fi
#
# Modify the prompt to indicate user is in 10X environment.
#
if [ ! -z $_10X_MODPROMPT ]; then
    ORIG_PROMPT=$PS1
    PREFIX=$TENX_PRODUCT
    if [ ! -z $_10XDEV_BRANCH_PROMPT ]; then
        command -v git >/dev/null 2>&1 || { echo >&2 "Error: git is required but not found."; exit 1; }
        PREFIX="10X:\`pushd $(echo $MROPATH | cut -d : -f1) > /dev/null;git rev-parse --abbrev-ref HEAD;popd > /dev/null\`"
    fi
    export PS1="\[\e[0;34m\]$PREFIX\[\e[m\]>$ORIG_PROMPT"
fi
#
# Set aside environment variables if they may conflict with 10X environment
#

if [ -z "$_TENX_LD_LIBRARY_PATH" ]; then
    export _TENX_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH=""
fi

#
# Unset environment variables if they may conflict with 10X environment
#

if [ ! -z "$PYTHONPATH" ]; then
    unset PYTHONPATH
fi

if [ ! -z "$PYTHONHOME" ]; then
    unset PYTHONHOME
fi

if [ ! -z "$MROPATH" ]; then
    unset MROPATH
fi

#
# Add module binary paths to PATH
#

if [ -z "$PATH" ]; then
    export PATH="$DIR/bin"
else
    export PATH="$DIR/bin:$PATH"
fi

if [ -z "$PATH" ]; then
    export PATH="$DIR/external/anaconda/bin"
else
    export PATH="$DIR/external/anaconda/bin:$PATH"
fi

if [ -z "$PATH" ]; then
    export PATH="$DIR/external/martian/bin"
else
    export PATH="$DIR/external/martian/bin:$PATH"
fi

if [ -z "$PATH" ]; then
    export PATH="$DIR/bin/rna"
else
    export PATH="$DIR/bin/rna:$PATH"
fi

if [ -z "$PATH" ]; then
    export PATH="$DIR/lib/bin"
else
    export PATH="$DIR/lib/bin:$PATH"
fi

if [ -z "$PATH" ]; then
    export PATH="$DIR/tenkit/bin"
else
    export PATH="$DIR/tenkit/bin:$PATH"
fi

if [ -z "$PATH" ]; then
    export PATH="$DIR/tenkit/lib/bin"
else
    export PATH="$DIR/tenkit/lib/bin:$PATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/tenkit/lib/python"
else
    export PYTHONPATH="$DIR/tenkit/lib/python:$PYTHONPATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/external/martian/adapters/python"
else
    export PYTHONPATH="$DIR/external/martian/adapters/python:$PYTHONPATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/external/illuminate"
else
    export PYTHONPATH="$DIR/external/illuminate:$PYTHONPATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/lib/python"
else
    export PYTHONPATH="$DIR/lib/python:$PYTHONPATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/external/tsne"
else
    export PYTHONPATH="$DIR/external/tsne:$PYTHONPATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/external/"
else
    export PYTHONPATH="$DIR/external/:$PYTHONPATH"
fi

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$DIR/external/umap"
else
    export PYTHONPATH="$DIR/external/umap:$PYTHONPATH"
fi

if [ -z "$MROPATH" ]; then
    export MROPATH="$DIR/tenkit/mro"
else
    export MROPATH="$DIR/tenkit/mro:$MROPATH"
fi

if [ -z "$MROPATH" ]; then
    export MROPATH="$DIR/mro"
else
    export MROPATH="$DIR/mro:$MROPATH"
fi

#
# Module-specific env vars
#
export LC_ALL="C"
export PYTHONNOUSERSITE="1"

