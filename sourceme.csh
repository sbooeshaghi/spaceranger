#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#
# Environment setup for package spaceranger.
# Source this file before running.
#
set SOURCE=($_)
if ( "$SOURCE" != "" ) then
    set SOURCE=`readlink -f "$SOURCE[2]"`
else
    set SOURCE=`readlink -f "$0"`
endif
set DIR=`dirname $SOURCE`
#
# Source user's own environment first.
#
# Note .login is for login shells only.
source ~/.cshrc
#
# Modify the prompt to indicate user is in 10X environment.
#

#
# Set aside environment variables if they may conflict with 10X environment
#

if ( ! $?_TENX_LD_LIBRARY_PATH ) then
    setenv _TENX_LD_LIBRARY_PATH "$LD_LIBRARY_PATH"
    setenv LD_LIBRARY_PATH ""
endif

#
# Unset environment variables if they may conflict with 10X environment
#

if ( $?PYTHONPATH ) then
    unsetenv PYTHONPATH
endif

if ( $?PYTHONHOME ) then
    unsetenv PYTHONHOME
endif

if ( $?MROPATH ) then
    unsetenv MROPATH
endif

#
# Add module binary paths to PATH
#

if ( ! $?PATH ) then
    setenv PATH "$DIR/bin"
else
    setenv PATH "$DIR/bin:$PATH"
endif

if ( ! $?PATH ) then
    setenv PATH "$DIR/external/anaconda/bin"
else
    setenv PATH "$DIR/external/anaconda/bin:$PATH"
endif

if ( ! $?PATH ) then
    setenv PATH "$DIR/external/martian/bin"
else
    setenv PATH "$DIR/external/martian/bin:$PATH"
endif

if ( ! $?PATH ) then
    setenv PATH "$DIR/bin/rna"
else
    setenv PATH "$DIR/bin/rna:$PATH"
endif

if ( ! $?PATH ) then
    setenv PATH "$DIR/lib/bin"
else
    setenv PATH "$DIR/lib/bin:$PATH"
endif

if ( ! $?PATH ) then
    setenv PATH "$DIR/tenkit/bin"
else
    setenv PATH "$DIR/tenkit/bin:$PATH"
endif

if ( ! $?PATH ) then
    setenv PATH "$DIR/tenkit/lib/bin"
else
    setenv PATH "$DIR/tenkit/lib/bin:$PATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/tenkit/lib/python"
else
    setenv PYTHONPATH "$DIR/tenkit/lib/python:$PYTHONPATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/external/martian/adapters/python"
else
    setenv PYTHONPATH "$DIR/external/martian/adapters/python:$PYTHONPATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/external/illuminate"
else
    setenv PYTHONPATH "$DIR/external/illuminate:$PYTHONPATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/lib/python"
else
    setenv PYTHONPATH "$DIR/lib/python:$PYTHONPATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/external/tsne"
else
    setenv PYTHONPATH "$DIR/external/tsne:$PYTHONPATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/external/"
else
    setenv PYTHONPATH "$DIR/external/:$PYTHONPATH"
endif

if ( ! $?PYTHONPATH ) then
    setenv PYTHONPATH "$DIR/external/umap"
else
    setenv PYTHONPATH "$DIR/external/umap:$PYTHONPATH"
endif

if ( ! $?MROPATH ) then
    setenv MROPATH "$DIR/tenkit/mro"
else
    setenv MROPATH "$DIR/tenkit/mro:$MROPATH"
endif

if ( ! $?MROPATH ) then
    setenv MROPATH "$DIR/mro"
else
    setenv MROPATH "$DIR/mro:$MROPATH"
endif

#
# Module-specific env vars
#
setenv LC_ALL "C"
setenv PYTHONNOUSERSITE "1"

