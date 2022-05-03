#!/bin/bash

if [ "" == "${NGS_HOME}" ]; then
    echo "NGS_HOME must be defined"
    exit -1
fi

if [ "" == "${SPATS_HOME}" ]; then
    echo "SPATS_HOME must be defined"
    exit -1
fi

NGS_MODULE="$1"
if [ "${NGS_MODULE}" == "" ]; then
   NGS_MODULE="ngs"
fi

PYTHONPATH="${NGS_HOME}:${SPATS_HOME}" python3 -m ${NGS_MODULE}
