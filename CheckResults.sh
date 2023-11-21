#!/bin/bash

file_load=$1

root -l << EOF
TString file = "${file_load}"
.L checkResults.C 
checkResults(file)
EOF
