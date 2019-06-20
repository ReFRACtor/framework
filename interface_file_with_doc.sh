#!/bin/bash

DOXY_CONFIG="$(dirname $0)/doxygen_single.cfg"
DOXY2SWIG="$(dirname $0)/doxy2swig.py"

doxygen_bin=$1

if [ ! -e $doxygen_bin ]; then
    echo "Doxygen binary not found: $doxygen_bin"
    exit 1
fi

out_dir=$2

mkdir -p $out_dir

input_swig_file=$3
output_swig_file=$4

# Find source code to run through Doxygen
source_header=$(echo $input_swig_file | sed 's/\.i$/.h/')
if [ ! -e "$source_header" ]; then
    echo "Could not find header file for $input_swig_file"
    cp $input_swig_file $output_swig_file
else
    OUTDIR=$out_dir INPUT=$source_header $doxygen_bin $DOXY_CONFIG
    $DOXY2SWIG ${out_dir}/xml/index.xml $output_swig_file
    cat $input_swig_file >> $output_swig_file
fi
