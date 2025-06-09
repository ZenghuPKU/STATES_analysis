#!/bin/bash

start=1
end=132

user_dir="/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata"
sample="B4_WT_mousebrain2"
subdir="04_stitch"

getOutpath() {
    echo "${user_dir}/${sample}/${subdir}"
}

seq $start $end > "$(getOutpath)/orderlist" 
