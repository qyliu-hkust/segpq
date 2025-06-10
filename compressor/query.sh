#!/bin/bash

#dataset="rb04-4096"
#dataset="cw12b-1M"
#dataset="ccnews-1M"
#dataset="cw12b-4096"

result_dir="/mnt/home/xyzhu/result"
data_dir="/mnt/home/xyzhu/datasets"
code_dir="/mnt/home/xyzhu/codes/Build"

read_only="t"
source_dir="../pgm_index"
index_type="pgm"

mkdir -p "$result_dir/index/$dataset/$index_type"
mkdir -p "$result_dir/log/$dataset/$index_type"
#
## build index
rm -r $code_dir
cmake -B $code_dir -S $source_dir
sleep 1
cd $code_dir
make -j
#16 32 62 63 64 126 127 128 254 255 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144
#64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072
#16 32 62 63 64 126 127 128 254 255 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144
#"simd" "normal" "simd_simple" "ccnews-1M"

for dataset in "cw12b-1M" "ccnews-1M"
do
  for query_type in "AND" "OR"
  do
    for query_num in 5 4 3 2
    do
      for epsilon in 1
      do
        for decode_type in "normal"
        do
          echo "————————————dataset : $dataset epsilon: $epsilon query_type: $query_type query_num: $query_num decode_type: $decode_type————————————"
          for repeat in 1 2
          do
            mkdir -p $result_dir/index/$dataset/$index_type/$dataset-$index_type-$epsilon/
            $code_dir/pgm_query $index_type $data_dir/$dataset/$dataset $result_dir/index/$dataset/$index_type/$dataset-$index_type-$epsilon/ $epsilon $read_only $decode_type $data_dir/$dataset/$dataset-$query_num.queries $query_type $result_dir/log/$dataset/$index_type/$dataset-$index_type-$epsilon.query-$query_num-log.txt
          done
          echo " "
        done
      done
    done
  done
done
