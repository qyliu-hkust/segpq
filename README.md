# SegPQ

## Install `faiss` library (conda is recommended)
```
conda create -n segpq python=3.11
conda install -c pytorch faiss-cpu (for CPU version)
conda install -c pytorch faiss-gpu (for GPU version)
```

## Prepare raw dataset

Four datasets are adopted:

- **Audio**: https://research.google.com/youtube8m/download.html
- **SIGT/GIST**: http://corpus-texmex.irisa.fr/
- **Deep1B**: https://research.yandex.com/blog/benchmarks-for-billion-scale-similarity-search

Statistics of the datasets are summarized as follows:
| Dataset      | #Vectors | #Dimensions     |  Raw Size     |
| :---        |    :----:   |          ---: |          ---: |
| GIST      | 1,000,000       | 960   | 3.7 GiB |
| Audio   | 438,229,156        | 128      | 211 GiB |
| SIFT | 1,000,000,000 | 128 | 123 GiB |
| Deep1B | 1,000,000,000 | 128 | 361 GiB |


## Compile the code to generate codebook
```
g++ faiss.cpp --std=c++20 -I./include -I./lib/sdsl-lite/include -I/opt/anaconda3/include -lfaiss -Wl,-rpath,/opt/anaconda3/lib -o faiss_gen
```
Note that: you need to pass your own path to `faiss.h` and `libfaiss.so` (on MacOS should be `libfaiss.dylib`). 

## Prepare codebooks
```
./faiss_gen RAW_VEC_FILE NUM_POINTS NUM_DIM PQ_TYPE PQ_DIM OUTPUT_FILE
```

For example: 
```
./faiss_gen ./sift_base.fvecs 1000000 128 OPQ 32 ./sift_1m_pq32
```

## Run SegPQ Framework

Compile SegPQ:
```
bash ./compressor/build.sh
```

Compress PQ codebook:
```
bash ./compressor/compress.sh
```

Decode the compressed PQ codebook:
```
bash ./compressor/decode.sh
```

Note that: change the `data_dir` in above shell files as the path containing the generated PQ codebooks.
