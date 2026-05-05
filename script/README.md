# Collection of helper script for downstream analysis
This folder contain several usefull script for downstream analysis support. I will keep adding usefull script such as visualization to support many downstream analysis process.

## 1. `trim_tir.py`
By default, retained GEVEs sequence by findGEVE are including TIR sequence in both sides. In case you don't want to remove TIR sequence before downstream analysis, you can use this script.
```
python trim_tir.py input.geve.fna
```

## 2. `extract_tir.py`
Sometimes, you may want to see wheter retained GEVE sequence has commot TIR sequence or not. Here is a script to extract TIR sequence for each GEVE and concatenate it.
```
extract_tir.py extract_tir.py
```

Ouput structure for each extracted TIR sequence will look like this
```
>GEVE1
[TIR-L][revcomp(TIR-R)]
>GEVE2
...
```