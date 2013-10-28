echo "sample"
python ChunkGenome.py --reference ../input/MSU_r7.fa --assembly assembly.fa --project HEG4.ALLPATHLG > log 2> log2

echo "HEG4_RAW.3"
python ChunkGenome.py --reference ../input/MSU_r7.fa --assembly ../input/HEG4_RAW.3.fa --project HEG4.ALLPATHLG > log 2> log2 &
