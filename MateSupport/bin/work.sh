echo "test"
python MateSupport.py --bam ../input/test.bam --type inward --insert 2385 --sd 490

echo "HEG4_RAW.3.3k"
python MateSupport.py --bam ../input/out.smaltmap.bam --type inward --insert 2385 --sd 490 --pair 12188424 --project HEG4_RAW.3.3k


