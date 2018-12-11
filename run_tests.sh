python em_pe/tests/imports_test.py
python em_pe/tests/models_test.py

### Set up test data and results directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/
python em_pe/tests/generate_test_data.py

### Do a normal test in one band with a single model
python em_pe/generate_posterior_samples.py --dat em_pe/tests/temp/ -v --m one_band_test 1 --f test_bandA.txt --out em_pe/tests/temp/posterior_samples.txt --cutoff 1e-8

### Do a test in two bands with a single model
#python em_pe/generate_posterior_samples.py --min 100 --max 100 --dat em_pe/tests/temp/ -v --m two_band_test 1 --f test_bandA.txt --f test_bandB.txt --out em_pe/tests/temp/posterior_samples_2.txt --cutoff 1e-12

### Make the corner plots
python em_pe/plot_utils/plot_corner.py --p a --p b --posterior_samples em_pe/tests/temp/posterior_samples.txt --out em_pe/tests/temp/fig.png --truth_file em_pe/tests/temp/test_truths.txt
