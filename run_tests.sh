python em_pe/tests/imports_test.py
python em_pe/tests/models_test.py

rm -rf em_pe/tests/Data/
mkdir em_pe/tests/Data/
python em_pe/tests/generate_test_data.py

python em_pe/generate_posterior_samples.py --dat em_pe/tests/Data/ -v --m one_band_test 1 --f test_bandA.txt --out em_pe/tests/Data/posterior_samples.txt --cutoff 1e-8
#python em_pe/generate_posterior_samples.py --min 100 --max 100 --dat em_pe/tests/Data/ -v --m two_band_test 1 --f test_bandA.txt --f test_bandB.txt --out em_pe/tests/Data/posterior_samples.txt --cutoff 1e-12
#python em_pe/generate_posterior_samples.py --dat em_pe/tests/Data/ -v --m one_band_test 0.5 --m one_band_test 0.5 --f test_bandA.txt --out em_pe/tests/Data/posterior_samples.txt --cutoff 1e-8

python em_pe/plot_utils/plot_corner.py --posterior_samples em_pe/tests/Data/posterior_samples.txt --out em_pe/tests/Data/fig.png --truth_file em_pe/tests/Data/test_truths.txt
