### Set up test data and results directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/
python em_pe/tests/generate_woko2017_data.py

### Do a normal test in one band with a single model
python em_pe/generate_posterior_samples.py --dat em_pe/tests/temp/ -v --m woko2017 1 \
                                           --f 4775.6.txt \
                                           --f 6129.5.txt \
                                           --f 7484.6.txt \
                                           --f 8657.8.txt \
                                           --f 9603.1.txt \
                                           --f 12350.txt \
                                           --f 16620.txt \
                                           --f 21590.txt \
                                           --out em_pe/tests/temp/posterior_samples.txt --cutoff 1e-8

### Make the corner plots
python em_pe/plot_utils/plot_corner.py --p mej --p vej --posterior_samples em_pe/tests/temp/posterior_samples.txt \
                                       --out em_pe/tests/temp/fig.png --truth_file em_pe/tests/temp/test_truths.txt
