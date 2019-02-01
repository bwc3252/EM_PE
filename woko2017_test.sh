### Set up test data and results directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/

### Make fake data
python em_pe/tests/generate_woko2017_data.py

### Do a normal test in one band with a single model
python em_pe/generate_posterior_samples.py --min 5 --max 5 --dat em_pe/tests/temp/ -v --m woko2017 1 \
                                           --out em_pe/tests/temp/posterior_samples.txt --cutoff 1e-1000 \
                                           --f g.txt \
                                           --f H.txt \
                                           --f i.txt

### Make the corner plots
python em_pe/plot_utils/plot_corner.py --p mej --p vej --p delta_t --posterior_samples em_pe/tests/temp/posterior_samples.txt \
                                       --out em_pe/tests/temp/fig.png --truth_file em_pe/tests/temp/test_truths.txt
