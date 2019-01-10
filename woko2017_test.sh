### Set up test data and results directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/
#python em_pe/tests/generate_woko2017_data.py

### parse JSON data
python em_pe/parser/parse_json.py --f Data/GW170817.json --out em_pe/tests/temp/ --t0 57982.5285231481 \
                                --b y \
                                --b g \
                                --b r \
                                --b z \
                                --b K \
                                --b J \
                                --b H \
                                --b i


### Do a normal test in one band with a single model
python em_pe/generate_posterior_samples.py --min 5 --max 5 --dat em_pe/tests/temp/ -v --m woko2017 1 \
                                           --out em_pe/tests/temp/posterior_samples.txt --cutoff 1e-400 \
                                           --f y.txt \
                                           --f g.txt \
                                           --f r.txt \
                                           --f z.txt \
                                           --f K.txt \
                                           --f J.txt \
                                           --f H.txt \
                                           --f i.txt

### Make the corner plots
python em_pe/plot_utils/plot_corner.py --p mej --p vej --p dist --posterior_samples em_pe/tests/temp/posterior_samples.txt \
                                       --out em_pe/tests/temp/fig.png #--truth_file em_pe/tests/temp/test_truths.txt
