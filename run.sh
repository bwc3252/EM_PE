### Set up test data and results directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/

### parse JSON data
python em_pe/parser/parse_json.py --f Data/GW170817.json --out em_pe/tests/temp/ \
                                  --t0 57982.5285231481 --maxpts 50  --tmax 8 \
                                --b g \
                                --b r \
                                --b i \
                                --b z \
                                --b y \
                                --b J \
                                --b H \
                                --b K

### Generate the posterior samples
python em_pe/generate_posterior_samples.py --min 200 --max 200 --dat em_pe/tests/temp/ \
                                           -v --m woko2017 1 --out em_pe/tests/temp/posterior_samples.txt \
                                           --cutoff 0 \
                                           --f J.txt

                                           #--f g.txt \
                                           #--f r.txt \
                                           #--f i.txt \
                                           #--f z.txt \
                                           #--f y.txt \
                                           #--f J.txt \
                                           #--f H.txt \
                                           #--f K.txt

### Make the corner plots
python em_pe/plot_utils/plot_corner.py --p mej --p vej --p delta_t --posterior_samples em_pe/tests/temp/posterior_samples.txt \
                                       --out em_pe/tests/temp/fig.png
