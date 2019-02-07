### basic script to run PE
# Usage: sh run.sh [test name]

### Set up test data and results directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/

mkdir results/$1/

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
python em_pe/generate_posterior_samples.py --min 5 --max 5 --dat em_pe/tests/temp/ \
                                           -v --m woko2017 1 --out results/$1/posterior_samples.txt \
                                           --cutoff 0 \
                                           --f H.txt \
                                           --f K.txt

                                           #--f g.txt \
                                           #--f r.txt \
                                           #--f i.txt \
                                           #--f z.txt \
                                           #--f y.txt \
                                           #--f J.txt \
                                           #--f H.txt \
                                           #--f K.txt

### Make the corner plots
python em_pe/plot_utils/plot_corner.py --p mej --p vej --posterior_samples results/$1/posterior_samples.txt \
                                       --out results/$1/corner.png
