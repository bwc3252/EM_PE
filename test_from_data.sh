### Set up temp directory
rm -rf em_pe/tests/temp/
mkdir em_pe/tests/temp/

### Parse JSON file
python em_pe/parser/parse_json.py --f Data/GRB060614.json --b R --out em_pe/tests/temp/ --t0 53902

### Generate posterior samples
python em_pe/generate_posterior_samples.py --dat em_pe/tests/temp/ -v --m linear 1 --f R.txt --out em_pe/tests/temp/posterior_samples.txt --cutoff 1e-118

### Create corner plot
python em_pe/plot_utils/plot_corner.py --p m --p y0 --posterior_samples em_pe/tests/temp/posterior_samples.txt --out em_pe/tests/temp/fig.png
