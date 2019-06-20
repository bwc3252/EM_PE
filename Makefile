directories:
	mkdir -p pe_runs/

TMIN = 0.1
TMAX = 8.0

### general injection parameters

INJECTION_PARAMS = --n 50 --err 0.2 --time_format mjd --tmin ${TMIN} --tmax ${TMAX}

### woko2017 injection test using z band

WOKO2017_INJECTION_PARAMS = --p mej 0.05 --p vej 0.2 --p dist 40.0

test_woko2017: directories
	mkdir -p pe_runs/$@/
	python em_pe/tests/generate_data.py --m woko2017 --out pe_runs/$@ ${INJECTION_PARAMS} ${WOKO2017_INJECTION_PARAMS}
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m woko2017 -v --f z.txt --min 20 --max 20 --out samples_z.txt" > pe_runs/$@/sample.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior_samples samples_z.txt --truth_file test_truths.txt --out corner_z.png --p mej --p vej --p dist" > pe_runs/$@/plot_corner.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior_samples samples_z.txt --out lc_z.png --m woko2017 --tmin ${TMIN} --tmax ${TMAX} --lc_file z.txt --b z" > pe_runs/$@/plot_lc.sh
