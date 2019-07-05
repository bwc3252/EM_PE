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
	python em_pe/tests/generate_data.py --m woko2017 --out pe_runs/$@/ ${INJECTION_PARAMS} ${WOKO2017_INJECTION_PARAMS}
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --epoch 5 --dat ./ --m woko2017 -v --f z.txt --min 20 --max 20 --out samples_z.txt --estimate_dist" > pe_runs/$@/sample.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior_samples samples_z.txt --truth_file test_truths.txt --out corner_z.png --p mej --p vej --p dist" > pe_runs/$@/plot_corner.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior_samples samples_z.txt --out lc_z.png --m woko2017 --tmin ${TMIN} --tmax ${TMAX} --lc_file z.txt --b z" > pe_runs/$@/plot_lc.sh

### GW170817 runs

GW170817_START = 1187008882.43

TMIN = 0.5

GW170817_woko2017_fixed_dist: directories
	mkdir -p pe_runs/$@/
	python ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --t0 ${GW170817_START} --f ${EM_PE_INSTALL_DIR}/Data/GW170817.json --b r --out pe_runs/$@/ --maxpts 100 --tmax ${TMAX}
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --epoch 5 --dat ./ --m woko2017 -v --f r.txt --min 20 --max 20 --out samples_r.txt --fixed_param dist 40" > pe_runs/$@/sample.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior_samples samples_r.txt --out corner_r.png --p mej --p vej" > pe_runs/$@/plot_corner.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior_samples samples_r.txt --out lc_r.png --m woko2017 --tmin ${TMIN} --tmax ${TMAX} --lc_file r.txt --b r --fixed_param dist 40" > pe_runs/$@/plot_lc.sh

GW170817_woko2017_fixed_dist_orientation: directories
	mkdir -p pe_runs/$@/
	python ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --t0 ${GW170817_START} --f ${EM_PE_INSTALL_DIR}/Data/GW170817.json --b r --out pe_runs/$@/ --maxpts 100 --tmax ${TMAX}
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --orientation gaussian --epoch 5 --dat ./ --m woko2017 -v --f r.txt --min 20 --max 20 --out samples_r.txt --fixed_param dist 40" > pe_runs/$@/sample.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior_samples samples_r.txt --out corner_r.png --p mej --p vej --p angle" > pe_runs/$@/plot_corner.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior_samples samples_r.txt --out lc_r.png --m woko2017 --tmin ${TMIN} --tmax ${TMAX} --lc_file r.txt --b r --fixed_param dist 40" > pe_runs/$@/plot_lc.sh

### GRB160821B

GRB160821B_START = 0.0

GRB160821B_woko2017: directories
	mkdir -p pe_runs/$@/
	python ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --time_format mjd --t0 ${GRB160821B_START} --f ${EM_PE_INSTALL_DIR}/Data/GRB160821B.json --b g --out pe_runs/$@/ --maxpts 100 --tmax ${TMAX}
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --epoch 5 --dat ./ --m woko2017 -v --f g.txt --min 20 --max 20 --out samples_g.txt --fixed_param dist 776" > pe_runs/$@/sample.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior_samples samples_g.txt --out corner_g.png --p mej --p vej" > pe_runs/$@/plot_corner.sh
	echo "python ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior_samples samples_g.txt --out lc_g.png --m woko2017 --tmin ${TMIN} --tmax ${TMAX} --lc_file g.txt --b g --fixed_param dist 776" > pe_runs/$@/plot_lc.sh

