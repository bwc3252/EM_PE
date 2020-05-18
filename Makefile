directories:
	mkdir -p pe_runs/

TMIN = 0.1
TMAX = 8.0
MEJ = 0.01
VEJ = 0.1
DIST = 40.0
M1 = 1.4
M2 = 1.35

### general injection parameters
INJECTION_PARAMS = --n 50 --err 0.2 --time-format mjd --tmin ${TMIN} --tmax ${TMAX} --p dist ${DIST}
EJECTA_PARAMS = --p mej ${MEJ} --p vej ${VEJ}
BNS_PARAMS = --p m1 ${M1} --p m2 ${M2}

test_woko2017: directories
	mkdir -p pe_runs/$@/
	python3 em_pe/tests/generate_data.py --m woko2017 --out pe_runs/$@/ ${INJECTION_PARAMS} ${EJECTA_PARAMS}
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m woko2017 -v --f z.txt --min 20 --max 20 --out samples_z.txt --fixed-param dist 40.0" > pe_runs/$@/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples_z.txt --truth-file test_truths.txt --out corner_z.png --p mej --p vej" > pe_runs/$@/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior-samples samples_z.txt --out lc_z.png --m woko2017 --tmin ${TMIN} --tmax ${TMAX} --lc-file z.txt --b z --fixed-param dist 40.0" > pe_runs/$@/plot_lc.sh
	chmod u+x pe_runs/$@/sample.sh
	chmod u+x pe_runs/$@/plot_corner.sh
	chmod u+x pe_runs/$@/plot_lc.sh

test_me2017_non_lanthanide: directories
	mkdir -p pe_runs/$@/
	python3 em_pe/tests/generate_data.py --m me2017_non_lanthanide --out pe_runs/$@/ ${INJECTION_PARAMS} ${EJECTA_PARAMS}
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m me2017_non_lanthanide -v --f z.txt --min 20 --max 20 --out samples_z.txt --fixed-param dist 40.0" > pe_runs/$@/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples_z.txt --truth-file test_truths.txt --out corner_z.png --p mej --p vej" > pe_runs/$@/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior-samples samples_z.txt --out lc_z.png --m me2017_non_lanthanide --tmin ${TMIN} --tmax ${TMAX} --lc-file z.txt --b z --fixed-param dist 40.0" > pe_runs/$@/plot_lc.sh
	chmod u+x pe_runs/$@/sample.sh
	chmod u+x pe_runs/$@/plot_corner.sh
	chmod u+x pe_runs/$@/plot_lc.sh
