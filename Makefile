directories:
	mkdir -p pe_runs/

TMIN = 0.1
TMAX = 8.0
MEJ = 0.01
VEJ = 0.1
LOG_MEJ_BLUE = -1.7
LOG_MEJ_PURPLE = -1.3
LOG_MEJ_RED = -2.0
VEJ_BLUE = 0.25
VEJ_PURPLE = 0.15
VEJ_RED = 0.15
TC_BLUE = 700.0
TC_PURPLE = 1300.0
TC_RED = 3700.0
DIST = 40.0
M1 = 1.4
M2 = 1.35

### general injection parameters
INJECTION_PARAMS = --n 50 --err 0.2 --time-format mjd --tmin ${TMIN} --tmax ${TMAX} --p dist ${DIST}
EJECTA_PARAMS = --p mej ${MEJ} --p vej ${VEJ}
EJECTA_PARAMS_3C = --p log_mej_red ${LOG_MEJ_RED} --p log_mej_purple ${LOG_MEJ_PURPLE} --p log_mej_blue ${LOG_MEJ_BLUE} --p vej_red ${VEJ_RED} --p vej_purple ${VEJ_PURPLE} --p vej_blue ${VEJ_BLUE} --p Tc_red ${TC_RED} --p Tc_purple ${TC_PURPLE} --p Tc_blue ${TC_BLUE}
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

test_kilonova_3c: directories
	mkdir -p pe_runs/$@/
	python3 em_pe/tests/generate_data.py --m kilonova_3c --out pe_runs/$@/ ${INJECTION_PARAMS} ${EJECTA_PARAMS_3C}
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m kilonova_3c -v --f g.txt --f r.txt --f i.txt --f z.txt --f y.txt --f J.txt --f H.txt --f K.txt --min 40 --max 40 --out samples.txt --estimate-dist --fixed-param Tc_red 3700.0 --fixed-param Tc_purple 1300.0 --fixed-param Tc_blue 700.0 --correlate-dims log_mej_red vej_red --correlate-dims log_mej_purple vej_purple --correlate-dims log_mej_blue vej_blue --nprocs 8 --burn-in 10 10 --beta-start 0.1 --keep-npts 100000" > pe_runs/$@/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples.txt --truth-file test_truths.txt --out corner.png --p log_mej_red --p log_mej_purple --p log_mej_blue --p vej_red --p vej_purple --p vej_blue --p dist" > pe_runs/$@/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior-samples samples.txt --out lc.png --m kilonova_3c --tmin ${TMIN} --tmax ${TMAX} --lc-file g.txt --b g --lc-file r.txt --b r --lc-file i.txt --b i --lc-file z.txt --b z --lc-file y.txt --b y --lc-file J.txt --b J --lc-file H.txt --b H --lc-file K.txt --b K --fixed-param Tc_red 3700.0 --fixed-param Tc_purple 1300.0 --fixed-param Tc_blue 700.0" > pe_runs/$@/plot_lc.sh
	chmod u+x pe_runs/$@/sample.sh
	chmod u+x pe_runs/$@/plot_corner.sh
	chmod u+x pe_runs/$@/plot_lc.sh

### GW170817

GW170817_START = 1187008882.43
TMAX = 20.0

BAND = K

GW170817_kilonova: directories
	mkdir -p pe_runs/$@/
	python3 ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --t0 ${GW170817_START} --f ${EM_PE_INSTALL_DIR}/Data/GW170817.json --b g --b r --b i --b z --b y --b J --b H --b K --out pe_runs/$@/
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m kilonova -v --f K.txt --min 20 --max 20 --out samples_${BAND}.txt --estimate-dist" > pe_runs/$@/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples_${BAND}.txt --out corner_${BAND}.png --p log_mej --p vej --p log_kappa --p Tc --p dist" > pe_runs/$@/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior-samples samples_${BAND}.txt --out lc_${BAND}.png --m kilonova --tmin ${TMIN} --tmax ${TMAX} --lc-file ${BAND}.txt --b ${BAND}" > pe_runs/$@/plot_lc.sh
	chmod u+x pe_runs/$@/sample.sh
	chmod u+x pe_runs/$@/plot_corner.sh
	chmod u+x pe_runs/$@/plot_lc.sh

GW170817_kilonova_2c: directories
	mkdir -p pe_runs/$@/
	python3 ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --t0 ${GW170817_START} --f ${EM_PE_INSTALL_DIR}/Data/GW170817.json --b g --b r --b i --b z --b y --b J --b H --b K --maxpts 20 --out pe_runs/$@/
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m kilonova_2c -v --f K.txt --min 20 --f r.txt --max 20 --out samples.txt --estimate-dist --fixed-param Tc_red 1000.0 --fixed-param Tc_blue 4000.0" > pe_runs/$@/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples.txt --out corner.png --p log_mej_red --p log_mej_blue --p vej_red --p vej_blue" > pe_runs/$@/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior-samples samples.txt --out lc.png --m kilonova_2c --tmin ${TMIN} --tmax ${TMAX} --lc-file K.txt --b K --lc-file r.txt --b r" > pe_runs/$@/plot_lc.sh
	chmod u+x pe_runs/$@/sample.sh
	chmod u+x pe_runs/$@/plot_corner.sh
	chmod u+x pe_runs/$@/plot_lc.sh

GW170817_kilonova_3c: directories
	mkdir -p pe_runs/$@/
	python3 ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --t0 ${GW170817_START} --f ${EM_PE_INSTALL_DIR}/Data/GW170817.json --b g --b r --b i --b z --b y --b J --b H --b K --out pe_runs/$@/
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m kilonova_3c -v --f g.txt --f r.txt --f i.txt --f z.txt --f y.txt --f J.txt --f H.txt --f K.txt --min 60 --max 60 --out samples.txt --estimate-dist --burn-in 20 20 --keep-npts 1000000 --fixed-param Tc_red 3700.0 --fixed-param Tc_purple 1300.0 --fixed-param Tc_blue 700.0 --correlate-dims log_mej_red vej_red --correlate-dims log_mej_purple vej_purple --correlate-dims log_mej_blue vej_blue" > pe_runs/$@/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples.txt --out corner.png --p log_mej_red --p log_mej_purple --p log_mej_blue --p vej_red --p vej_purple --p vej_blue --p dist" > pe_runs/$@/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --posterior-samples samples.txt --out lc.png --m kilonova_3c --tmin ${TMIN} --tmax ${TMAX} --lc-file g.txt --b g --lc-file r.txt --b r --lc-file i.txt --b i --lc-file z.txt --b z --lc-file y.txt --b y --lc-file J.txt --b J --lc-file H.txt --b H --lc-file K.txt --b K --fixed-param Tc_red 3700.0 --fixed-param Tc_purple 1300.0 --fixed-param Tc_blue 700.0" > pe_runs/$@/plot_lc.sh
	chmod u+x pe_runs/$@/sample.sh
	chmod u+x pe_runs/$@/plot_corner.sh
	chmod u+x pe_runs/$@/plot_lc.sh
