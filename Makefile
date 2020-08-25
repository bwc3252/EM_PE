directories:
	mkdir -p pe_runs/

TMIN = 0.1
TMAX = 8.0
MEJ = 0.01
VEJ = 0.1
MEJ_BLUE = 0.020
MEJ_PURPLE = 0.050
MEJ_RED = 0.010
VEJ_BLUE = 0.25
VEJ_PURPLE = 0.15
VEJ_RED = 0.15
TC_BLUE = 700.0
TC_PURPLE = 1300.0
TC_RED = 3700.0
SIGMA = 0.1
DIST = 40.0
M1 = 1.4
M2 = 1.35

### general injection parameters
INJECTION_PARAMS = --n 50 --err 0.2 --sigma ${SIGMA} --time-format mjd --tmin ${TMIN} --tmax ${TMAX} --p dist ${DIST}
EJECTA_PARAMS = --p mej ${MEJ} --p vej ${VEJ}
EJECTA_PARAMS_3C = --p mej_red ${MEJ_RED} --p mej_purple ${MEJ_PURPLE} --p mej_blue ${MEJ_BLUE} --p vej_red ${VEJ_RED} --p vej_purple ${VEJ_PURPLE} --p vej_blue ${VEJ_BLUE} --p Tc_red ${TC_RED} --p Tc_purple ${TC_PURPLE} --p Tc_blue ${TC_BLUE} --p sigma ${SIGMA}
BNS_PARAMS = --p m1 ${M1} --p m2 ${M2}

test_kilonova_3c: directories
	mkdir -p pe_runs/$@_$(shell date +%Y%m%d)/
	python3 em_pe/tests/generate_data.py --m kilonova_3c --out pe_runs/$@_$(shell date +%Y%m%d)/ ${INJECTION_PARAMS} ${EJECTA_PARAMS_3C}
	echo "time python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m kilonova_3c -v --f g.txt --f r.txt --f i.txt --f z.txt --f y.txt --f J.txt --f H.txt --f K.txt --min 80 --max 80 --out samples.txt --fixed-param dist 40.0 --correlate-dims mej_red vej_red --correlate-dims mej_purple vej_purple --correlate-dims mej_blue vej_blue --burn-in 10 --beta-start 0.01 --keep-npts 2000000 --nprocs 8 --fixed-param Tc_red ${TC_RED} --fixed-param Tc_purple ${TC_PURPLE} --fixed-param Tc_blue ${TC_BLUE}" > pe_runs/$@_$(shell date +%Y%m%d)/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples.txt --truth-file test_truths.txt --out corner.png --p mej_red --p mej_purple --p mej_blue --p vej_red --p vej_purple --p vej_blue --p sigma" > pe_runs/$@_$(shell date +%Y%m%d)/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --log-time --posterior-samples samples.txt --out lc.png --m kilonova_3c --tmin ${TMIN} --tmax ${TMAX} --lc-file g.txt --b g --lc-file r.txt --b r --lc-file i.txt --b i --lc-file z.txt --b z --lc-file y.txt --b y --lc-file J.txt --b J --lc-file H.txt --b H --lc-file K.txt --b K --fixed-param dist 40.0 --fixed-param Tc_red ${TC_RED} --fixed-param Tc_purple ${TC_PURPLE} --fixed-param Tc_blue ${TC_BLUE}" > pe_runs/$@_$(shell date +%Y%m%d)/plot_lc.sh
	chmod u+x pe_runs/$@_$(shell date +%Y%m%d)/sample.sh
	chmod u+x pe_runs/$@_$(shell date +%Y%m%d)/plot_corner.sh
	chmod u+x pe_runs/$@_$(shell date +%Y%m%d)/plot_lc.sh

### GW170817

GW170817_START = 1187008882.43
TMAX = 30.0

GW170817_kilonova_3c: directories
	mkdir -p pe_runs/$@_$(shell date +%Y%m%d)/
	python3 ${EM_PE_INSTALL_DIR}/em_pe/parser/parse_json.py --t0 ${GW170817_START} --f ${EM_PE_INSTALL_DIR}/Data/GW170817.json --b g --b r --b i --b z --b y --b J --b H --b K --out pe_runs/$@_$(shell date +%Y%m%d)/
	echo "time python3 ${EM_PE_INSTALL_DIR}/em_pe/sampler.py --dat ./ --m kilonova_3c -v --f g.txt --f r.txt --f i.txt --f z.txt --f y.txt --f J.txt --f H.txt --f K.txt --min 60 --max 60 --out samples.txt --fixed-param dist 40.0 --correlate-dims mej_red vej_red --correlate-dims mej_purple vej_purple --correlate-dims mej_blue vej_blue --burn-in 20 --beta-start 0.01 --nprocs 8 --keep-npts 500000 --set-limit mej_blue 0.01 0.02 --set-limit mej_purple 0.037 0.043 --set-limit mej_red 0.009 0.015 --set-limit vej_red 0.1 0.18 --set-limit vej_purple 0.122 0.138 --set-limit vej_blue 0.25 0.28 --set-limit Tc_red 2500.0 3700.0 --set-limit Tc_purple 1000.0 1250.0 --set-limit Tc_blue 400.0 1200.0 --ncomp Tc_blue 2 " > pe_runs/$@_$(shell date +%Y%m%d)/sample.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_corner.py --posterior-samples samples.txt --out corner.png --p mej_red --p mej_purple --p mej_blue --p vej_red --p vej_purple --p vej_blue --p sigma --p Tc_red --p Tc_purple --p Tc_blue" > pe_runs/$@_$(shell date +%Y%m%d)/plot_corner.sh
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/plot_lc.py --log-time --posterior-samples samples.txt --out lc.png --m kilonova_3c --tmin ${TMIN} --tmax ${TMAX} --lc-file g.txt --b g --lc-file r.txt --b r --lc-file i.txt --b i --lc-file z.txt --b z --lc-file y.txt --b y --lc-file J.txt --b J --lc-file H.txt --b H --lc-file K.txt --b K --fixed-param dist 40.0" > pe_runs/$@_$(shell date +%Y%m%d)/plot_lc.sh
	chmod u+x pe_runs/$@_$(shell date +%Y%m%d)/sample.sh
	chmod u+x pe_runs/$@_$(shell date +%Y%m%d)/plot_corner.sh
	chmod u+x pe_runs/$@_$(shell date +%Y%m%d)/plot_lc.sh

kilonova_3c_pp_plot: directories
	chmod u+x pe_runs/pp_plot/run.sh
	python3 ${EM_PE_INSTALL_DIR}/scripts/pp_plot_helper.py --m kilonova_3c --directory pe_runs/ --name pp_plot_$(shell date +%Y%m%d) --npts 100 --sigma 0.1 --sampler-args "--correlate-dims mej_red vej_red --correlate-dims mej_purple vej_purple --correlate-dims mej_blue vej_blue" --fixed-param dist ${DIST} --fixed-param Tc_red ${TC_RED} --fixed-param Tc_purple ${TC_PURPLE} --fixed-param Tc_blue ${TC_BLUE}
	echo "python3 ${EM_PE_INSTALL_DIR}/em_pe/plot_utils/pp_plot.py --m kilonova_3c --name \"\" --directory ./" > pe_runs/pp_plot_$(shell date +%Y%m%d)/plot.sh
	chmod u+x pe_runs/pp_plot_$(shell date +%Y%m%d)/run.sh
	chmod u+x pe_runs/pp_plot_$(shell date +%Y%m%d)/plot.sh
