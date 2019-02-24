from em_pe import sampler
from em_pe.parser import parse_json
from em_pe.plot_utils import plot_corner

### Parse JSON data and convert to necessary format
#parse_json.parse_json(57982.5285231481, "Data/GW170817.json", ["H"],
#                      "em_pe/tests/temp/", maxpts=50, tmax=8)

### Initialize sampler
s = sampler.sampler("em_pe/tests/temp/", "two_comp", ["H.txt"],
                    "results/two_comp/posterior_samples.txt", min_iter=100, max_iter=100)

### Generate samples
s.generate_samples()

### Make a corner plot
plot_corner.generate_plot(["results/two_comp/posterior_samples.txt"],
                           "results/two_comp/corner.png", ["mej1", "vej1", "mej2", "vej2", "frac"],
                           leg=["test"])
