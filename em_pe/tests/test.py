from em_pe import sampler
from em_pe.parser import parse_json
from em_pe.plot_utils import plot_corner

### Parse JSON data and convert to necessary format
parse_json.parse_json(57982.5285231481, "Data/GW170817.json", ["H", "K"],
                      "em_pe/tests/temp/", maxpts=50, tmax=8)

### Initialize sampler
s = sampler.sampler("em_pe/tests/temp/", [["woko2017", 1]], ["H.txt", "K.txt"],
                    "results/woko2017/posterior_samples.txt", min_iter=5, max_iter=5)

### Generate samples
s.generate_samples()

### Make a corner plot
plot_corner.generate_plot(["results/woko2017/posterior_samples.txt"],
                           "results/woko2017/corner.png", ["mej", "vej"],
                           leg=["test"])
