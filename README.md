# Inference of differential splicing factor (SF) activity with SF-exon networks and VIPER

## R
```{R}
require(tidyverse)
require(viper)
source("scripts/util.R")

# load exon inclusion table and metadata

# compute signature
signature = read_tsv(signature_file)

# load SF-exon networks
sf_networks = load_networks(regulons_path, n_tails)

# compute activities
protein_activities = viper(signature, regulons, verbose=FALSE)

print(protein_activities)
```

## Python
```{python}

```