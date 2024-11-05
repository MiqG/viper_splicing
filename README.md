# Inference of differential splicing factor (SF) activity with SF-exon networks and VIPER

This repository contains the SF-exon networks obtained from experiments to estimate differential splicing factor activity from exon inclusion matrices (with [VastDB](https://vastdb.crg.eu/wiki/Downloads) exon identifiers, for the moment).


## Installation
Just download SF-exon networks by cloning the repository and run regular VIPER in either R or Python.

```{shell}
git clone https://github.com/MiqG/viper_splicing.git
```

## R
### Requirements
- [`viper`](https://www.bioconductor.org/packages/release/bioc/html/viper.html)

```{shell}
mamba install bioconda::bioconductor-viper
```

### Usage
```{R}
require(tidyverse)
require(viper)
source("scripts/util.R")

# load example exon inclusion table and metadata
psi = read_tsv("data/event_psi/Nijhuis2020-EX.tsv.gz")
metadata = read_tsv("data/metadata/Nijhuis2020.tsv.gz")

# load SF-exon networks
sf_networks = load_networks("data/empirical_sf_networks-EX")

# compute signature by subtracting average PSI in control samples 
psi = psi %>% column_to_rownames("EVENT")
ctl_samples = metadata %>% filter(condition == "DMSO") %>% pull(sampleID)
psi_ctl = psi[,ctl_samples] %>% rowMeans()
signature = sweep(psi, MARGIN=1, STATS=psi_ctl, FUN="-")

# compute activities
protein_activities = viper(signature, sf_networks, verbose=FALSE)

print("Done!")
```

## Python (`pyviper` is still under development)
### Requirements
- [`pyviper`](https://github.com/alevax/pyviper)

```{shell}
pip install viper-in-python
```

### Usage
```{python}
import os
import pandas as pd
import scanpy
import pyviper

# load eexample xon inclusion table and metadata
psi = pd.read_table("data/event_psi/Nijhuis2020-EX.tsv.gz")
metadata = pd.read_table("data/metadata/Nijhuis2020.tsv.gz")

# load SF-exon networks
networks_dir = "data/empirical_sf_networks-EX"
sf_networks = pd.concat([
    pd.read_table(os.path.join(networks_dir,f))
    for f in os.listdir(networks_dir)
])
sf_networks = sf_networks.rename(columns={"tfmode":"mor"})
sf_networks = sf_networks[["regulator","target","mor","likelihood"]]
interactome = pyviper.Interactome("net", sf_networks)

# compute signature by subtracting average PSI in control samples
psi = psi.set_index("EVENT")
ctl_samples = metadata.loc[metadata["condition"]=="DMSO","sampleID"].to_list()
psi_ctl = psi[ctl_samples].mean(axis=1).values.reshape(-1,1)
signature = psi - psi_ctl
signature = scanpy.AnnData(signature.T.fillna(0)) # It does not work without imputation

# compute activities
interactome.filter_targets(signature.var_names)
protein_activities = pyviper.viper(gex_data=signature, interactome=interactome, enrichment="area", verbose=False)
protein_activities = protein_activities.to_df().T.reset_index().rename(columns={"index":"regulator"})

print("Done!")
```

## Citation
Anglada-Girotto, M., Moakley, D. F., Zhang, C., Miravet-Verde, S., Califano, A., & Serrano, L. (2024). Disentangling the splicing factor programs underlying complex molecular phenotypes. bioRxiv, 2024-06.