# SPNET

Survival P-NET, which implements P-NET (a biologically informed deep neural network)

## Prerequisites

Ensure you have the following dependencies installed:

- Python 3.10
- `networkx`
- `sklearn`
- `torch`
- `sksurv`
- `captum`

Install missing packages using:
```bash
pip install networkx scikit-learn torch scikit-survival captum
```

## Database and Dataset

### Reactome2021
The `Reactome2021` dataset is sourced directly from [Elmarakeby's work](https://github.com/marakeby/pnet_prostate_paper). It provides hierarchical relationships between genes and different pathway levels.

Assign the correct path to `reactome_base_dir` to access these hierarchical relationships.

### TCGA-LIHC
The [TCGA-LIHC dataset](https://drive.google.com/drive/folders/1yCvFoF-xzNpdkzhXaXRrxapqvhaSHSVP?usp=sharing) is used in this study. Data processing steps are recorded in `makedata.R`.


## License
This project is licensed under the MIT License.
