# SingleCell

SingleCell is a python class available in the singlecelldata package for managing single-cell RNA-seq data. It contains three pandas dataframes; `data` for holding gene expression values (counts/normalized counts), `genedata` for holding more information about the genes e.g., gene names, and `celldata` which contains more information about cells such as cell types, labels etc.

## Installation

The singlecelldata package can be easily installed using the following command:

`pip install singlecelldata`

## Documentation

The SingleCell class reference manual can be found [here](https://github.com/edwinv87/singlecelldata/blob/master/docs/_build/latex/singlecelldata.pdf)

## Usage

### Basic Usage

The SingleCell class can be used to create an object which stores single-cell gene expression data and additional data about genes and cells in their respective dataframes. To create a SingleCell object sc, the following python code can be used:

```python
import pandas as pd
from singlecelldata import SingleCell

dataset = 'biase'

data_path = "data/" + dataset + '/' + dataset + "_data.csv"
celldata_path = "data/" + dataset + '/' + dataset + "_celldata.csv"
genedata_path = "data/" + dataset + '/' + dataset + "_genedata.csv"

# Create pandas dataframes by reading data from files
data = pd.read_csv(data_path, index_col=0)
celldata = pd.read_csv(celldata_path, index_col=0)
genedata = pd.read_csv(genedata_path, index_col = 0)

# Create a single cell object
sc = SingleCell(dataset, data, celldata, genedata)
```

In the above example, a SingleCell object, sc, was ceated by passing the dataset name and the main data, the cell data and gene data as Pandas dataframes. [Pandas](https://pandas.pydata.org/) is a powerpul python library for creating data structures from a variety of sources. Pandas can open and read data from numerous differernt file types such as csv files and creating dataframes from it. This enables the user to create SingleCell objects from a variety of different data file types.

### Detailed Usage

Detailed usage information can be found at

## Contact

Contact the author on vans.edw@gmail.com to give feedback/suggestions for further improvements and to report issues.
