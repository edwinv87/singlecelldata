# SingleCell

SingleCell is a python class available in the singlecelldata package for managing single-cell RNA-seq data. It contains three pandas dataframes; `data` for holding gene expression values (counts/normalized counts), `genedata` for holding more information about the genes e.g., gene names, and `celldata` which contains more information about cells such as cell types, labels etc.

## Installation

The singlecelldata package can be easily installed using the following command:

`pip install singlecelldata`

## Documentation

The SingleCell class reference manual can be found [here](https://singlecelldata.readthedocs.io/en/latest/index.html)

## Using the SingleCell class

### Basic Example

The SingleCell class can be used to create an object which stores single-cell gene expression data and additional data about genes and cells in their respective dataframes. To create a SingleCell object sc, the following python code can be used:

```python
# Import SingleCell
from singlecelldata import SingleCell

# Create the path and the dataset variables
path = '../data/293t_jurkat/'
dataset = '293t'

# Create the SingleCell object
sc = SingleCell(path, dataset, "10X")
sc.print()
```

In the above example, a SingleCell object, sc, was ceated by passing the dataset name and the main data, the cell data and gene data as Pandas dataframes. [Pandas](https://pandas.pydata.org/) is a powerpul python library for creating data structures from a variety of sources. Pandas can open and read data from numerous differernt file types such as csv files and creating dataframes from it. This enables the user to create SingleCell objects from a variety of different data file types.

### Detailed Example

More detailed example can be found [here](https://anaconda.org/edwinvans/creating_singlecell_object)

## Contact

Contact the author on <info@edwinvans.com> to give feedback/suggestions for further improvements and to report issues.
