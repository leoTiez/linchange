# LinChange

Investigating the linear correlation between two variables is a standard procedure
for analysing genomic data. This helper script in Python provides a command line
interface for producing pretty plots for the linear dependency.

## Requirements and Installation
The script is written in Python3 (Python 3.6 or larger is recommended). Moreover, the 
necessary libraries are installed with pip. Run

```commandline
python3 -m pip install requirements.txt
```

## Usage
Run
```commandline
python3 linchange.py [data_paths] --annotation=path/to/bed/or/gtf [-n name1 -n name2]  --title=Title --distinguish_strand --percentile=100 --save_fig --save_path=path/to/saved/figure.png
```

with
- `data_paths`: 2 or 4 data paths to the bigwig files. If 4 paths are passed, the first two and the second two 
are added together (e.g. + strand and - strand of WT and mutant).
- `--bed`: path to the annotation bed file.
- `-n`: that are printed on the plots axes. The first names is used for the x-axis, whilst
the second name is used for the y-axis.
- `--title` or `-t`: plot title. 
- `--distinguish_strand`: if set, two consecutive strands are considered as + and - (in that order)
- `--percentile`: Only keep the lower percentile.
- `--save_fig`: If set, figure is set to file at save_path
- `--save_path`: Figures is saved under this file name if save_fig is set.


