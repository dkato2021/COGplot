# COGplot
**dependency**
- matplotlib-venn  

## Installation
```
$ pip install matplotlib-venn
$ git clone https://github.com/dkato2021/COGplot.git
$ chmod u+x *.py
```
## Usage
```
$ Lossplot.py -AA gene1.fasta gene2.fasta ...
$ COGplot.py -AA gene1.fasta gene2.fasta ...
```

## tips
```
$ COGplot.py -AA ./X/* ./Y/* -e 1e-4 1e-12 1e-20
```
## How to determine the Evalue
![](./_/lossver18.png)
![](./_/anime1.png)
![](./_/animeTA.gif)
![](./_/LossGraphA.png)
![](./_/allPCA_ratioA.gif)
![](./_/LossGraphY.png)
![](./_/allPCA_ratioY.gif)

![](./_/LossGraphTver3.png)
![](./_/animeT.gif)
## How to detect unique genes
![](./_/unique_ver2.png)
## How to interpret the PCA diagram
![](./_/X.png)
Reference
- https://statistics.co.jp/reference/software_R/statR_9_principal.pdf
## Output Example
![](./_/COG_count.png)
![](./_/COG_ratio.png)
![](./_/venn3Diagram.png)
![](./_/COGvenn3Diagrams.png)
![](./_/w.png)
![](./_/q.png)
![](./_/1.png)
![](./_/3.png)

