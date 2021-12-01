# COGplot
**遺伝子のアミノ酸データをCOG分類して棒グラフとベン図を出力するスクリプト**

## Usage
```
$ python3 COGplot.py -AA  [genes1.fasta [genes2.fasta ...]]

#rpsblastの結果を入力データとして扱いたい場合
$ python3 COGplot.py -AA  [genes1.txt [genes2.txt ...]]
```

## 出力例
![](./images/COG_count.png)
![](./images/COG_ratio.png)
![](./images/venn3Diagram.png)
![](./images/COGvenn3Diagrams.png)
