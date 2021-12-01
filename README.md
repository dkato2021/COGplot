# COGplot
**遺伝子のアミノ酸データをCOG分類して棒グラフとベン図を出力するスクリプト**

## Usage
```
$ python3 COGplot.py -AA  [genes1.fasta [genes2.fasta ...]]

#rpsblastの結果を入力データとして扱いたい場合
$ python3 COGplot.py -AA  [genes1.txt [genes2.txt ...]]
```

## 補遺
- rpsblastを実行する際のe-valueは偽遺伝子がカウントされることを防ぐためにデフォルトの値を1e-15にしています。

- ベン図を出力できる入力データ数の上限は6です。

- 入力データ数が４以上の場合はベン図の積集合の要素数とその面積を一致させることができていません。

## 出力例
![](./images/COG_count.png)
![](./images/COG_ratio.png)
![](./images/venn3Diagram.png)
![](./images/COGvenn3Diagrams.png)
![](./images/1.png)
![](./images/2.png)
![](./images/3.png)
![](./images/4.png)
