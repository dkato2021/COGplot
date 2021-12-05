# COGplot
**遺伝子のアミノ酸データをCOG分類して棒グラフとベン図を出力するスクリプト**
- rpsblastを実行する際のe-valueは偽遺伝子がカウントされることを防ぐためにデフォルトの値を1e-25にしています。
- ベン図を出力できる入力データ数の上限は6です。
- 入力データ数が４以上の場合はベン図の積集合の要素数とその面積を一致させることができていません。
- 比率データのCLR -> PCAの図も出力します。


**依存**
- matplotlib-venn

## Usage
```
$ python3 COGplot.py -AA  [genes1.fasta [genes2.fasta ...]]

#rpsblastの結果を入力データとして扱いたい場合
$ python3 COGplot.py -rps  [genes1.txt [genes2.txt ...]]
```
## optional arguments
```
-h, --help            show this help message and exit
-rps [RPS [RPS ...]]  path_to_rpsRes
-AA [AA [AA ...]]     paths　to your amino acid file of genes(Venn diagram is
                      not output if there are 6 or more files)
-e EVALUE             evalue in rpsblast(default:1e-25)
-cogdb COGDB          path to your cogdb(default/home/tmp/db/COG/Cog:)
-cddid CDDID          path to your
                      cddid_COG.tbl(default:/home/tmp/db/COG/Cog)
-cog COG              path to your
                      cog-20.def.tsv(default:/home/tmp/db/COG/Cog)
```
## 出力例
![](./images/COG_count.png)
![](./images/COG_ratio.png)
![](./images/venn3Diagram.png)
![](./images/COGvenn3Diagrams.png)
![](./images/1.png)
![](./images/2.png)
![](./images/3.png)
![](./images/4.png)
