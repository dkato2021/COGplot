# COGplot
**遺伝子のアミノ酸データをCOG分類して棒グラフとベン図、PCAの図を出力するスクリプト**
- rpsblastを実行する際のe-valueは偽遺伝子がカウントされることを防ぐためにデフォルトの値を1e-25にしています。
- ベン図を出力できる入力データ数の上限は6です。
- 入力データ数が４以上の場合はベン図の積集合の要素数とその面積は一致しません。
- PCAの入力にはCOG分類結果の比率データをCLR変換したものを使用しています。
- png -> pdfにしたい

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
  -rps [RPS [RPS ...]]  path to your results of rpsblast
  -AA [AA [AA ...]]     paths　to your amino acids files of genes(Venn diagram is not output if there are 6 or more files)
  -e EVALUE             evalue in rpsblast(default:1e-25)
  -cogdb COGDB          path to your cogdb to run rpsblast(default/home/tmp/db/COG/Cog)
  -cddid CDDID          path to your cddid_COG.tbl(default:/home/tmp/db/COG/cdd2cog/cddid_COG.tbl)
  -cog COG              path to your cog-20.def.tsv(default:/home/tmp/db/COG/cdd2cog/cog-20.def.tsv)
```
## 出力例
|  GCF_000024905.1_ASM2490v1_translated_cds  |  CDD  |COG  |Group  |
| ---- | ---- |
|  NC_013522.1_prot_WP_164925053.1_1  |  CDD:223666  | COG0593  | L |
|  NC_013522.1_prot_WP_012868759.1_2  |  CDD:223665  | COG0592  | L  |

![image](https://user-images.githubusercontent.com/78598272/144970110-b199816e-fec0-471f-891a-18ea0daaf967.png)

![](./images/COG_count.png)
![](./images/COG_ratio.png)
![](./images/venn3Diagram.png)
![](./images/COGvenn3Diagrams.png)
![](./images/1.png)
![](./images/2.png)
![](./images/3.png)
![](./images/4.png)
