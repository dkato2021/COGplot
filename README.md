# COGplot
**遺伝子のアミノ酸データをCOG分類して棒グラフとベン図を出力するスクリプト**
- rpsblastを実行する際のe-valueは偽遺伝子がカウントされることを防ぐためにデフォルトの値を1e-25にしています。
- ベン図を出力できる入力データ数の上限は6です。
- 入力データ数が４以上の場合はベン図の積集合の要素数とその面積を一致させることができていません。
- CLR -> PCAの図も出力したい
from skbio.stats.composition import clr
from sklearn.decomposition import PCA

**依存**
- matplotlib-venn
- scikit-bio
## Usage
```
$ python3 COGplot.py -AA  [genes1.fasta [genes2.fasta ...]]

#rpsblastの結果を入力データとして扱いたい場合
$ python3 COGplot_rpsblast.py -rps  [genes1.txt [genes2.txt ...]]
```

```
def CLR_PCA(df = None):#各行にCOG
    #assert len(ratio_data.index)==4, print(1)
    df_clr = clr(df.T).T
    pca = PCA(n_components=2)
    _ = pca.fit_transform(df_clr.T)
    df_pca = pd.DataFrame(_, columns = ["PCA1", "PCA2"])
    def plot_PCA(df_pca, pca):
        fig = plt.figure(figsize=(6, 6))
        plt.scatter(df_pca.PCA1, df_pca.PCA2, alpha=0.8)
        #for x, y, name in zip(pca.components_[0], pca.components_[1], df.index):
        #    plt.text(x, y, name)
        #plt.scatter(pca.components_[0], pca.components_[1], alpha=0.8)
        plt.grid()
        plt.xlabel(f"PC1({(pca.explained_variance_ratio_[0]*100).round(2)}%)")
        plt.ylabel(f"PC2({(pca.explained_variance_ratio_[1]*100).round(2)}%)")
        plt.show()
    plot_PCA(df_pca, pca)
    return fig, df_pca
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
