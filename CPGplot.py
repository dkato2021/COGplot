import os, sys, argparse, warnings, csv
sys.path.append('./') ;warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib_venn import venn2
from matplotlib_venn import venn3
from collections import Counter
from MyLibrary import venn

def get_args():
    parser = argparse.ArgumentParser(description='dkato. November, 2021')
    #parser.add_argument('-rps' , dest ='rps', nargs='*',
    #                    help = 'path_to_rpsRes', required=True)
    parser.add_argument('-AA' , dest ='AA', nargs='*',
                        help = 'paths　to your amino acid file of genes', required=True)  
    parser.add_argument('-e' , dest ='evalue',
                        default= 1e-25, 
                        help = 'evalue in rpsblast(default:1e-25)')
    parser.add_argument('-cogdb' , dest ='cogdb',
                        default= '/home/tmp/db/COG/Cog', 
                       help = 'path to your cogdb(default/home/tmp/db/COG/Cog:)')    
    parser.add_argument('-cddid' , dest ='cddid',
                        default= '/home/tmp/db/COG/cdd2cog/cddid_COG.tbl',
                        help = 'path to your cddid_COG.tbl(default:/home/tmp/db/COG/Cog)')
    parser.add_argument('-cog', dest='cog',
                        default='/home/tmp/db/COG/cdd2cog/cog-20.def.tsv',
                        help = 'path to your cog-20.def.tsv(default:/home/tmp/db/COG/Cog)')
    return parser.parse_args()

  def run_rpsblast(paths_to_proteins = None, 
                 path_to_cogdb = None, 
                 evalue = None):
    
    error1 = "specify the path to your Cog database with cogdb option. (default:/home/tmp/db/COG/Cog)"
    #assert os.path.exists('/home/tmp/db/COG/Cog/'), error1
    
    if 'res_rpsblast' not in os.listdir(path='./'):
        os.system('mkdir res_rpsblast')
    
    path_to_rpsRes = []
    for path in paths_to_proteins:
        name = os.path.splitext(os.path.basename(path))[0]
        subprocess.run(f"rpsblast -query {path} -db {path_to_cogdb} -out ./res_rpsblast/{name}_rpsblastout.txt -evalue {evalue} -outfmt 6"
                       , shell=True)
        path_to_rpsRes.append(f"./res_rpsblast/{name}_rpsblastout.txt")
    return path_to_rpsRes
  
  def preprocess(rps = None,
               cddid = None,
               cog = None):

    cddid["CDD"] = "CDD:" + cddid["CDD"].astype(str)
    _ = pd.merge(rps, cddid, on = ["CDD"]).iloc[:, [0, 1, 12]]
    _df = pd.merge(_, cog, on = ["COG"]).iloc[:, [0, 1, 2, 3]]
    return _df, dict(Counter("".join(_df['Group'])))

def sorter(df_i = None,
           A2Z = None):
    tmp = []
    for i in range(len(A2Z)):
        if A2Z[i] in list(df_i.keys()):
            tmp += [df_i[A2Z[i]]]
        else:
            tmp += [0]
    return tmp

def get_main_dataset(path_to_rpsRes = None,
                     path_to_cddid = None,
                     path_to_cog = None):
    if 'out' not in os.listdir(path='./'):
        os.system('mkdir ./out/')
    if 'COGdata' not in os.listdir(path='./out/'):
        os.system('mkdir ./out/COGdata/') 
    
    df_i = {}
    out_COG_i = {}
    out = pd.DataFrame()
    out_COG = pd.DataFrame()
    A2Z = [chr(i) for i in range(65, 65+26)]
    for i, path in enumerate(path_to_rpsRes):
        #load data
        rps_i = pd.read_table(path_to_rpsRes[i], names=["cdd_id", "CDD", "a", "b","c","d",
                                        "e","f","g","h","eValue","j"]).drop_duplicates(['cdd_id'])

        cddid = pd.read_table(path_to_cddid, names=["CDD", "COG", "a", "b", "c"])

        cog = pd.read_table(path_to_cog, names=["COG", "Group", "C3",
                                         "D3", "E3", "F3", "G3"], encoding='cp1252')
        
        col_name = os.path.splitext(os.path.basename(path))[0]#[:-len('_rpsblastout')]
        #processing
        out_COG_i[col_name], df_i[col_name] = preprocess(rps = rps_i, 
                                  cddid = cddid,
                                  cog = cog)
        
        
        COG_i = out_COG_i[col_name].rename(columns={'cdd_id': f"{col_name}"}).iloc[:, [0, 3]]
        out_i = pd.DataFrame(sorter(df_i = df_i[col_name], A2Z = A2Z), columns=[f"{col_name}"])

        out = pd.concat([out, out_i], axis = 1)
        out_COG = pd.concat([out_COG, COG_i], axis = 1)
    count_data = pd.concat([pd.DataFrame(A2Z, columns=['COG']), out], axis = 1)
    

    _count_data = count_data.iloc[:,1:len(count_data.columns)]
    _ = _count_data/_count_data.sum()
 
    ratio_data = pd.concat([pd.DataFrame(A2Z, columns=['COG']), _], axis = 1)
    count_data.to_csv("./out/COGdata/COG_count.csv") ;ratio_data.to_csv("./out/COGdata/COG_ratio.csv")
    out_COG.to_csv("./out/COGdata/COG_annotation.csv")
    
    return count_data, ratio_data, out_COG_i

def plot_bar(df = None, name = None):
    # 棒の配置位置、ラベルを用意
    labels = list(df['COG'])
    x = np.array(range(len(labels)))

    # 各系列のデータを用意
    data, legend = [], []
    for col_name in df.columns[1:len(df.columns)]:
        data.append(df[col_name])
        legend.append(col_name) 

    # マージンを設定
    margin = 0.2  #0 <margin< 1
    totoal_width = 1 - margin
    fig = plt.figure(figsize=(15,10))
    # 棒グラフをプロット
    for i, h in enumerate(data):
        pos = x - totoal_width *( 1- (2*i+1)/len(data) )/2
        plt.bar(pos, h, width = totoal_width/len(data))

    plt.legend(legend)
    plt.xticks(x, labels)
    #plt.show()
    fig.savefig(f"./out/COG_{name}.png")

def plot_or_not(unique_COGs):
    return sum([unique_COG==set() for unique_COG in unique_COGs]) !=len(unique_COGs)

def venn_func(unique_COG, labels, ax):
    subsets = venn.get_labels(unique_COG, fill=['number', 'logic'])
    if len(list(dataset.keys()))==2:
        return venn2(subsets=unique_COG, set_labels = labels)
    elif len(list(dataset.keys()))==3:
        return venn3(subsets=unique_COG, set_labels = labels)
    elif len(list(dataset.keys()))==4:
        return venn.venn4(subsets, ax, names = labels)
    elif len(list(dataset.keys()))==5:
        return venn.venn5(subsets, ax, names = labels)
    elif len(list(dataset.keys()))==6:
        return venn.venn6(subsets, ax, names = labels)
    else:
        sys.exit()
            
def plot_venn(dataset = None):
    A2Z = [chr(i) for i in range(65, 65+26)]

    #1
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, aspect='equal')

    unique_COG = []
    for j in range(len(dataset.keys())):
        x = dataset[list(dataset.keys())[j]]
        unique_COG.append(set(x['COG'].unique()))

    venn_func(unique_COG, list(dataset.keys()), ax)
    ax.set_title('All genes')  
    fig.savefig(f"./out/venn{len(list(dataset.keys()))}Diagrams.png")

    #2
    fig = plt.figure(figsize=(30,40))
    for i, alphabet in enumerate(A2Z):
        unique_COG = []
        for j in range(len(dataset.keys())):
            x = dataset[list(dataset.keys())[j]]
            unique_COG.append(set(x[x['Group']==f'{alphabet}']['COG'].unique()))

        if plot_or_not(unique_COG):
            ax = fig.add_subplot(6, 5, i+1)
            venn_func(unique_COG, list(dataset.keys()), ax)
            ax.set_title(f'{alphabet}')
            plt.tight_layout()
            fig.savefig(f"./out/COGvenn{len(list(dataset.keys()))}Diagrams.png")

            
if __name__ == "__main__":
    print('1.rpsblast now..')
    path_to_rpsRes = run_rpsblast(paths_to_proteins = get_args().AA, 
                                  path_to_cogdb = get_args().cogdb, 
                                 evalue = get_args().evalue)
    print('2.creating barplot..')
    count_data, ratio_data, dataset = get_main_dataset(path_to_rpsRes = path_to_rpsRes,
                                                      path_to_cddid = get_args().cddid,
                                                      path_to_cog = get_args().cog)

    plot_bar(df = count_data, name ='count')
    plot_bar(df = ratio_data, name ='ratio')
    print(f'==>COG_count.png and COG_ratio.png are created.')
    print('3.creating venn diagrams..')
    plot_venn(dataset = dataset)
    print(f'==>venn diagrams are created.')


