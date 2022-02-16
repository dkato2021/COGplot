#!/usr/bin/env python3
import os, sys, argparse, warnings, csv
warnings.filterwarnings('ignore')
from subprocess import Popen
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sklearn.decomposition import PCA
import math
from itertools import chain
from collections import Counter, Iterable


def get_args():
    parser = argparse.ArgumentParser(description='dkato. Feb. 2022')
    parser.add_argument('-AA' , dest ='AA', nargs='*',
                        help = 'paths　to your amino acids files of genes')
    parser.add_argument('-rl' , dest ='rpsloss', nargs='*',
                        help = 'specify ./rps_forLoss/* when you omit rpsblast')
    parser.add_argument('-df' , dest ='df', 
                        help = 'specify ./LossGraph.csv')
    parser.add_argument('-t', dest='num_threads',
                       default=42,type = int, help = 'num_threads(default:42)')        
    parser.add_argument('-l' , dest ='loss_size',
                        default= 6, type = int, help = 'specify a integer value: graph size of loss graph(default:6)')
    parser.add_argument('-p' , dest ='points',
                        default= 150, type = int, help = 'number of points in loss graph')
    parser.add_argument('-d' , dest ='delta',
                        default= 4, type = int, help = 'search interval of E-value(value 4 is recommended)') 
    
    parser.add_argument('-PCA' , dest ='PCA_size',
                        default= 5, type = int, help = 'specify a integer value: graph size of PCA plot(default:5)')
    parser.add_argument('-o', dest='n_orange',
                        default=0,type = int, help = 'Number of points dyed in orange in a PCA plot(default:0)')
    
    parser.add_argument('-cogdb' , dest ='cogdb',
                        default= '/home/tmp/db/COG/Cog', 
                       help = 'path to your cogdb to run rpsblast(default:/home/tmp/db/COG/Cog)')    
    parser.add_argument('-cddid' , dest ='cddid',
                        default= '/home/tmp/db/COG/cdd2cog/cddid_COG.tbl',
                        help = 'path to your cddid_COG.tbl(default:/home/tmp/db/COG/cdd2cog/cddid_COG.tbl)')
    parser.add_argument('-cog', dest='cog',
                        default='/home/tmp/db/COG/cdd2cog/cog-20.def.tsv',
                        help = 'path to your cog-20.def.tsv(default:/home/tmp/db/COG/cdd2cog/cog-20.def.tsv)')


    return parser.parse_args()
#'/Users/daiki/Python/M2/rpsblast/data/cddid_COG.tbl',
#'/home/tmp/db/COG/cdd2cog/cddid_COG.tbl'
#'/Users/daiki/Python/M2/rpsblast/data/cog-20.def.tsv',
#'/home/tmp/db/COG/cdd2cog/cog-20.def.tsv'

def run_rpsblast_ForLossGraph(paths_to_proteins = None, 
                 path_to_cogdb = None, num_threads = None):
    
    def split_list(lst, n):  
        for i in range(0, len(lst), n): 
            yield lst[i:i + n] 
           
    procs_ForLoss = []
    path_to_rpsRes_ForLoss = []
    if f'rps_forLoss' not in os.listdir(path='./'):
        os.system(f'mkdir rps_forLoss')
    
    block_paths = list(split_list(paths_to_proteins, num_threads))
    for block_path in block_paths:
        for path in block_path:
            name = os.path.splitext(os.path.basename(path))[0]
            procs_ForLoss += [Popen(f"rpsblast -query {path} -db {path_to_cogdb} -out ./rps_forLoss/{name}.txt -evalue .1 -outfmt 6"
                               , shell=True)]
            path_to_rpsRes_ForLoss.append(f"./rps_forLoss/{name}.txt")

        [p.wait() for p in procs_ForLoss]
    return path_to_rpsRes_ForLoss


def preprocess(rps = None,
               cddid = None,
               cog = None):

    cddid["CDD"] = "CDD:" + cddid["CDD"].astype(str)
    _ = pd.merge(rps, cddid, on = ["CDD"]).iloc[:, [0, 1, 12]]
    _df = pd.merge(_, cog, on = ["COG"]).iloc[:, [0, 1, 2, 3, 4, 5]]
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

def get_rps_i_forLoss(path_to_rpsRes_i, evalue = None):
    tmp_i = pd.read_table(path_to_rpsRes_i, names=["cdd_id", "CDD", "a", "b","c","d",
                                        "e","f","g","h","eValue","j"])
    return tmp_i[tmp_i.eValue <= evalue].drop_duplicates(['cdd_id'])

def get_main_dataset(path_to_rpsRes = None,
                     path_to_cddid = None,
                     path_to_cog = None,  Log = None, evalue = None):

    df_i = {}
    out_COG_i = {}
    out = pd.DataFrame()
    out_COG = pd.DataFrame()
    A2Z = [chr(i) for i in range(65, 65+26)]
    for i, path in enumerate(path_to_rpsRes):
        #load data

        rps_i = get_rps_i_forLoss(path_to_rpsRes[i], evalue =  evalue)

        cddid = pd.read_table(path_to_cddid, names=["CDD", "COG", "a", "b", "c"])

        cog = pd.read_table(path_to_cog, names=["COG", "Group", "gene_name",
                                         "gene", "E3", "F3", "G3"], encoding='cp1252')
        
        col_name = os.path.splitext(os.path.basename(path))[0]
        #processing
        out_COG_i[col_name], df_i[col_name] = preprocess(rps = rps_i, 
                                  cddid = cddid,
                                  cog = cog)
        
        COG_i = out_COG_i[col_name].rename(columns={'cdd_id': f"{col_name}"}).iloc[:, [0,1,2,3,4,5]]
        out_i = pd.DataFrame(sorter(df_i = df_i[col_name], A2Z = A2Z), columns=[f"{col_name}"])

        out = pd.concat([out, out_i], axis = 1)
        out_COG = pd.concat([out_COG, COG_i], axis = 1)
    count_data = pd.concat([pd.DataFrame(A2Z, columns=['COG']), out], axis = 1)
    

    _count_data = count_data.iloc[:,1:len(count_data.columns)]
    _ = _count_data/_count_data.sum()
 
    ratio_data = pd.concat([pd.DataFrame(A2Z, columns=['COG']), _], axis = 1)
    if Log:
        if f'out_{evalue}' not in os.listdir(path='./'):
            os.system(f'mkdir ./out_{evalue}/')
        if 'COGdata' not in os.listdir(path=f"./out_{evalue}/"):
            os.system(f'mkdir ./out_{evalue}/COGdata/') 
        #count_data.to_csv(f"./out_{evalue}/COGdata/COG_count.csv") 
        #ratio_data.to_csv(f"./out_{evalue}/COGdata/COG_ratio.csv")
        #out_COG.to_csv(f"./out_{evalue}/COGdata/COG_annotation.csv")

    return count_data, ratio_data, out_COG_i
 
def get_loss_data(path_to_rpsRes_forLoss = None,
              evalue = None,
              delta = None,
              points = None,
              size = None, 
              path_to_cddid = None,
              path_to_cog = None):

    _k={}
    ratio_data={}
    for i, e in enumerate(evalue):
        _k[i], ratio_data[i] = get_main_dataset(path_to_rpsRes =  path_to_rpsRes_forLoss,
                                 path_to_cddid = path_to_cddid,
                                 path_to_cog = path_to_cog,
                                 Log = False,
                                 evalue = float(e))[:2]
    def get_bray(k, k_dash):
        return (abs(k - k_dash).sum()+1)/(k + k_dash).sum()
    
    loss=[]
    for i in range(delta, len(evalue)-delta):
        _k_set_i = [_k[i-delta], _k[i], _k[i+delta]]
        k_minus  = np.array(_k_set_i[0].iloc[:,1:])
        k        = np.array(_k_set_i[1].iloc[:,1:])
        k_plus   = np.array(_k_set_i[2].iloc[:,1:])
        loss    += [math.log((len(path_to_rpsRes_forLoss)/(k.sum()+1))*(get_bray(k, k_minus)*get_bray(k, k_plus)))]

    if f'LossGraph' not in os.listdir(path='./'):
        os.system(f'mkdir ./LossGraph/')
    e  = [float(_) for _ in evalue]
    df = pd.DataFrame(zip(e[delta:-delta], loss), columns = ['Evalue', 'Loss'])
    #optimal_e = format(df[df.Loss == min(loss)].iloc[0,0],'.0e')
    df.to_csv(f"./LossGraph/LossGraph_delta{delta}_points{points}.csv")
    return df, ratio_data


def plot_loss(df, delta = None, points = None, size = None):
    optimal_e = format(df[df.Loss == min(df.Loss)].iloc[0,0],'.0e')
    fig = plt.figure(figsize=(4*size,2*size))
    ax1 = fig.subplots()
    ax1.set_title(f'Optimal E-value is {optimal_e}', fontsize= 4*size)
    ax1.set_xscale('log'); ax1.grid(which="both")
    ax1.set_xlabel('E-value', size = 4*size)
    plt.xticks(fontsize= 4*size); plt.yticks(fontsize= 3*size)
    ax1.set_ylabel('Loss', size = 4*size)
    plt.plot(df.Evalue, df.Loss, marker="D", c="m")
    
    fig.savefig(f"./LossGraph/LossGraph_delta{delta}_points{points}.pdf") 
    
def CLR_PCA(df = None, size = None, delta = None, tag = None, n_orange = None, CLR = None, evalue = None):#各行にCOG。
    if f'allPCA_{tag}' not in os.listdir(path='./'):
        os.system(f'mkdir ./allPCA_{tag}/')
    if f'{evalue}' not in os.listdir(path=f"./allPCA_{tag}/"):
        os.system(f'mkdir ./allPCA_{tag}/{evalue}/') 
    #if f'PCA_{tag}_{delta}' not in os.listdir(path=f"./out_{evalue}/PCA_{tag}/"):
    #    os.system(f'mkdir ./out_{evalue}/PCA_{tag}/PCA_{tag}_{delta}') 
        
    def Myclr(df):
        def geo_mean(iterable):
            a = np.array(iterable).astype(float)
            return a.prod()**(1.0/len(a))
        df_clr = pd.DataFrame()
        for i in range(len(df.index)):
            tmp = df.iloc[i,:]/geo_mean(df.iloc[i,:])
            df_clr = pd.concat([df_clr, tmp.map(math.log)], axis=1)
        return df_clr.T
    #ゼロ値の補完
    clr_in = df.iloc[:, 1:] + 1
    
    #CLR
    if CLR:
        df_clr = Myclr(clr_in.T).T#clr関数は行方向に和が１ものしか受け付けない
    else:
        df_clr = clr_in
    
    #PCA
    pca = PCA(n_components=2, random_state=42)
    _ = pca.fit_transform(df_clr.T)
    df_pca = pd.DataFrame(_, columns = ["PCA1", "PCA2"])
    
    def plot_PCA(df_pca, pca, df, evalue):
        fig = plt.figure(figsize=(size *2, size * 2))
        ax1 = fig.subplots()
        ax1.scatter(df_pca.PCA1[:n_orange], df_pca.PCA2[:n_orange], alpha=0.8, c='darkorange')
        ax1.scatter(df_pca.PCA1[n_orange:], df_pca.PCA2[n_orange:], alpha=0.8)
        for x, y, name in zip(df_pca.PCA1, df_pca.PCA2, df.columns[1:]):
            ax1.text(x, y, name)
        
        ax1.grid()
        ax1.set_xlabel(f"PC1({(pca.explained_variance_ratio_[0]*100).round(2)}%)")
        ax1.set_ylabel(f"PC2({(pca.explained_variance_ratio_[1]*100).round(2)}%)")
        fig.savefig(f"./allPCA_{tag}/{evalue}/PCA_COG_{tag}_{evalue}.pdf")

        ax2 = ax1.twiny().twinx()
        for x, y, name in zip(pca.components_[0], pca.components_[1], df.COG):
            ax2.text(x, y, name)
            ax2.arrow(x=0,y=0, dx=x, dy=y,
                     width=.0001, length_includes_head=True,color='m')
        ax2.scatter(pca.components_[0],  pca.components_[1], alpha=0, color='m')
        ax1.set_title(f"PC1 Loading", fontsize=20/size*2, color='m')
        ax2.set_ylabel(f"PC2 Loading", fontsize=20/size*2, color='m')
        fig.savefig(f"./allPCA_{tag}/{evalue}/PCA_COG_{tag}_{evalue}_withLoadingFactor.pdf")

    plot_PCA(df_pca, pca, df, evalue)
    
    #コードが冗長
    def plot_PCA_NoName(df_pca, pca, df, evalue):
        fig = plt.figure(figsize=(size *2, size * 2))
        ax1 = fig.subplots()
        ax1.scatter(df_pca.PCA1[:n_orange], df_pca.PCA2[:n_orange], alpha=0.8, c='darkorange')
        ax1.scatter(df_pca.PCA1[n_orange:], df_pca.PCA2[n_orange:], alpha=0.8)
        ax1.grid()
        ax1.set_xlabel(f"PC1({(pca.explained_variance_ratio_[0]*100).round(2)}%)")
        ax1.set_ylabel(f"PC2({(pca.explained_variance_ratio_[1]*100).round(2)}%)")

        ax2 = ax1.twiny().twinx()
        for x, y, name in zip(pca.components_[0], pca.components_[1], df.COG):
            ax2.text(x, y, name)
            ax2.arrow(x=0,y=0, dx=x, dy=y,
                     width=.0001, length_includes_head=True,color='m')
        ax2.scatter(pca.components_[0],  pca.components_[1], alpha=0, color='m')
        ax1.set_title(f"PC1 Loading", fontsize=20/size*2, color='m')
        ax2.set_ylabel(f"PC2 Loading", fontsize=20/size*2, color='m')
        fig.savefig(f"./allPCA_{tag}/{evalue}/PCA_COG_{tag}_{evalue}_NoName.pdf")

    plot_PCA_NoName(df_pca, pca, df, evalue)
    
def main():
    num = get_args().points
    evalue_ForLoss = (np.ones(num+2)*float(1e-1) ** np.arange(num+2))[1:]
    
    if get_args().AA is not None:
        print(f'- rpsblast for loss graph..')
        path_to_rpsRes_ForLoss = run_rpsblast_ForLossGraph(paths_to_proteins = get_args().AA, 
                                                           path_to_cogdb = get_args().cogdb,
                                                           num_threads = get_args().num_threads)

        print('- loss graph..')
        df, ratio_data = get_loss_data(path_to_rpsRes_forLoss = path_to_rpsRes_ForLoss,
                                  evalue = evalue_ForLoss,
                                  points = get_args().points,
                                  size = get_args().loss_size, 
                                  delta = get_args().delta, 
                                  path_to_cddid = get_args().cddid,
                                  path_to_cog = get_args().cog)
        
        plot_loss(df, delta = get_args().delta, points = get_args().points, size = get_args().loss_size)
        print('- PCA..')
        e = df.Evalue
        for i in range(get_args().delta, len(ratio_data)-get_args().delta):
            CLR_PCA(df = ratio_data[i], size = get_args().PCA_size,
                    delta =1, tag = "ratio", n_orange = get_args().n_orange, CLR = True, 
                    evalue = format(e[i],'.0e'))
        
    elif get_args().df is not None:
        print('- loss graph..')
        df = pd.read_csv(get_args().df, index_col=0)
        plot_loss(df, delta = get_args().delta, points = get_args().points, size = get_args().loss_size)
        
    elif get_args().rpsloss is not None:
        print('- loss graph..')
        df, ratio_data = get_loss_data(path_to_rpsRes_forLoss = get_args().rpsloss,
                                  evalue = evalue_ForLoss,
                                  points = get_args().points,
                                  size = get_args().loss_size, 
                                  delta = get_args().delta, 
                                  path_to_cddid = get_args().cddid,
                                  path_to_cog = get_args().cog)
        
        plot_loss(df, delta = get_args().delta, points = get_args().points, size = get_args().loss_size)
        print('- PCA..')
        e = df.Evalue
        for i in range(get_args().delta, len(ratio_data)-get_args().delta):
            CLR_PCA(df = ratio_data[i], size = get_args().PCA_size,
                    delta =1, tag = "ratio", n_orange = get_args().n_orange, CLR = True, 
                    evalue = format(e[i-get_args().delta],'.0e'))
        
if __name__ == "__main__":
    main()







