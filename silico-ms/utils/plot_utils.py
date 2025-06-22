import pandas as pd  
import anndata as ad
import matplotlib.pyplot as plt   
import seaborn as sns 


def load_results_data(data_file: str,
                      group_file: str) -> ad.AnnData:
    data_df = pd.read_csv(data_file)
    group_df = pd.read_csv(group_file)
    
    vars = pd.DataFrame(data_df["compound_name"])
    X = data_df.drop(["compound_name"], axis=1).to_numpy().T
    adata = ad.AnnData(X, var=vars, obs=group_df) # X (n_obs * n_var)
    return adata





def plot_network(adata: ad.AnnData):
    pass


def plot_stack_bar_lipid(adata: ad.AnnData):
    compounds = list(adata.var["compound_name"])

    adata.var["lipid_class"] = [compound.split(' ')[0] for compound in compounds]
    df = pd.DataFrame(adata.X, columns=adata.var["compound_name"].to_list(), index=adata.obs["sample"].to_list()).T
    df = df.reset_index(names="compound_name")
    df = pd.melt(df, id_vars=["compound_name"], var_name="sample", value_name="intensity")
    df = pd.merge(df, adata.var, how="inner", on=["compound_name"])
    #df = pd.merge(df, adata.obs, how="inner", on=["sample"])

    df = df.groupby(["sample", "lipid_class"])["intensity"].sum().reset_index()
    df["intensity"] = df["intensity"] / df.groupby(["sample"])["intensity"].transform("sum") # normalization 

    df = pd.pivot(df, index='sample', columns='lipid_class',values='intensity')
    ax = df.plot(kind='bar', stacked=True)
    ax.legend(loc="center right", bbox_transform=ax.transAxes)
    plt.show()
    return ax

def plot_heatmap(adata: ad.AnnData):
    df = pd.DataFrame(adata.X, columns=adata.var["compound_name"].to_list(), index=adata.obs["sample"].to_list()).T
    sns.heatmap(data=df,square=True) 


def plot_barplot(data, 
            x = None, 
            y = None, 
            hue = None,
            title = '',
            errorbar = ('ci', 0.95), 
            capsize = 0.1,
            save_fname = None,
            dpi = 400):
    plt.title(title)
    sns.barplot(data=data, x=x, y=y, hue=hue, errorbar=errorbar, capsize=capsize)
    if save_fname is not None:
        plt.savefig(save_fname, dpi=dpi)
    plt.show()    


def plot_boxplot(data, x, y, hue, title = '', fliersize = 1, dpi = 300, save_fname = None):
    plt.figure(dpi=dpi)
    plt.title(title)
    ax = sns.boxplot(x=x, y=y, hue=hue, data=data, fliersize=fliersize)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.3, 1))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    if save_fname is not None:
        plt.savefig(save_fname, dpi=dpi)
    plt.show()
    #return ax 



def plot_composition_heatmap(data,
                             title = "",
                             save_fname = None,
                             dpi = 400):
    
    plt.figure(figsize=(10, 10))
    plt.title(title)
    ax = sns.heatmap(data, cmap="Blues", annot=True, 
                     fmt=".2f",linewidths=0.5, linecolor='black')
    ax.set(xlabel="", ylabel="")
    #ax.xaxis.tick_top()
    ax.spines["bottom"].set_visible(True)
    ax.spines["top"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["right"].set_visible(True)
    if save_fname is not None:
        plt.savefig(save_fname, dpi=dpi)
    plt.show()


def plot_fc_barplot(data,
                    title = "",
                    save_fname = None,
                    dpi = 400):
    plt.figure(figsize=(10, 10))
    plt.title(title)
    ax = sns.barplot(data = data, y ="full_fatty_acid", x = "log2FC")
    ax.set(xlabel="log2FC", ylabel="")
    ax.axvline(x=0, color="black", linewidth=1.0)
    if save_fname is not None:
        plt.savefig(save_fname, dpi=dpi)
    plt.show()


def plot_relative_composition_barplot(data,
                                      title = "",
                                      save_fname = None,
                                      dpi = 400):
    
    from plotnine import ggplot, geom_bar, aes, ggsave, ggtitle
    ax = (
        ggplot(data, aes("chain_level", fill="double_bond"))
        + geom_bar(position="fill") +
        ggtitle(title)
    )
    if save_fname is not None:
        ggsave(ax, filename=save_fname, dpi=dpi)
    return ax


def plot_valcanoplot(data,
                     title = "",
                     fc_threshold = 1.2,
                     pval_threshold = 0.05,
                     save_fname = None,
                     dpi = 400):
    import math 
    log2_fc_thres = math.log2(fc_threshold)
    log10_pval_thres = -math.log10(pval_threshold)
    
    plt.figure(figsize=(10, 10))
    plt.title(title)
    ax = sns.scatterplot(data=data, x="log2_fold_change", y = "log10_fdr", hue="significant")
    ax.axhline(y=log10_pval_thres, color='grey', linestyle='--')
    ax.axvline(x=-log2_fc_thres, color='grey', linestyle='--')
    ax.axvline(x=log2_fc_thres, color='grey', linestyle='--')
    if save_fname is not None:
        plt.savefig(save_fname, dpi=dpi)
    plt.show()