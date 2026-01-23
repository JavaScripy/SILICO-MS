import os
import itertools
from typing import Tuple, List


import yaml
import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter


from silico_ms.pipeline_runner import LipidAnnotatorRunner


def dict_to_yaml(
    file: str,
    data: dict
) -> None:
    """Write a dictionary to a YAML file."""
    with open(file, "w", encoding="utf-8") as f:
        yaml.dump(data, f, allow_unicode=True, sort_keys=False)
    

def yaml_to_dict(file: str) -> dict:
    """Load a YAML file into a dictionary."""
    with open(file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


class HyperparamOptimizer:
    
    def __init__(
        self,
        config: dict,
        params: dict
    ) -> None:
        self.config = config
        self.params = params
    
    
    def optimize(
        self,
        feature_id_col: str ="feature_id", 
        score_col: str ="total_score",
        mode: str = "weights"
    ) -> pd.DataFrame:
        """"""
        
        match mode:
            case "weights":
                params_list = self._prepare_params_list_weights(params=self.params)
                df_metrics_list = self._optimize_base(
                                        params_list=params_list,
                                        feature_id_col=feature_id_col,
                                        score_col=score_col
                                    )
                df_metrics_total = self._process_total_metrics(
                                    df_list=df_metrics_list,
                                    params=params_list[0],
                                    mode=mode
                                )
            case "spec_sim_type":
                params_list = self._prepare_params_list_spec_sim_type(params=self.params)
                df_metrics_list = self._optimize_base(
                                        params_list=params_list,
                                        feature_id_col=feature_id_col,
                                        score_col=score_col
                                    )
                df_metrics_total = self._process_total_metrics(
                                    df_list=df_metrics_list,
                                    params=params_list[0],
                                    mode=mode
                                )
            case _:
                raise ValueError(f"Unknown hyperparams: {mode}")
        
        return df_metrics_total


    def _prepare_params_list_weights(
        self,
        params: dict
    ) -> List[dict]:
        """Prepare a list of parameter dictionaries for grid search."""
        weight_step = params["weight_step"]
        equalize_weights = params["equalize_weights"]
        n = max(1, int(np.round(1 / weight_step)))
        spec_weights = np.linspace(0, 1, n + 1) 
        
        if equalize_weights:
            rt_weights = (1 - spec_weights) / 2
            precursor_mz_weights = rt_weights
            df_weights = pd.DataFrame({
                "rt_weight": rt_weights,
                "precursor_mz_weight": precursor_mz_weights,
                "spec_weight": spec_weights
            })
            out_dir_suffix = "equal_param_set"
        else:
            rt_weights  = np.linspace(0, 1, n + 1)
            combinations = list(itertools.product(rt_weights, spec_weights))
            df_weights = pd.DataFrame(combinations, columns=["rt_weight", "spec_weight"])
            df_weights["precursor_mz_weight"] = 1 - (df_weights["rt_weight"] + df_weights["spec_weight"])        
            out_dir_suffix = "unequal_param_set"
        
        params_list = []
        for idx, row in df_weights.iterrows():
            params_copy = params.copy()
            #params_copy["rt_weight"] = float(rt_weights[idx])
            #params_copy["precursor_mz_weight"] = float(precursor_mz_weights[idx])
            #params_copy["spec_weight"] = float(spec_weights[idx])
            params_copy["rt_weight"] = float(row["rt_weight"])
            params_copy["precursor_mz_weight"] = float(row["precursor_mz_weight"])
            params_copy["spec_weight"] = float(row["spec_weight"])
            tmp_out_dir = os.path.join(
                                        params["base_out_dir"],
                                        f"{out_dir_suffix}_{idx+1}"
                                    )
            params_copy["out_dir"] = str(tmp_out_dir).replace("\\", "/")
            params_list.append(params_copy)
        
        return params_list
    
    def _prepare_params_list_spec_sim_type(
        self,
        params: dict
    ) -> List[dict]:
        """"""
        out_dir_suffix = "spec_similarity_type"
        spec_sim_types = params["ms2_spectrum_similarity_type"]
        
        params_list = []
        for spec_sim_type in spec_sim_types:
            params_copy = params.copy()
            params_copy["ms2_spectrum_similarity_type"] = spec_sim_type
            tmp_out_dir = os.path.join(
                                        params["base_out_dir"],
                                        f"{out_dir_suffix}_{spec_sim_type}"
                                    )
            params_copy["out_dir"] = str(tmp_out_dir).replace("\\", "/")
            params_list.append(params_copy)
        
        return params_list
    
    def _optimize_base(
        self,
        params_list: List[dict],
        feature_id_col="feature_id", 
        score_col="total_score"
    ) -> List[pd.DataFrame]:
        """Perform a grid search over hyperparameters."""
        df_results_total_list = self.annotate_features(params_list=params_list)
        
        zip_list = list(zip(df_results_total_list, params_list))
        df_metrics_list = [
            self._generate_metrics_dataframe(
                data=df_results, 
                params=params,
                feature_id_col=feature_id_col,
                score_col=score_col
            )
            for df_results, params in tqdm(zip_list, desc="Generate metrics dataframe ...")
        ]
    
        return df_metrics_list
    
    def annotate_features(
        self,
        params_list: List[dict]
    ) -> Tuple[List[pd.DataFrame], List[dict]]:
        """"""
        df_results_total_list = [
            self._run_lipid_annotator(
                config=self.config,
                params=params
            ) for params in tqdm(params_list, desc="Hyperparam Optimize")
        ]
        
        return df_results_total_list
    
    def _run_lipid_annotator(
        self,
        config: dict,
        params: dict
    ) -> pd.DataFrame:
        """"""
        runner = LipidAnnotatorRunner(config=config, params=params)
        df_results_total = runner.annotate_lipids()
        
        return df_results_total
    
    def _generate_metrics_dataframe(
        self,
        data: pd.DataFrame,
        params: dict,
        feature_id_col: str = "feature_id",
        score_col: str = "total_score"
    ) -> pd.DataFrame:
        """"""
        df_results = data.copy()
        df_summary = self._compute_metrics(
                        data=df_results, 
                        feature_id_col=feature_id_col, 
                        score_col=score_col
                    )
        df_summary = df_summary.reset_index(drop=True)
        
        df_params = pd.DataFrame({
            "work_dir": params["out_dir"],
            "rt_weight": params["rt_weight"],
            "precursor_mz_weight": params["precursor_mz_weight"],
            "spec_weight": params["spec_weight"],
            "spec_sim_type": params["ms2_spectrum_similarity_type"]
        }, index=[0])
        
        df_params_multi = pd.concat(
                            [df_params] * len(df_summary), 
                            axis=0, 
                            ignore_index=True
                        ).reset_index(drop=True)
        df_metrics = pd.concat([df_params_multi, df_summary], axis=1)
        
        # write dataframe
        metrics_file = os.path.join(params["out_dir"], "metrics.csv")
        df_metrics.to_csv(metrics_file, index=None)

        return df_metrics
    
    def _compute_metrics(
        self,
        data: pd.DataFrame, 
        feature_id_col: str = "feature_id",
        score_col: str = "score"
    ) -> pd.DataFrame:
        """"""
        df = data[[feature_id_col, score_col]].copy()
        df[feature_id_col] = df[feature_id_col].astype(str)
        df[score_col] = df[score_col].astype(float)
        df = df.sort_values([feature_id_col, score_col], ascending=[True, False])

        results = []    
        top1_score = np.nan
        top2_avg_score = np.nan
        top3_avg_score = np.nan
        top1_gap = np.nan
        top1_margin = np.nan
        top2_margin = np.nan
        top3_margin = np.nan
        H_norm = np.nan
            
        for feature_id, group in df.groupby(feature_id_col):
            scores = group[score_col]
            top1_score = scores.iloc[0]
            if len(group) >= 2:
                top2_score = scores.iloc[1]
                top2_avg_score = scores.iloc[:2].mean()
                top1_gap = top1_score - top2_score
                top1_margin = top1_score - scores.iloc[2:].mean()
                H_norm = self._compute_normalized_entropy(scores.dropna().values)
                if len(group) >= 3:
                    top3_avg_score = scores.iloc[:3].mean()
                    top2_margin = top2_avg_score - scores.iloc[3:].mean()
                    if len(group) >= 4:
                        top3_margin = top3_avg_score - scores.iloc[4:].mean()

            results.append({
                "feature_id": feature_id,
                "top1_score": top1_score,
                "top2_avg_score": top2_avg_score,
                "top3_avg_score": top3_avg_score,
                "top1_gap": top1_gap,
                "top1_margin": top1_margin,
                "top2_margin": top2_margin,
                "top3_margin": top3_margin,
                "H_norm": H_norm
            })
        
        df_result = pd.DataFrame(results)
        
        return df_result

    def _compute_normalized_entropy(
        self,
        scores: pd.Series
    ) -> float:
        """"""
        probs = scores / scores.sum()
        entropy = -np.sum(probs * np.log2(probs))
        k = len(scores)
        H_norm = entropy / np.log2(k) if k > 1 else 0
    
        return H_norm
    
    def _process_total_metrics(
        self,
        df_list: List[pd.DataFrame],
        params: dict,
        mode: str = "weights"
    ) -> pd.DataFrame:
        ""
        df_metrics_list = df_list.copy()
        match mode:
            case "weights":
                equalize_weights = params["equalize_weights"]
                dir_suffix = "equal_param_set" if equalize_weights else "unequal_param_set"
            case "spec_sim_type":
                dir_suffix = "spec_sim_type"
            case _:
                raise ValueError(f"Unknown hyperparams: {mode}")

        out_dir = params["base_out_dir"]
        total_metrics_file = os.path.join(out_dir, f"{dir_suffix}_metrics.csv")
        total_metrics_file = str(total_metrics_file).replace("\\", "/")
        
        df_total = pd.concat(df_metrics_list, axis=0, ignore_index=True)
        df_total.to_csv(total_metrics_file, index=None)
        
        print(f"write total metrics dataframe {total_metrics_file} done!")
    
        return df_total
   

class MetricsPlotter:
    def __init__(
        self, 
        out_dir: str = ""
    ):
        self.out_dir = out_dir
        self.metrics_list = [
            "top1_score",
            "top2_avg_score",
            "top3_avg_score",
            "top1_gap",
            "top1_margin",
            "top2_margin",
            "top3_margin",
            "H_norm",
        ]
    
    def metrics_boxplot(
        self, 
        df_metrics: pd.DataFrame,
        out_dir: str = None,
        figsize: Tuple[float, float] = (5, 5),
        mode: str = "weights"
    ) -> None:
        if out_dir is None:
            out_dir = self.out_dir
        
        if not os.path.isdir(out_dir): 
            os.makedirs(out_dir)
        
        match mode:
            case "weights":
                self.weights_params_metrics_boxplot(
                        df_metrics=df_metrics,
                        out_dir=out_dir,
                        figsize=figsize
                    )
            case "spec_sim_type":
                    self.spec_sim_metrics_boxplot(
                        df_metrics=df_metrics,
                        out_dir=out_dir,
                        figsize=figsize
                    )
            case _:
                raise ValueError(f"Unknown type: {mode}")
    
    def weights_params_metrics_boxplot(
        self,
        df_metrics: pd.DataFrame,
        out_dir: str,
        figsize: Tuple[float, float] = (5, 5)
    ) -> None:

        df_tmp = df_metrics.copy()
        removed_columns = [
            "work_dir", "feature_id", "spec_sim_type",
            "rt_weight", "precursor_mz_weight",
        ]
        df_filtered = df_tmp.drop(columns=removed_columns, inplace=False)
        df_long = df_filtered.melt(
            id_vars=["spec_weight"], 
            var_name="score_type", 
            value_name="score"
        )

        for metrics in self.metrics_list:
            df = df_long.copy()
            df_plot = df[df["score_type"] == metrics]
            fig, ax = self._weights_params_metrics_boxplot_single(
                        data=df_plot, 
                        x="spec_weight",
                        y="score",
                        hue="spec_weight",
                        metrics=metrics,
                        figsize=figsize
                    )

            out_file = os.path.join(out_dir, f"equal_params_{metrics}.pdf")
            out_file = out_file.replace("\\", "/")
            fig.savefig(fname=out_file, dpi=400)
            plt.close(fig)
    
    def _weights_params_metrics_boxplot_single(
        self,
        data: pd.DataFrame,
        x: str = None,
        y: str = None,
        hue: str = None,
        metrics: str = "",
        figsize: Tuple[float, float] = (6, 5)
    ) -> Tuple[Figure, Axes]:
        fig, ax = plt.subplots(figsize=figsize)
        ax: Axes = ax
        ax = sns.boxplot(
            x=x,
            y=y, 
            hue=hue, 
            data=data, 
            ax=ax
        )
        ax.set_title(f"Boxplot of {metrics} by Spec Sim Type", fontsize=15)
        ax.set_xlabel("Spectrum weight", fontsize=14)
        ax.set_ylabel("Score",  fontsize=14)
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(axis="x", labelrotation=45, labelright=True, labelsize=10)
        ax.legend(title="Spectrum weight", loc="center left", bbox_to_anchor=(1, 0.5))
    
        fig.tight_layout()
    
        return fig, ax
        
    def spec_sim_metrics_boxplot(
        self,
        df_metrics: pd.DataFrame,
        out_dir: str,
        figsize: Tuple[float, float] = (5, 5)
    ) -> None:
        """"""
        removed_columns = [
            "work_dir", "feature_id", 
            "rt_weight", "precursor_mz_weight", "spec_weight", 
        ]
        df_tmp = df_metrics.copy()
        df_filtered = df_tmp.drop(columns=removed_columns, inplace=False)
        df_long = df_filtered.melt(
            id_vars=["spec_sim_type"], 
            var_name="score_type", 
            value_name="score"
        )
        
        for metrics in self.metrics_list:
            df = df_long.copy()
            df_plot = df[df["score_type"] == metrics]
            fig, ax = self._spec_sim_metrics_boxplot_single(
                        data=df_plot, 
                        x="spec_sim_type",
                        y="score",
                        metrics=metrics,
                        figsize=figsize
                    )

            out_file = os.path.join(out_dir, f"spec_sim_{metrics}.pdf")
            out_file = out_file.replace("\\", "/")
            fig.savefig(fname=out_file, dpi=400)
            plt.close(fig)
    
    def _spec_sim_metrics_boxplot_single(
        self,
        data: pd.DataFrame,
        x: str ="spec_sim_type",
        y: str = "score",
        metrics: str = "",
        figsize: Tuple[float, float] = (5, 5)
    ) -> Tuple[Figure, Axes]:
        fig, ax = plt.subplots(figsize=figsize)
        ax: Axes = ax
        ax = sns.boxplot(
            x=x,
            y=y, 
            hue=x, 
            data=data, 
            ax=ax
        )
        ax.set_title(f"Boxplot of {metrics} by Spec Sim Type", fontsize=15)
        ax.set_xlabel("Spec Sim Type",  fontsize=14)
        ax.set_ylabel("Score",  fontsize=14)
        ax.tick_params(axis="x", labelrotation=45, 
                   labelright=True, labelsize=10)
        fig.tight_layout()
    
        return fig, ax




def run_optimize_weights(
    config: dict,
    weights_params: dict
) -> pd.DataFrame:
    optimizer_equal = HyperparamOptimizer(
                            config=config,
                            params=weights_params
                        )
    df_metrics  = optimizer_equal.optimize(mode="weights")
    
    return df_metrics
 

def run_optimize_spec_sim_type(
    config: dict,
    spec_params: dict
) -> pd.DataFrame:
    optimizer_spec_sim = HyperparamOptimizer(
                            config=config,
                            params=spec_params
                        )
    df_metrics  = optimizer_spec_sim.optimize(mode="spec_sim_type")
    
    return df_metrics


def run_optimize_hela() -> None:
    """"""
    HELA_CONFIG = {
        "pos":{
            "ozid_database_file": "database/clean_data/ozonolysis_delta_mass.json",
            "ozid": {
                "ms1_file": "example/data/Hela/pos-Hela-O3_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/Hela/pos-Hela-O3.mgf",
                "ms2_file_type": "mgf",
            },
            "cid": {
                "ms1_file": "example/data/Hela/pos-Hela-N2_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/Hela/pos-Hela-N2.mgf",
                "ms2_file_type": "mgf",
            }
        },
        "neg": {
            "ozid_database_file": "database/clean_data/ozonolysis_delta_mass.json",
            "ozid": {
                "ms1_file": "example/data/Hela/neg-Hela-O3_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/Hela/neg-Hela-O3.mgf",
                "ms2_file_type": "mgf",
            },
            "cid": {
                "ms1_file": "example/data/Hela/neg-Hela-N2_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/Hela/neg-Hela-N2.mgf",
                "ms2_file_type": "mgf",
            }
        },
    }
    WEIGHTS_PARAMS = {
        "rt_tol": 0.5,
        "rt_tol_mode": "absolute",
        "mz_tol": 1.0,
        "mz_tol_mode": "Da",
        "ms2_spectrum_similarity_type": "None",
        "score_threshold": 0.0,
        "sn_threshold": 2.0,
        "top_n": 10,
        "remove_false_positive": True,
        "base_out_dir": "example/hyperparams_opt/Hela",
        "equalize_weights": True,
        "weight_step": 0.1
    }
    SPEC_SIM_PARAMS = {
        "rt_tol": 0.5,
        "rt_tol_mode": "absolute",
        "mz_tol": 1.0,
        "mz_tol_mode": "Da",
        "ms2_spectrum_similarity_type": ["None", "NIST-LC", "NIST-GC", "SQRT", "MassBank",
                                         "ModifiedCosine", "NeutralLossesCosine", "CosineHungarian"],
        "score_threshold": 0.0,
        "sn_threshold": 2.0,
        "top_n": 10,
        "remove_false_positive": True,
        "base_out_dir": "example/hyperparams_opt/Hela",
        "rt_weight": 0.0,
        "precursor_mz_weight": 0.0,
        "spec_weight": 1.0,
    }
   
    config_file = "example/hyperparams_opt/Hela/file_config.yaml"
    weights_params_file = "example/hyperparams_opt/Hela/weights_params_config.yaml"
    spec_sim_params_file = "example/hyperparams_opt/Hela/spec_sim_params_config.yaml"
    
    dict_to_yaml(
        file=config_file,
        data=HELA_CONFIG
    )
    dict_to_yaml(
        file=weights_params_file,
        data=WEIGHTS_PARAMS
    )
    dict_to_yaml(
        file=spec_sim_params_file,
        data=SPEC_SIM_PARAMS
    )
    config = yaml_to_dict(file=config_file)
    weights_params = yaml_to_dict(file=weights_params_file)
    spec_sim_params = yaml_to_dict(file=spec_sim_params_file)
    
    # get metrics
    optimizer_spec_sim = HyperparamOptimizer(
                            config=config,
                            params=spec_sim_params
                        )
    df_metrics_spec_sim  = optimizer_spec_sim.optimize(mode="spec_sim_type")
    
    optimizer_weights = HyperparamOptimizer(
                            config=config,
                            params=weights_params
                        )
    df_metrics_weights  = optimizer_weights.optimize(mode="weights")

    # plot
    metrics_plotter = MetricsPlotter(out_dir="example/hyperparams_opt/Hela-figures")
    metrics_plotter.metrics_boxplot(
        df_metrics=df_metrics_spec_sim, 
        mode="spec_sim_type",
        out_dir="example/hyperparams_opt/Hela-figures",
        figsize=(5, 5)
    )    
    metrics_plotter.metrics_boxplot(
        df_metrics=df_metrics_weights, 
        mode="weights",
        out_dir="example/hyperparams_opt/Hela-figures",
        figsize=(5, 5)
    )


def run_optimize_nist_plasma() -> None:
    """"""
    NIST_PLASMA_CONFIG = {
        "pos":{
            "ozid_database_file": "database/clean_data/ozonolysis_delta_mass.json",
            "ozid": {
                "ms1_file": "example/data/nist_plasma/pos-nist_plasma-O3_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/pos-nist_plasma-O3.mgf",
                "ms2_file_type": "mgf",
            },
            "cid": {
                "ms1_file": "example/data/nist_plasma/pos-nist_plasma-N2_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/pos-nist_plasma-N2.mgf",
                "ms2_file_type": "mgf",
            }
        },
        "neg": {
            "ozid_database_file": "database/clean_data/ozonolysis_delta_mass.json",
            "ozid": {
                "ms1_file": "example/data/nist_plasma/neg-nist_plasma-O3_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/neg-nist_plasma-O3.mgf",
                "ms2_file_type": "mgf",
            },
            "cid": {
                "ms1_file": "example/data/nist_plasma/neg-nist_plasma-N2_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/neg-nist_plasma-N2.mgf",
                "ms2_file_type": "mgf",
            }
        },
    }
    
    SPEC_SIM_PARAMS = {
        "rt_tol": 0.5,
        "rt_tol_mode": "absolute",
        "mz_tol": 1.0,
        "mz_tol_mode": "Da",
        "ms2_spectrum_similarity_type": ["None", "NIST-LC", "NIST-GC", "SQRT", "MassBank",
                                         "ModifiedCosine", "NeutralLossesCosine", "CosineHungarian"],
        "score_threshold": 0.0,
        "sn_threshold": 2.0,
        "top_n": 10,
        "remove_false_positive": True,
        "base_out_dir": "example/hyperparams_opt/nist_plasma",
        "rt_weight": 0.0,
        "precursor_mz_weight": 0.0,
        "spec_weight": 1.0,
    }
   
    WEIGHTS_PARAMS = {
        "rt_tol": 0.5,
        "rt_tol_mode": "absolute",
        "mz_tol": 1.0,
        "mz_tol_mode": "Da",
        "ms2_spectrum_similarity_type": "None",
        "score_threshold": 0.0,
        "sn_threshold": 2.0,
        "top_n": 10,
        "remove_false_positive": True,
        "base_out_dir": "example/hyperparams_opt/nist_plasma",
        "equalize_weights": True,
        "weight_step": 0.1
    }

    config_dir = "example/hyperparams_opt/nist_plasma"
    if not os.path.isdir(config_dir): 
        os.makedirs(config_dir)
    config_file = os.path.join(config_dir, "file_config.yaml")
    config_file = config_file.replace("\\", "/")
    weights_params_file = os.path.join(config_dir, "weights_params_config.yaml")
    weights_params_file = weights_params_file.replace("\\", "/")
    spec_sim_params_file = os.path.join(config_dir, "spec_sim_params_config.yaml")
    spec_sim_params_file = spec_sim_params_file.replace("\\", "/")
    
    dict_to_yaml(
        file=config_file,
        data=NIST_PLASMA_CONFIG
    )
    dict_to_yaml(
        file=weights_params_file,
        data=WEIGHTS_PARAMS
    )
    dict_to_yaml(
        file=spec_sim_params_file,
        data=SPEC_SIM_PARAMS
    )
    config = yaml_to_dict(file=config_file)
    weights_params = yaml_to_dict(file=weights_params_file)
    spec_sim_params = yaml_to_dict(file=spec_sim_params_file)
    
    # get metrics
    optimizer_spec_sim = HyperparamOptimizer(
                            config=config,
                            params=spec_sim_params
                        )
    df_metrics_spec_sim  = optimizer_spec_sim.optimize(mode="spec_sim_type")
    
    optimizer_weights = HyperparamOptimizer(
                            config=config,
                            params=weights_params
                        )
    df_metrics_weights  = optimizer_weights.optimize(mode="weights")

    # plot
    
    metrics_plotter = MetricsPlotter(out_dir="example/hyperparams_opt/nist-plasma-figures")
    metrics_plotter.metrics_boxplot(
        df_metrics=df_metrics_spec_sim, 
        mode="spec_sim_type",
        out_dir="example/hyperparams_opt/nist-plasma-figures",
        figsize=(5, 5)
    )    
    metrics_plotter.metrics_boxplot(
        df_metrics=df_metrics_weights, 
        mode="weights",
        out_dir="example/hyperparams_opt/nist-plasma-figures",
        figsize=(5, 5)
    )


if __name__ == "__main__":
    
    #run_optimize_hela()
    run_optimize_nist_plasma()
    """
    metrics_weights_file = "example/hyperparams_opt/nist_plasma/equal_param_set_metrics.csv"
    metrics_spec_sim_file = "example/hyperparams_opt/nist_plasma/spec_sim_type_metrics.csv"
    df_metrics_spec_sim = pd.read_csv(
                            metrics_weights_file,
                            keep_default_na=False, 
                            na_values=["", "NULL", "NA"]
                        )
    df_metrics_weights = pd.read_csv(metrics_spec_sim_file)
    
    metrics_plotter = MetricsPlotter(out_dir="example/hyperparams_opt/nist-plasma-figures")
    metrics_plotter.metrics_boxplot(
        df_metrics=df_metrics_spec_sim, 
        mode="spec_sim_type",
        out_dir="example/hyperparams_opt/nist-plasma-figures",
        figsize=(5, 5)
    )    
    metrics_plotter.metrics_boxplot(
        df_metrics=df_metrics_weights, 
        mode="weights",
        out_dir="example/hyperparams_opt/nist-plasma-figures",
        figsize=(5, 5)
    )
    """

