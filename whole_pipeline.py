import utils as utils



if __name__ == "__main__":
    ozonolysis_delta_mass_file = "database/clean_data/ozonolysis_delta_mass.json"
    
    # 1. Hela 
    ## 1.1 pos
    structure_file = "database/clean_data/reference_ms2_pos.tsv"
    ms1_file = "data/Hela/pos-Hela_quant_full.csv"
    ms2_file = "data/Hela/pos-Hela.mgf"
    out_file = "results/Hela/pos-Hela-ms1.csv"
    utils.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )
    
    ## 1.2 neg 
    structure_file = "database/clean_data/reference_ms2_neg.tsv"
    ms1_file = "data/Hela/neg-Hela_quant_full.csv"
    ms2_file = "data/Hela/neg-Hela.mgf"
    out_file = "results/Hela/neg-Hela-ms1.csv"
    utils.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )
    
    # 2. Liver
    ## 2.1 pos
    structure_file = "database/clean_data/reference_ms2_pos.tsv"
    ms1_file = "data/Liver/pos-liver_quant_full.csv"
    ms2_file = "data/Liver/pos-liver.mgf"
    out_file = "results/Liver/pos-Liver-ms1.csv"
    utils.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )

    ## 2.2 neg
    structure_file = "database/clean_data/reference_ms2_neg.tsv"
    ms1_file = "data/Liver/neg-liver_quant_full.csv"
    ms2_file = "data/Liver/neg-liver.mgf"
    out_file = "results/Liver/neg-Liver-ms1.csv"
    utils.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )




