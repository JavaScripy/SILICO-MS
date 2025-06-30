import silico_ms

if __name__ == "__main__":
    ozonolysis_delta_mass_file = "database/clean_data/ozonolysis_delta_mass.json"
    
    ## pos
    structure_file = "database/clean_data/reference_ms2_pos.tsv"
    ms1_file = "example/data/Hela/pos-Hela_quant_full.csv"
    ms2_file = "example/data/Hela/pos-Hela.mgf"
    out_file = "example/results/Hela/pos-Hela-ms1.csv"
    silico_ms.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )
    
    ## neg 
    structure_file = "database/clean_data/reference_ms2_neg.tsv"
    ms1_file = "example/data/Hela/neg-Hela_quant_full.csv"
    ms2_file = "example/data/Hela/neg-Hela.mgf"
    out_file = "example/results/Hela/neg-Hela-ms1.csv"
    silico_ms.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )