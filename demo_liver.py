import silico_ms


if __name__ == "__main__":
    ozonolysis_delta_mass_file = "database/clean_data/ozonolysis_delta_mass.json"
    ## pos
    structure_file = "database/clean_data/reference_ms2_pos.tsv"
    ms1_file = "example/data/Liver/pos-liver_quant_full.csv"
    ms2_file = "example/data/Liver/pos-liver.mgf"
    out_file = "example/results/Liver/pos-Liver-ms1.csv"
    silico_ms.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )

    ## neg
    structure_file = "database/clean_data/reference_ms2_neg.tsv"
    ms1_file = "example/data/Liver/neg-liver_quant_full.csv"
    ms2_file = "example/data/Liver/neg-liver.mgf"
    out_file = "example/results/Liver/neg-Liver-ms1.csv"
    silico_ms.algorithm_whole_pipeline(
        ms1_file=ms1_file,
        ms2_file=ms2_file,
        structure_file=structure_file,
        delta_mass_file=ozonolysis_delta_mass_file,
        out_file=out_file
    )