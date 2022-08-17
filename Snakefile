# fichiers visés : (les graphes du fractionation biais + ) un fichier de résultats des tests stats ???
rule targets:
    input:
        "results/python/fractionation_stat.txt"


rule measure_frac_biais:
    input: 

    output:

    shell:
        "cd scripts/python"
        "python main2.py > ../../results/python/fractionation_stat.txt"


rule cut_tabulation:
    input:
        "results/iadhore/all/multiplicon_pairs.txt"
    output:
        "results/iadhore/all/multiplicon_pairs_modified.txt"
    shell:
        "sed 's:\t\t*:\t:g' {input} > {output}"


rule iadhore:
    input: 
        "submit/iadhore_job.sh"
        