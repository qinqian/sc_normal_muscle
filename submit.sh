
snakemake --cluster-config bsub.json -s Snakefile.cytotrace -j 8 --cluster 'bsub -q rerunnable -J dive -m "cn031 cn036 cn033 cn035 cn038 cn037 cn034 cn039 cn032 cn034 cn038 cn036 cn032"'

