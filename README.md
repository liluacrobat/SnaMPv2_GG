# Snakemake Microbiome Pipeline

## Download pipeline
```bash
git clone https://github.com/liluacrobat/SnaMPv2_GG.git
```

## Download reference database

download 'gg' folder from google drive:
https://drive.google.com/open?id=1SBTR3SdiVMyoTIsL8PxVcVqvw2IGzQ_T

put folder 'gg' under the 'database' directory

## Create environment

```bash
module python/py38-anaconda-2020.11
conda create -c bioconda -m -p pyenvs/py35-snakemake python=3.5 pandas snakemake
```

## Load environment

```bash
module python/py38-anaconda-2020.11
eval "$(/util/common/python/py38/anaconda-2020.11/bin/conda shell.bash hook)"
source activate pyenvs/py35-snakemake
```

## Remove environment

```
source deactivate
```

## Workflow

  * Closed reference OTU picking
  * BLAST agaist HOMD database

![workflow_0](./misc/dag.png)

### Usage

1. Prepare sequencing data

  Repalce the `$SOURCE_FILES` with the __zipped__ sequencing result, e.g. `WHI_Repo/RT530_Batch2/*.gz`

  ```bash
  cd input
  ln -s $SOURCE_FILES .
  cd ..
  ```

2. Launch jobs

  The pipeline will utilize CCR resource to parallel execution.
  OTU table and statisics about merge rate, filter rate, hit rate wiil be placed under _table_

  ```bash
  snakemake -p -j 100 --cluster-config cluster.json --cluster "sbatch --partition {cluster.partition} --time {cluster.time} --nodes {cluster.nodes} --ntasks-per-node {cluster.ntasks-per-node}" --latency-wait 60
  ```

3. Results

 Three result files are placed under `table` directory.
 
 * QC_table.txt  
 * raw_OTU_table_collapsed.txt  
 * raw_OTU_table_uncollapsed.txt
 

## Misc

  To remove generated files:
  ```
  snakemake clean
  ```

  Check snakemake version:
  ```
  snakemake -v
  # 3.13.3
  ```

## Known Issues

* cluster time limit: if any of the jobs run of the time limit, it will be silently killed by CCR and the pipeline will keep waiting. To solve this, you need to terminate `Snakemake`, increase time limit specified in `cluster.json` accordingly and rerun the pipeline with extra option `--rerun-incomplete`.
