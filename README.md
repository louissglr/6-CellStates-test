## Run the pipeline from the workflow

```bash
nohup snakemake \
  --configfile ../config/config-glicid-scratch.yml \
  --profile profiles/default/ \
  --rerun-triggers mtime \
  >> snake.out &
```

## Generate the sample summary report

```bash
nohup snakemake \
  --configfile ../config/config-glicid-scratch.yml \
  --profile profiles/default/ \
  --rerun-triggers mtime \
  --report /scratch/nautilus/users/sanglier-l@univ-nantes.fr/6-CellStates-test/results/report.zip \
  >> snake.out &
```
