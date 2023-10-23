## Environment Setup

```
mamba env create -f conda-env.yml -n mod-site
conda activate mod-site
```

## Update Environment

```
mamba env update -n mod-site --file conda-env.yml

```

<!-- compound arguments json explained -->
## arguments description
### compound arguments, description, type and default value
- `mz_tolerance`: mass tolerance for matching | float | 0.01