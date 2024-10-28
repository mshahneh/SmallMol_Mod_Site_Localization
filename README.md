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


# documentation
to make the documentation:
`pip install sphinx-rtd-theme`


run the server for local developement of the readme
`python -m http.server 3000`

make the docs:
`make html`

auto generate docs:
`sphinx-apidoc -o source/utilities_doc ../modifinder/utilities/ -e --tocfile index`

# testing
`python -m unittest discover -s tests`