---
subtitle: A tool for inferring symbiont transmission mode from metagenomic data
title: transmission
---

# Generating simulated summary statistics

This is a big job and I find it helpful to farm this job out to a computing
cluster. The main program is written for interactive use in a notebook or shell,
so to that end, the cli tool `transmission-priorgen` is included with
transmission. Use `transmission-priorgen --help` for details on available
options. 

## Natively

For example:

```
transmission-priorgen -n 10 -d 5 -M 2 -s 100 \
    -p '{"sigma": (0, 0.1), "tau":(1, 1), "rho"(10, 10)}' \
    --h_opts '{"bias": False}' outfile.pickle
```

## Using the Docker image

```
docker run --rm -v </path/to/host/directory>:/home/jovyan/work \
    mpjuers/transmission:<version> \
    transmission-priorgen [options] /home/jovyan/work/<outfile.pickle>
```

`--rm` removes containers you are finished with while `-v` ("volume") binds
a directory on the host machine to one on the container. <version> should be
prefixed by 'v', e.g. v0.0.3

If you are working on a compute cluster, you might have access to Singularity
rather than Docker. The transition is straightforward:

```
singularity exec --contain -B </path/to/host/directory>:/home/jovyan/work \
    docker://mpjuers/transmission:<version> \
    transmission-priorgen [options] /home/jovyan/work/<outfile.pickle>
```

You may or man not need the `--contain` flag; I needed because I had some
locally installed python modules that conflicted with those in the container.
