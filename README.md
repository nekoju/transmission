# Transmission: A tool for inferring symbiont transmission mode from metagenomic data

## Installation

Install r package abc in the usual way, then install Transmission with

`pip install transmission-popgen`


## Generating simulated summary statistics

This is a big job and I find it helpful to farm this job out to a computing
cluster. The main program is written for interactive use in a notebook or shell,
so to that end, the cli tool `transmission-priorgen` is included with
transmission. Use `transmission-priorgen --help` for details on available
options. 

### Natively

For example:

```
transmission-priorgen -n 10 -d 5 -M 2 -s 100 \
    -p '{"eta": (0, 0.1), "tau":(1, 1), "rho"(10, 10)}' \
    --h_opts '{"bias": False}' outfile.pickle
```

### Using the Docker image

There is a docker image available to use the command line tools as well as
to launch a Jupyter notebook server. It can be found [here](https://cloud.docker.com/repository/docker/mpjuers/transmission).
The most recent push to master is tagged `latest`, however it is probably
desirable to use explicit version numbers for reproducibility.

#### To generate priors

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

You may or may not need the `--contain` flag; I needed because I had some
locally installed python modules that conflicted with those in the container.

#### Using the Docker image to run a Jupyter notebook server.

I have had some difficulty compiling rpy2 on my platform (MacOS Mojave).
This may or may not still be the case, but if you wish to avoid any problems
associated with setting up transmission on your machine, it is possible to
run a Jupyter server from the Docker image. First, launch Docker on your
machine and then run, e.g.,

```
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes \
    -v </path/to/local/dir>a:/home/jovyan/work mpjuers/transmission:<version>
```

Omit `-e JUPYTER_ENEABLE_LAB=yes` if you are going to use a vanilla notebook. 
Then you can just open `<hostname>:8888/` in a browser and enter the token that
appears in the terminal and anything saved to `~/work` 
on the container will appear in `</path/to/local/dir>` (and vice-versa).
For more information on running Jupyter Docker containers, go to [https://github.com/jupyter/docker-stacks].

