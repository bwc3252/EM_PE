# EM PE

Joint GW/EM parameter estimation.

[Documentation](https://em-pe.readthedocs.io/en/latest/)

## Installation

First, clone this repository:

```bash
$ git clone https://github.com/bwc3252/EM_PE
```

Then, switch to the top-level directory, and install:

```bash
$ cd em_pe
$ python setup.py install
```

`setup.py` may ask for root privileges -- if so, you can either run the above
command as root (to install at the system level), or use the `--user` option to
just install to your account:

```bash
$ python setup.py install --user
```

Now you should be ready to import and use it:

```bash
$ python

>>> import em_pe
```
