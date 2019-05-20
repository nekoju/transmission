---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.0'
      jupytext_version: 0.8.6
  kernelspec:
    display_name: transmission
    language: python
    name: transmission3
---

# Demonstration of parameter estimation using Approximate Bayesian Computation






```python
from __future__ import print_function, division
import transmission as txmn
import pandas as pd

%matplotlib inline
```


 First, simulate a data set. The `sim()` function in Transmission is the
 workhorse function that simulates a geneaology given a set of parameters.
 Generally, it is called by `ms_simulate()` and there is no reason to call
 it directly.
