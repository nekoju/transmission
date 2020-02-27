---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.0'
      jupytext_version: 0.8.6
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
from sympy import *

init_printing()
```

```python
fn = Function('fn')
fm = Function('fm')
a, b, t, r0, C0, C1 = symbols('a b t r0 C0 C1')
Fm = a + (1 - a) * fm(t - 1) - fm(t)
Fn = b + (1 - b) * fn(t - 1) - fn(t)
mtClosed=rsolve(Fm, fm(t))
nucClosed=rsolve(Fn, fn(t)).subs(C0, C1)
```

```python

```
