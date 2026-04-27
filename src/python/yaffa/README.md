# The `Chain` class
It can be used to implicitly iterate over lists of objects including ROOT objects.
Traditional loops like:

``` python
hists = [h1, h2, h3]
for hist in hists:
    hist.Scale(2)
```
are replaced by
``` python
hists = Chain([h1, h2, h3])
hists.Scale(2) # The loop over each histogram is implicit
```

Member functions are broadcasted, this means that this class allows to use any member function as long as it is defined for the corresponding type or object

Methods can be chained too. T

``` python
hists = Chain([h1, h2, h3])
hists.Scale(1. / hists.Integral('width')) # implicit expansion
```
