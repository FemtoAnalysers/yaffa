'''
Class to enable looping over lists of objects. See README.md for documentation
'''

class Chain:
    def __init__(self, items):
        self.items = list(items)
        self.results = {}

    def __getattr__(self, name):
        def wrapper(*args, **kwargs):
            out = []
            for i, obj in enumerate(self.items):
                new_args = [arg.items[i] if isinstance(arg, Chain) else arg for arg in args]
                new_kwargs = {
                    k: (v.items[i] if isinstance(v, Chain) else v)
                    for k, v in kwargs.items()
                }
                out.append(getattr(obj, name)(*new_args, **new_kwargs))
            return Chain(out)
        return wrapper

    def apply(self, func, *args, id=None, **kwargs):
        out = [func(obj, *args, **kwargs) for obj in self.items]

        if id:
            self.results[id] = out

        return Chain(out)

    def collect(self, func, id):
        self.results[id] = [func(obj) for obj in self.items]
        return self

    def __iter__(self):
        return iter(self.items)

    def __getitem__(self, i):
        return self.items[i]

    def __truediv__(self, other):
        if isinstance(other, Chain):
            return Chain([a / b for a, b in zip(self.items, other.items)])
        else:
            return Chain([a / other for a in self.items])

    def __rtruediv__(self, other):
        return Chain([other / a for a in self.items])
