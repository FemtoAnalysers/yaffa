class Chain:
    def __init__(self, items):
        self.items = list(items)
        self.results = {}

    def __getattr__(self, name):
        def wrapper(*args, **kwargs):
            out = []
            for obj in self.items:
                out.append(getattr(obj, name)(*args, **kwargs))
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
