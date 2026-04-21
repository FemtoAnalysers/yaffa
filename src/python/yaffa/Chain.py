class Chain:
    def __init__(self, chain):
        self._chain = chain
        self.results = {}

    def do(self, func, *args, id=None):
        if id and self.results.get(id):
            raise RuntimeError('results with id="{}" already exists'.format(id))
        
        results = []
        for ring in self._chain:
            results.append(func(ring, *args))

        if id:
            self.results[id] = results

        return self
