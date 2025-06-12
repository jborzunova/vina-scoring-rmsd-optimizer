import os
import pickle


def log(msg):
    #os.system(f'echo {msg}')
    print(msg)


class History:
    def __init__(self):
        self.hist = {}

    def append(self, params, raw_rmsd, mean_leaky_relu_rmsd):
        self.hist[tuple(params)] = (raw_rmsd, mean_leaky_relu_rmsd)
    
    def all(self):
        return self.hist

    def get_raw_rmsd(self, params):
        return self.hist.get(tuple(params))[0]

    def save(self, file_name):
        with open(f'../{file_name}', 'wb') as file:
            pickle.dump(self.hist, file)


