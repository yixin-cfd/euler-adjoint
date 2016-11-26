import yaml

def defaults():
    return { "mach" : 0.3,
             "ilhs" : 1,
             "resid" : 10,
             "steps" : 100,
             "cfl" : 10.0,
             "aoa" : 0.0,
             "bcs" : [
                 { "face": "jmin", "type": "periodic", "j": [0,    0], "k": [0,  -1] },
                 { "face": "jmax", "type": "periodic", "j": [-1,  -1], "k": [0,  -1] },
                 { "face": "kmin", "type": "wall"    , "j": [0,   -1], "k": [0,   0] },
                 { "face": "kmax", "type": "farfield", "j": [0,   -1], "k": [-1, -1] }
             ]
             }

def read_inputs(filename):
    inputs = defaults()
    with open(filename, "r") as f:
        inputs.update(yaml.load(f))
    return inputs

