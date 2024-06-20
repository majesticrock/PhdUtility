from itertools import product

def iterate_containers(*containers):
    # Adjust the range based on whether the container is a single value or an array
    index_combinations = product(*[range(1) if not hasattr(container, '__len__') else range(len(container)) for container in containers])
    
    for indices in index_combinations:
        values = [container if not hasattr(container, '__len__') else container[index] for container, index in zip(containers, indices)]
        yield values
        
def total_size(*containers):
    size = 1
    for container in containers:
        if hasattr(container, '__len__'):
            size *= len(container)
    return size

def naming_scheme(Ts, Us, Vs):
    for T, U, V in iterate_containers(Ts, Us, Vs):
        yield f"T={T}/U={U}/V={V}"
        
def naming_scheme_tuples(tuples):
    for T, U, V in tuples:
        yield f"T={T}/U={U}/V={V}"
## Example usage
#import numpy as np
#
#Ts = 1
#Us = 4
#Vs = 8
#
#for T, U, V in iterate_containers(Ts, Us, Vs):
#    # Do something with T, U, V
#    print(f"T: {T}, U: {U}, V: {V}")
#
#for name in naming_scheme(Ts, Us, Vs):
#    print(name)
#
#print(total_size(Ts, Us, Vs))