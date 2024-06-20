# Basis usage of the classes

## Preface

If there is no namespace specified in this file, we imply the namespace ```SymbolicOperators```.
Each relevent class defines and overload of ```std::ostream& operator<<```, i.e., they can easily be used as ```std::cout << class_obj```.
The output is formated to be used within an ```align```-environment within LaTeX.
Your code only needs to include ```Wick.hpp```.

## The ```Momentum``` class
This class represents momenta. It includes addition and substraction operators as well as a ```bool add_Q```. _Q_ is defined as _(pi, pi, ...)_, i.e., _n * Q = 0_ for all even _n_.
Besides the normal operators in which you specify the class members, you can also pass a string to the constructor like ```"3k+l-p"``` to create that specific momentum. If you want to add Q here, you can do so by passing ```true``` to the same constructor as a second argument.

## The ```Operator``` class
This class represents the standard fermionic creation and annihilation operators. You can specify its momentum, its indizes and whether it is supposed to be daggered (a creation operator) or not (an annihilation operator).

## The ```Term``` class
This class has various kind of constructors that allow setting coefficient(s), sums, operators and deltas.

### Sums
Sums are contained within the ```SumContainer``` class. It hosts both sums of momenta and sums of spins, each one is accessible via the appropriate class member and its ```operator[]```, e.g., ```container.momenta[i]```.

### Deltas
Deltas are represented via the ```KroneckerDelta``` class template.
The template argument specifies what kind of Delta is to be used, e.g., ```Momentum```.
Similar to how ```std::pair<T1,T2>``` implements ```std::make_pair(T1, T2)``` a function ```make_delta(T,T)``` is provided.

## Creating a Hamiltonian

A Hamiltonian is characterized as ```std::vector<Term>```.
It can consist of any number of individual terms.

## Creating a basis
A basis in this sense refers to the set of operators that need to be commuted.
Each individual basis term can consist of multiple different summands.
Therefore, it is to be represented by an ```std::vector<Term>```. The straightforward idea would be to represent the entire basis as ```std::vector<std::vector<Term>>```.
After you have created a basis, you can obtain the Hermitian conjugate of your basis by means of this code

```
std::vector<term_vec> basis_daggered(basis);
for (auto& t : basis_daggered) {
    hermitianConjugate(t);
}
```

However, often it is advisible to rename the momenta of the conjugate basis because otherwise expressions like [b_k, b_k] will actually assume both k to be the same.
This can be achieved via

```
for (auto& t : basis_daggered) {
	hermitianConjugate(t);
	renameMomenta(t, 'k', 'l');
}
```

# How to commute

First, create an ```std::vector<Term>``` that will hold your result.
Then, you can simply call ```commutator(result, left, right)``` to computed [left, right]. Afterwards, it is recommended to call ```cleanUp(result)``` which will remove unnecessary terms, e.g., those that will be 0 anyways.
Additionally, this function groups identical terms.

A double commutator is equally simple:
Create a buffer for the inner commutator, then perform the inner commutation (and clean up). Afterwards simply perform the outer commutation.


# How to use Wick's theorem

Prerequisite: The terms you want to apply Wick's theorem on are saved in an ```std::vector<Term>```.

## Creating Wick templates

Applying Wick's theorem often involves omiting certain expecation values because you know them to be 0 for symmtry reasons. Therefore the class ```WickOperatorTemplate``` exists. Here, you specify, which kind of expectation values will be finite.
In the following, the meaning of the different attributes is listed:

### ```std::vector<IndexComparison> indexComparison```
Holds various ```IndexComparision```

#### ```struct IndexComparison``` 
If ```any_identical``` is ```true```, any two identical indizes are considered valid. An example would be in the number operator *c_(k, sigma)^+ c_(k, sigma')*: No matter what sigma is, as long as sigma = sigma' the expectation value will be finite.
If ```any_identical``` is ```false```, the members ```base``` and ```other``` become relevant: They define what the indizes need to be, e.g., for a pair annihlation operator *c_(-k down) c_(k up)* one would set ```base``` to _down_ and ```other``` to _up_.

Note, once one operator is set as a template, it is not necessary to set its Hermitian conjugate.

### ```Momentum momentum_difference```
Defines the allowed difference in momentum, e.g., for a number operator, this would be 0.
Note, this also applies to a standard pair creation/annihilation operator, because in total, these operators create/annihilate a particle with *-k* and one with *k*, resulting in 0 net momentum.

### ```OperatorType type```
Specifies what kind of WickOperator will be the result, see ```enum OperatorType``` in ```WickOperator.hpp```.

### ```bool is_sc_type```
Specifies whether the operator is a pair creation/annihilation operator or a standard *c^+ c* type term.


## Apply Wick's theorem
Create an instance of ```WickTermCollector```. Then simply call 

```
wicks_theorem(terms, templates, wicks);
cleanWicks(wicks);
```

## Applying symmetries to the result
There may be some symmetries that simplify your results, e.g. <O^+> = <O>.
These symmetries can be implemented by inheriting from the  ```WickSymmetry``` class and defining the member function ```virtual void apply_to(WickTerm& term) const```.
Then create a ```std::vector<std::unique_ptr<WickSymmetry>> symmetries``` and make use of polymorphism by calling ```cleanWicks(wicks, symmetries)```
There are the following predefined symmetry operations:

#### ```SpinSymmetry```
Changes all spins of the operators in ```term``` to _up_.

#### ```TranslationalSymmetry```
Flips the momenta in such a way, that the first momentum in a term is always positive, i.e., _-k+l_ is changed to _k-l_ while _k-l_ would stay unmodified.

#### ```PhaseSymmetry```
Takes a list of ```OperatorType``` as template arguments.
Removes any dagger from all operators with a type from the list.
Example:
```PhaseSymmetry<SC_Type, CDW_Type>``` removes the dagger from ```SC_Type``` and ```CDW_Type``` operators.

## Output and use the result

You can print the the result to the console or utilize boost's serialization to load it later (or within another program).
Serialization can be achieved via this code

```
std::ofstream ofs("path/to/file.txt");
boost::archive::text_oarchive oa(ofs);
oa << wicks;
ofs.close();
```

To later on load the output use

```
std::ifstream ifs("path/to/file.txt");
boost::archive::text_iarchive ia(ifs);
target.clear();
ia >> target;
ifs.close();
```

or if you want to use this code of the iEoM, there is the class ```TermLoader``` for easy use.
It loads the terms for the matrices M and N and saves same as class members.