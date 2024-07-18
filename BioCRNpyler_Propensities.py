# File with all the propensities that can be imported

from biocrnpyler import *

# Create Propensity Class
# Logistic Positive Propensity 

class LogisticPositive(Propensity):
    def __init__(self, k: float, c: Species, c_max: float):
        """Logistic positive propensity is a nonlinear propensity with the following formula.
            p(c; k, c_max) = k * c * (1 - c/c_max)
            This propensity usually helps only if used for logistic growth function. In the model it is used for modeling the cell population
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        """
        Propensity.__init__(self)
        self.k = k
        self.c = c
        self.c_max = c_max
        self.name = 'logisticpositive'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k {self.c.pretty_print(**kwargs)}* ( 1 - {self.c.pretty_print(**kwargs)}/c_max)'

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        rate_formula = f"{k}*{c} * (1 -  {c}/{c_max})"

        return rate_formula
    
# Create Propensity Class
# decrease_linear_growth Propensity 

class DilutionRate(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, x : Species):
        """decrease_linear_growth is a propensity with the following formula.
            p(c,x; k, c_max) = k * x * (1 - c/c_max)
            Useful in cases where a rate decreases linearly with the specific growth rate of a cell.
            Equal to dilution rate if rate constant is the logistic growth rate, as used in the model 
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param x: species (chemical_reaction_network.species) on which the rate depends on 
        """
        Propensity.__init__(self)
        self.k = k
        self.c = c
        self.c_max = c_max
        self.x = x
        self.name = 'dilution_rate'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c
        
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, new_x):
        self._x = self._check_species(new_x)
        self.propensity_dict['species']['x'] = self.x
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * {self.x.pretty_print(**kwargs)} * (1 - {self.c.pretty_print(**kwargs)}/c_max )'

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        x = propensity_dict['species']['x']
        rate_formula = f"{k} * {x} * (1 - {c}/{c_max})"

        return rate_formula

    

# Create Propensity Class
# hill_growth Propensity 

class Hill_with_growth(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, x : Species, n : float, positive = True):
        """Hill_with_growth is a propensity with the following formula.
            p(c,x; k, c_max) = k * x * (f**n/(1 + f**n))
            or
            p(c,x; k, c_max) = k*x*(1/(1 + f**n) + b), where f = c/c_max
            depending on the condition 'positive' passed (can be True or False respectivley)
            This propensity can be used when there is a switching behaviour from exponential to stationary phase 
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param b : basal fraction of rate if needed (usually 0<b<1, 0 by default)
        """
        Propensity.__init__(self)
        self.k = k
        self.n = n
        self.positive = positive
        self.c = c
        self.c_max = c_max
        self.x = x
        self.name = 'Hill_with_growth'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def n(self):
        if isinstance(self._n, Parameter):
            return self._n.value
        else:
            return self._n

    @n.setter
    def n(self, new_n):
        self._n = self._check_parameter(new_n)
        self.propensity_dict['parameters']['n'] = self._n
    
    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c
        
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, new_x):
        self._x = self._check_species(new_x)
        self.propensity_dict['species']['x'] = self.x
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        if self.positive:
            return f' Kf = k * {self.x.pretty_print(**kwargs)} * (({self.c.pretty_print(**kwargs)}/c_max)^n/(1 + ({self.c.pretty_print(**kwargs)}/c_max)^n))'
        else:
            return f' Kf = k * {self.x.pretty_print(**kwargs)} *(1/(1 + ({self.c.pretty_print(**kwargs)}/c_max)^n))'    
    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        x = propensity_dict['species']['x']

        if self.positive:
            rate_formula = f"{k} * {x} * (({c}/{c_max})^{n}/(1 + ({c}/{c_max})^{n}))"
            
        else:
            rate_formula = f"{k} * {x} * (1/(1 + ({c}/{c_max})^{n}))"
             
        return rate_formula
    
    

    
# Create Propensity Class
# two_species_source Propensity 

class two_species_source(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, S1 : Species, S2 : Species, S_min : float):
        """two_species_source is a propensity with the following formula.
            p(c,x; k, c_max) = k * (S1+S2 - S_min) * (1 - 2 * c/c_max)
            These propensities are designed to model the reaction rates where we consider sum of different species 
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param S1 : species (species bound or unbound)
        :param S2 : species (species bound or unbound)
        :param Smin : minimum possible concentration of the total bound+unbound species
        """
        Propensity.__init__(self)
        self.k = k
        self.S1 = S1
        self.S2 = S2
        self.S_min = S_min
        self.c = c
        self.c_max = c_max
        self.name = 'two_species_source'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c

    @property
    def S1(self):
        return self._S1

    @S1.setter
    def S1(self, new_S1):
        self._S1 = self._check_species(new_S1)
        self.propensity_dict['species']['S1'] = self.S1

    @property
    def S2(self):
        return self._S2

    @S2.setter
    def S2(self, new_S2):
        self._S2 = self._check_species(new_S2)
        self.propensity_dict['species']['S2'] = self.S2

    @property
    def S_min(self):
        if isinstance(self._S_min, Parameter):
            return self._S_min.value
        else:
            return self._S_min

    @S_min.setter
    def S_min(self, new_S_min):
        self._S_min = self._check_parameter(new_S_min)
        self.propensity_dict['parameters']['S_min'] = self._S_min

    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * ({self.S1.pretty_print(**kwargs)} + {self.S2.pretty_print(**kwargs)} - S_min) *(1 - 2*{self.c.pretty_print(**kwargs)}/c_max)'
    
    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for these types of Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw
    
    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        S1 = propensity_dict['species']['S1']
        S2 = propensity_dict['species']['S2']
        S_min = propensity_dict['parameters']['S_min']

        rate_formula = f"{k} * ({S1} + {S2} - {S_min}) * (1 - 2 * {c}/{c_max})"
        
        return rate_formula
    
# Create Propensity Class
# three_species_source Propensity 

class three_species_source(Propensity):
    def __init__(self, k: float, c: Species, c_max: float,S1 : Species, S2 : Species, 
                 S3 : Species, S_min : float):
        """three_species_source is a propensity with the following formula.
            p(c,x; k, c_max) = k * (S1+S2+S3 - S_min) * (1 - 2 * c/c_max)
            
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param S1 : species (species bound or unbound)
        :param S2 : species (species bound or unbound)
        :param S3 : species (species bound or unbound)
        :param Smin : minimum possible concentration of the total bound+unbound species
        """
        Propensity.__init__(self)
        self.k = k
        self.S1 = S1
        self.S2 = S2
        self.S3 = S3
        self.S_min = S_min
        self.c = c
        self.c_max = c_max
        self.name = 'three_species_source'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c

    @property
    def S1(self):
        return self._S1

    @S1.setter
    def S1(self, new_S1):
        self._S1 = self._check_species(new_S1)
        self.propensity_dict['species']['S1'] = self.S1

    @property
    def S2(self):
        return self._S2

    @S2.setter
    def S2(self, new_S2):
        self._S2 = self._check_species(new_S2)
        self.propensity_dict['species']['S2'] = self.S2

    @property
    def S3(self):
        return self._S3

    @S3.setter
    def S3(self, new_S3):
        self._S3 = self._check_species(new_S3)
        self.propensity_dict['species']['S3'] = self.S3

    @property
    def S_min(self):
        if isinstance(self._S_min, Parameter):
            return self._S_min.value
        else:
            return self._S_min

    @S_min.setter
    def S_min(self, new_S_min):
        self._S_min = self._check_parameter(new_S_min)
        self.propensity_dict['parameters']['S_min'] = self._S_min

    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * ({self.S1.pretty_print(**kwargs)} + {self.S2.pretty_print(**kwargs)}  + {self.S3.pretty_print(**kwargs)} - S_min) * (1 - 2 * {self.c.pretty_print(**kwargs)}/c_max)'
    
    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for these types of Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw
    
    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        S1 = propensity_dict['species']['S1']
        S2 = propensity_dict['species']['S2']
        S3 = propensity_dict['species']['S3']
        S_min = propensity_dict['parameters']['S_min']

        rate_formula = f"{k} * ({S1} + {S2} + {S3} - {S_min}) * (1 - 2 * {c}/{c_max})"
        
        return rate_formula
    
# Create Propensity Class
# five_species_source Propensity 

class five_species_source(Propensity):
    def __init__(self, k: float, c: Species, c_max: float,S1 : Species, S2 : Species, 
                 S3 : Species, S4 : Species, S5 : Species, S_min : float):
        """five_species_source is a propensity with the following formula.
            p(c,x; k, c_max) = k * (S1+S2+S3 - S_min) * (1 - 2 * c/c_max)
            
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param S1 : species (species bound or unbound)
        :param S2 : species (species bound or unbound)
        :param S3 : species (species bound or unbound)
        :param S4 : species (species bound or unbound)
        :param S5 : species (species bound or unbound)
        :param Smin : minimum possible concentration of the total bound+unbound species
        """
        Propensity.__init__(self)
        self.k = k
        self.S1 = S1
        self.S2 = S2
        self.S3 = S3
        self.S4 = S4
        self.S5 = S5
        self.S_min = S_min
        self.c = c
        self.c_max = c_max
        self.name = 'five_species_source'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c

    @property
    def S1(self):
        return self._S1

    @S1.setter
    def S1(self, new_S1):
        self._S1 = self._check_species(new_S1)
        self.propensity_dict['species']['S1'] = self.S1

    @property
    def S2(self):
        return self._S2

    @S2.setter
    def S2(self, new_S2):
        self._S2 = self._check_species(new_S2)
        self.propensity_dict['species']['S2'] = self.S2

    @property
    def S3(self):
        return self._S3

    @S3.setter
    def S3(self, new_S3):
        self._S3 = self._check_species(new_S3)
        self.propensity_dict['species']['S3'] = self.S3

    @property
    def S4(self):
        return self._S4

    @S4.setter
    def S4(self, new_S4):
        self._S4 = self._check_species(new_S4)
        self.propensity_dict['species']['S4'] = self.S4

    @property
    def S5(self):
        return self._S5

    @S5.setter
    def S5(self, new_S5):
        self._S5 = self._check_species(new_S5)
        self.propensity_dict['species']['S5'] = self.S5

    @property
    def S_min(self):
        if isinstance(self._S_min, Parameter):
            return self._S_min.value
        else:
            return self._S_min

    @S_min.setter
    def S_min(self, new_S_min):
        self._S_min = self._check_parameter(new_S_min)
        self.propensity_dict['parameters']['S_min'] = self._S_min

    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * ({self.S1.pretty_print(**kwargs)} + {self.S2.pretty_print(**kwargs)}  + {self.S3.pretty_print(**kwargs)} + {self.S4.pretty_print(**kwargs)} + {self.S5.pretty_print(**kwargs)} - S_min) * (1 - 2 * {self.c.pretty_print(**kwargs)}/c_max)'
    
    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for these types of Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw
    
    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        S1 = propensity_dict['species']['S1']
        S2 = propensity_dict['species']['S2']
        S3 = propensity_dict['species']['S3']
        S4 = propensity_dict['species']['S4']
        S5 = propensity_dict['species']['S5']
        S_min = propensity_dict['parameters']['S_min']

        rate_formula = f"{k} * ({S1} + {S2} + {S3} + {S4} + {S5} - {S_min}) * (1 - 2 * {c}/{c_max})"
        
        return rate_formula
    

    

# Create Propensity Class
# decrease_linear_growth_basal Propensity 

class DilutionRate_basal(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, x : Species, b : float):
        """decrease_linear_growth is a propensity with the following formula.
            p(c,x; k, c_max) = k * x * (1 - c/c_max + b)
            This propensity was designed to model reactions that decrease linearly with increasing cell population but do not go to zero so added a basal rate
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param b : basal fraction of rate if needed (usually 0<b<1, 0 by default)
        """
        Propensity.__init__(self)
        self.k = k
        self.b = b
        self.c = c
        self.c_max = c_max
        self.x = x
        self.name = 'DilutionRate_basal'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def b(self):
        if isinstance(self._b, Parameter):
            return self._b.value
        else:
            return self._b

    @b.setter
    def b(self, new_b):
        self._b = self._check_parameter(new_b)
        self.propensity_dict['parameters']['b'] = self._b

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c
        
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, new_x):
        self._x = self._check_species(new_x)
        self.propensity_dict['species']['x'] = self.x
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * {self.x.pretty_print(**kwargs)} * (1 - {self.c.pretty_print(**kwargs)}/c_max + b)'

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        b = propensity_dict['parameters']['b']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        x = propensity_dict['species']['x']
        rate_formula = f"{k} * {x} * (1 - {c}/{c_max} + {b})"
            
        return rate_formula
    


# Create Propensity Class
# Hill_with_growth_basal Propensity 

class Hill_with_growth_basal(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, x : Species, b : float, n : float, positive = True):
        """hill_growth is a propensity with the following formula.
            p(c,x; k, c_max) = k * x * (f**n/(1 + f**n) + b)
            or
            p(c,x; k, c_max) = k*x*(1/(1 + f**n) + b), where f = c/c_max
            depending on the condition passed
            As name suggests, used for modeling reactinos that follow Hill type propensity but have a basal rate too 
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        :param b : basal fraction of rate if needed (usually 0<b<1, 0 by default)
        """
        Propensity.__init__(self)
        self.k = k
        self.b = b
        self.n = n
        self.positive = positive
        self.c = c
        self.c_max = c_max
        self.x = x
        self.name = 'Hill_with_growth_basal'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def b(self):
        if isinstance(self._b, Parameter):
            return self._b.value
        else:
            return self._b

    @b.setter
    def b(self, new_b):
        self._b = self._check_parameter(new_b)
        self.propensity_dict['parameters']['b'] = self._b

    @property
    def n(self):
        if isinstance(self._n, Parameter):
            return self._n.value
        else:
            return self._n

    @n.setter
    def n(self, new_n):
        self._n = self._check_parameter(new_n)
        self.propensity_dict['parameters']['n'] = self._n
    
    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c
        
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, new_x):
        self._x = self._check_species(new_x)
        self.propensity_dict['species']['x'] = self.x
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        if self.positive:
            return f' Kf = k * {self.x.pretty_print(**kwargs)} *(({self.c.pretty_print(**kwargs)}/c_max)^n/(1 + ({self.c.pretty_print(**kwargs)}/c_max)^n) + b)'
        else:
            return f' Kf = k * {self.x.pretty_print(**kwargs)} *(1/(1 + ({self.c.pretty_print(**kwargs)}/c_max)^n) + b)'    
    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        b = propensity_dict['parameters']['b']
        n = propensity_dict['parameters']['n']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        x = propensity_dict['species']['x']


        if self.positive:
            rate_formula = f"{k} * {x} * (({c}/{c_max})^{n}/(1 + ({c}/{c_max})^{n}) + {b})"
            
        else:
            rate_formula = f"{k} * {x} * (1/(1 + ({c}/{c_max})^{n}) + {b})"
             
        return rate_formula
    
class gamma(Propensity):
    def __init__(self, k: float, c: Species, c_max: float):
        """gamma is a propensity with the following formula.
            p(c,x; k, c_max) = k * (1 - c/c_max) * (c/c_max)
            This propensity is used to model amino acid replenishment that increases with increase in growth but also decreases as cell reaches stationary phase
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        """
        Propensity.__init__(self)
        self.k = k
        self.c = c
        self.c_max = c_max
        self.name = 'gamma'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c
        
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * (1 - {self.c.pretty_print(**kwargs)}/c_max ) * ({self.c.pretty_print(**kwargs)}/c_max)'

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        rate_formula = f"{k} * (1 - {c}/{c_max}) * ({c}/{c_max})"

        return rate_formula


# Protein folding propensity 
class protein_folding_propensity(Propensity):
    def __init__(self, k: float, x : Species, P_max : float, b : float, c : Species, c_max : float, p : Species):
        """protein_folding_propensity is a propensity with the following formula.
            p(c,x; k, c_max) = k * P * ((c/c_max) * (1 - c/c_max) * (1 - (P_total/Pmax) + b)
        :param k: rate constant (float)
        :param P_max: maximum protein capacity (float)
        :param b : basal fraction of rate if needed (usually 0<b<1, 0 by default)
        """
        Propensity.__init__(self)
        self.k = k
        self.P_max = P_max
        self.c = c
        self.c_max = c_max
        self.x = x
        self.b = b
        self.p = p
        self.name = 'protein_folding_propensity'

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def b(self):
        if isinstance(self._b, Parameter):
            return self._b.value
        else:
            return self._b

    @b.setter
    def b(self, new_b):
        self._b = self._check_parameter(new_b)
        self.propensity_dict['parameters']['b'] = self._b

    @property
    def P_max(self):
        if isinstance(self._P_max, Parameter):
            return self._P_max.value
        else:
            return self._P_max

    @P_max.setter
    def P_max(self, new_P_max):
        self._P_max = self._check_parameter(new_P_max)
        self.propensity_dict['parameters']['P_max'] = self._P_max

    @property
    def c_max(self):
        if isinstance(self._c_max, Parameter):
            return self._c_max.value
        else:
            return self._c_max

    @c_max.setter
    def c_max(self, new_c_max):
        self._c_max = self._check_parameter(new_c_max)
        self.propensity_dict['parameters']['c_max'] = self._c_max

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new_c):
        self._c = self._check_species(new_c)
        self.propensity_dict['species']['c'] = self.c


    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, new_p):
        self._p = self._check_species(new_p)
        self.propensity_dict['species']['p'] = self.p
    
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, new_x):
        self._x = self._check_species(new_x)
        self.propensity_dict['species']['x'] = self.x
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * {self.x.pretty_print(**kwargs)} *(({self.c.pretty_print(**kwargs)}/C_max) * (1 - {self.c.pretty_print(**kwargs)}/C_max) * (1 - ({self.x.pretty_print(**kwargs)}+{self.p.pretty_print(**kwargs)})/P_max) + b)'

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Creates KineticLaw object for SBML."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError('reverse reactions cannot exist for Logistic type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError("Could not write the rate law for reaction to SBML.\
                             Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        b = propensity_dict['parameters']['b']
        P_max = propensity_dict['parameters']['P_max']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        x = propensity_dict['species']['x']
        p = propensity_dict['species']['p']
        rate_formula = f"{k} * {x} * (({c}/{c_max}) * (1 - ({c}/{c_max})) * (1 - ({x} + {p})/{P_max}) + {b})"
       

        return rate_formula