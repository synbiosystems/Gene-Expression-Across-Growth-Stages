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
# alpha propensity 

class alpha_propensity(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, x : Species):
        """apply_alpha_propensity is a propensity with the following formula.
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
        self.name = 'alpha_propensity'

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
# delta Propensity 

class delta_propensity(Propensity):
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
        self.name = 'delta_propensity'

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
# delta_with_basal Propensity # Special case where there is a basal rate as well

class delta_with_basal(Propensity):
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
    

# Create Propensity Class
# gamma propensity applied for 0 order reactions
    

class gamma(Propensity):
    def __init__(self, k: float, c: Species, c_max: float, n : float):
        """gamma is a propensity with the following formula.
            p(c,x; k, c_max) = k * (1 - c/c_max) * (c/c_max)^n
            This propensity is used to model amino acid replenishment that increases with increase in growth but also decreases as cell reaches stationary phase
        :param k: rate constant (float)
        :param c: species (chemical_reaction_network.species)
        :param c_max: maximum capacity (float)
        """
        Propensity.__init__(self)
        self.k = k
        self.c = c
        self.c_max = c_max
        self.n = n
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
        
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f' Kf = k * ((1 - {self.c.pretty_print(**kwargs)}/c_max ) * ({self.c.pretty_print(**kwargs)}/c_max))^n'

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
        rate_formula = f"{k} * ((1 - {c}/{c_max}) * ({c}/{c_max}))^{n}"

        return rate_formula


# Create Propensity Class
# gamma propensity applied to 1st order reactions with a basal rate
    
class apply_gamma_propensity(Propensity):
    def __init__(self, k: float, x : Species, b : float, n : float, c : Species, c_max : float):
        """apply_gamma_propensity is a propensity with the following formula.
            p(c,x; k, c_max) = k * x * (y^n + b)
        :param k : rate constant (float)
        :param y^n: gamma to the power n (float) where y = f * (1 - f)
        :param b : basal fraction of rate if needed (usually 0<b<1, 0 by default)
        """
        Propensity.__init__(self)
        self.k = k
        self.c = c
        self.c_max = c_max
        self.x = x
        self.b = b
        self.n = n
        self.name = 'apply_gamma_propensity'

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
        return f' Kf = k * {self.x.pretty_print(**kwargs)} *(({self.c.pretty_print(**kwargs)}/C_max) * (1 - {self.c.pretty_print(**kwargs)}/C_max)^n  + b)'

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
        n = propensity_dict['parameters']['n']
        rate_formula = f"{k} * {x} * (({c}/{c_max}) * (1 - ({c}/{c_max}))^{n}  + {b})"
       

        return rate_formula
    
# Resource allocation propensity 
class resource_allocation_propensity(Propensity):
    def __init__(self, k_f : float, k_r : float, reactant : Species, R_max : float, n : float, 
                 c : Species, c_max : float, product : Species, reversible = True):
        """resource_allocation_reversible_propensity is a propensity with the following formula.
            p(c,x; k, c_max) = k_f * R_free * reactant - k_r * product
            where R_free = R_total - product and R_total = R_max * y^n where y = f * (1 - f)
        :param k_forward: forward rate constant (float)
        :param k_reverse: reverse rate constant (float)
        :param R_max: maximum resource capacity (float)
        """
        Propensity.__init__(self)
        self.k_f = k_f
        self.k_r = k_r
        self.R_max = R_max
        self.n = n
        self.c = c
        self.c_max = c_max
        self.reactant = reactant
        self.product = product
        self.reversible = reversible
        self.name = 'resource_allocation_reversible_propensity'

    @property
    def k_f(self):
        if isinstance(self._k_f, Parameter):
            return self._k_f.value
        else:
            return self._k_f

    @k_f.setter
    def k_f(self, new_k_f):
        self._k_f = self._check_parameter(new_k_f)
        self.propensity_dict['parameters']['k_f'] = self._k_f

    @property
    def k_r(self):
        if isinstance(self._k_r, Parameter):
            return self._k_r.value
        else:
            return self._k_r

    @k_r.setter
    def k_r(self, new_k_r):
        self._k_r = self._check_parameter(new_k_r)
        self.propensity_dict['parameters']['k_r'] = self._k_r

    @property
    def R_max(self):
        if isinstance(self._R_max, Parameter):
            return self._R_max.value
        else:
            return self._R_max

    @R_max.setter
    def R_max(self, new_R_max):
        self._R_max = self._check_parameter(new_R_max)
        self.propensity_dict['parameters']['R_max'] = self._R_max

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
    def product(self):
        return self._product

    @product.setter
    def product(self, new_product):
        self._product = self._check_species(new_product)
        self.propensity_dict['species']['product'] = self.product
    
    @property
    def reactant(self):
        return self._reactant

    @reactant.setter
    def reactant(self, new_reactant):
        self._reactant = self._check_species(new_reactant)
        self.propensity_dict['species']['reactant'] = self.reactant
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):
        if self.reversible:
            return f'rate = k_f * {self.reactant.pretty_print(**kwargs)} * (R_max * ({self.c.pretty_print(**kwargs)}/c_max * (1 - {self.c.pretty_print(**kwargs)}/c_max))^n - {self.product.pretty_print(**kwargs)}) - k_r * {self.product.pretty_print(**kwargs)}'
        else:
            return f'rate = k_f * {self.reactant.pretty_print(**kwargs)} * (R_max * ({self.c.pretty_print(**kwargs)}/c_max * (1 - {self.c.pretty_print(**kwargs)}/c_max))^n - {self.product.pretty_print(**kwargs)})'
    
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
        k_f = propensity_dict['parameters']['k_f']
        k_r = propensity_dict['parameters']['k_r']
        R_max = propensity_dict['parameters']['R_max']
        n = propensity_dict['parameters']['n']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        reactant = propensity_dict['species']['reactant']
        product = propensity_dict['species']['product']
        
        if self.reversible:
            rate_formula = f"{k_f} * {reactant} * ({R_max} * ({c}/{c_max} * (1 - {c}/{c_max}))^{n} - {product}) - {k_r} * {product}"
        
        else:
            rate_formula = f"{k_f} * {reactant} * ({R_max} * ({c}/{c_max} * (1 - {c}/{c_max}))^{n} - {product})"
       

        return rate_formula
    
# Resource allocation propensity with multiple products
class resource_allocation_propensity_2_products(Propensity):
    def __init__(self, k_f : float, k_r : float, reactant : Species, R_max : float, n : float, 
                 c : Species, c_max : float, product_1 : Species, product_2 : Species, product_case : Species):
        """resource_allocation_propensity_2_products is a propensity with the following formula.
            p(c,x; k, c_max) = k_f * R_free * reactant - k_r * product
            where R_free = R_total - products and R_total = R_max * y^n where y = f * (1 - f)
        :param k_forward: forward rate constant (float)
        :param k_reverse: reverse rate constant (float)
        :param R_max: maximum resource capacity (float)
        """
        Propensity.__init__(self)
        self.k_f = k_f
        self.k_r = k_r
        self.R_max = R_max
        self.n = n
        self.c = c
        self.c_max = c_max
        self.reactant = reactant
        self.product_1 = product_1
        self.product_2 = product_2
        self.product_case = product_case
        self.name = 'resource_allocation_propensity_2_products'

    @property
    def k_f(self):
        if isinstance(self._k_f, Parameter):
            return self._k_f.value
        else:
            return self._k_f

    @k_f.setter
    def k_f(self, new_k_f):
        self._k_f = self._check_parameter(new_k_f)
        self.propensity_dict['parameters']['k_f'] = self._k_f

    @property
    def k_r(self):
        if isinstance(self._k_r, Parameter):
            return self._k_r.value
        else:
            return self._k_r

    @k_r.setter
    def k_r(self, new_k_r):
        self._k_r = self._check_parameter(new_k_r)
        self.propensity_dict['parameters']['k_r'] = self._k_r

    @property
    def R_max(self):
        if isinstance(self._R_max, Parameter):
            return self._R_max.value
        else:
            return self._R_max

    @R_max.setter
    def R_max(self, new_R_max):
        self._R_max = self._check_parameter(new_R_max)
        self.propensity_dict['parameters']['R_max'] = self._R_max

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
    def product_1(self):
        return self._product_1

    @product_1.setter
    def product_1(self, new_product_1):
        self._product_1 = self._check_species(new_product_1)
        self.propensity_dict['species']['product_1'] = self.product_1

    @property
    def product_2(self):
        return self._product_2

    @product_2.setter
    def product_2(self, new_product_2):
        self._product_2 = self._check_species(new_product_2)
        self.propensity_dict['species']['product_2'] = self.product_2

    @property
    def product_case(self):
        return self._product_case

    @product_case.setter
    def product_case(self, new_product_case):
        self._product_case = self._check_species(new_product_case)
        self.propensity_dict['species']['product_case'] = self.product_case
    
    @property
    def reactant(self):
        return self._reactant

    @reactant.setter
    def reactant(self, new_reactant):
        self._reactant = self._check_species(new_reactant)
        self.propensity_dict['species']['reactant'] = self.reactant
    
    def pretty_print_rate(self, show_parameters = True, **kwargs):

        return f'rate = k_f * {self.reactant.pretty_print(**kwargs)} * (R_max * ({self.c.pretty_print(**kwargs)}/c_max * (1 - {self.c.pretty_print(**kwargs)}/c_max))^n - ({self.product_1.pretty_print(**kwargs)} + {self.product_2.pretty_print(**kwargs)})) - k_r * {self.product_case.pretty_print(**kwargs)}'
        
    
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
        k_f = propensity_dict['parameters']['k_f']
        k_r = propensity_dict['parameters']['k_r']
        R_max = propensity_dict['parameters']['R_max']
        n = propensity_dict['parameters']['n']
        c_max = propensity_dict['parameters']['c_max']
        c = propensity_dict['species']['c']
        reactant = propensity_dict['species']['reactant']
        product_1 = propensity_dict['species']['product_1']
        product_2 = propensity_dict['species']['product_2']
        product_case = propensity_dict['species']['product_case']
        
        rate_formula = f"{k_f} * {reactant} * ({R_max} * ({c}/{c_max} * (1 - {c}/{c_max}))^{n} - ({product_1} + {product_2})) - {k_r} * {product_case}"
        


       

        return rate_formula