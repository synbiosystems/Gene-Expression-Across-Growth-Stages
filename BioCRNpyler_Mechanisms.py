from biocrnpyler import *

from BioCRNpyler_Propensities import *

# Class for defining all the reactions in the transcription mechanism in the CRN
class BacterialTranscription(Mechanism):
    """Transcription with sigma factor and holoenzyme binding.
        sigma + RNAP <--> holoenzyme
        holoenzyme + DNA --> C_open (Open promoter complex)
        C_open + nucleotide --> DNA + RNAP + sigma + T (Transcription consumes nucleotides)
        ø <--> RNAP
    """

    def __init__(self, sigma: Species, holoenzyme : Species, nucleotide: Species, 
                 name = "bacterial_transcription", cell_count = None, **keywords):
        """Initializes a Transcription_MM instance.
        :param sigma: Species instance that represents sigma factor
        :param nucleotide: Species instance that represents nucleotides
        :param rnap: Species instance that is representing RNA polymerase
        :param name: name of the Mechanism, default: transcription_mm
        :param cell_count: Species representing cell count
        """
        
        if isinstance(sigma, Species):
            self.sigma = sigma
        else:
            raise ValueError("'sigma' argument must be a Species.")
        
        if isinstance(holoenzyme, Species):
            self.holoenzyme = holoenzyme
        else:
            raise ValueError("'holoenzyme' argument must be a Species.")
            
        if isinstance(nucleotide, Species):
            self.nucleotide = nucleotide
        else:
            raise ValueError("'NT' parameter must be a Species.")

        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")

        Mechanism.__init__(self=self, name=name,
                           mechanism_type="transcription")

    def update_species(self, dna, transcript = None, protein = None, **keywords):
        species = [dna, self.sigma, self.holoenzyme, self.nucleotide, transcript]
        open_complex = Complex([dna, self.holoenzyme], attribute = ["open"])
        species.append(open_complex)
        return species

    def update_reactions(self, dna, component, part_id = None, complex = None, transcript = None, protein = None,
                         **keywords):
        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name


        k_tx_1b = component.get_parameter("k_tx_1b", part_id = part_id, mechanism = self)
        k_tx_1u = component.get_parameter("k_tx_1u", part_id = part_id, mechanism = self)
        RNAP_max = component.get_parameter("RNAP_max", part_id = part_id, mechanism = self)
        n_gamma_RNAP = component.get_parameter("n_gamma_RNAP", part_id = part_id, mechanism = self)
        k_tx_2b = component.get_parameter("k_tx_2b", part_id = part_id, mechanism = self)
        k_tx_2u = component.get_parameter("k_tx_2u", part_id = part_id, mechanism = self)
        k_tx_3 = component.get_parameter("k_tx_3", part_id = part_id, mechanism = self)
        b_tx_4u = component.get_parameter("b_tx_4u", mechanism = "mrna_degradation")
        n_delta = component.get_parameter("n_delta", part_id = part_id, mechanism = self)
        c_max = component.get_parameter("c_max", part_id = part_id, mechanism = "logistic_cell_growth")
    
     
        open_complex = Complex([dna, self.holoenzyme], attributes = ["open"])
        
        rnap_sigma_binding_propensity = resource_allocation_propensity(k_f = k_tx_1b, k_r = k_tx_1u, reactant = self.sigma, 
                                                                                  R_max = RNAP_max, n = n_gamma_RNAP,
                                                                                  c = self.cell_count, c_max = c_max, 
                                                                                  product = self.holoenzyme)
        
        r1 = Reaction([self.sigma], [self.holoenzyme], propensity_type = rnap_sigma_binding_propensity)

        # holoenzyme and dna binding
        r2 = Reaction.from_massaction([self.holoenzyme, dna], [open_complex], k_forward = k_tx_2b)
        
        
        # reverse reaction (unbinding) of Open promoter to DNA + Holoenzyme. This reaction is favoured when cells 
        # are in stationary phase, as the Open promoter complex is destabilized

        propensity_open_complex = delta_with_basal(k = k_tx_2u, c = self.cell_count, c_max = c_max, x = open_complex, 
            b = b_tx_4u, positive = True, n = n_delta)
        r3 = Reaction([open_complex], [self.holoenzyme, dna], propensity_type = propensity_open_complex)
        
        
        # Transcription elongation
        r4 = Reaction.from_massaction([self.nucleotide, open_complex], [dna, transcript, self.sigma],
                                      k_forward = k_tx_3)
        
        
        return [r1, r2, r3, r4]
    

# Class for defining all the reactions for mRNA degradation mechanism in the CRN    
class BacterialDegradation_mRNA(GlobalMechanism, MichaelisMenten):
    """Michaelis Menten mRNA Degredation by Endonucleases that results in nucleotides being free.
       mRNA + Endo <--> mRNA:Endo --> Endo + NT
       All species of type "rna" are degraded by this mechanisms, including those inside of a ComplexSpecies.
       ComplexSpecies are seperated by this process, including embedded ComplexSpecies. 
       OrderedPolymerSpecies are ignored.
    """
    def __init__(self, nuclease, nucleotide, 
                 name = "mrna_degradation", 
                 cell_count = None, mechanism_type = "mrna_degradation", 
                 default_on = False, recursive_species_filtering = True, 
                 filter_dict = None, **keywords):

        if isinstance(nuclease, Species):
            self.nuclease = nuclease
        else:
            raise ValueError("'nuclease' must be a Species.")
            
        if isinstance(nucleotide, Species):
            self.nucleotide = nucleotide
        else:
            raise ValueError("'nucleotide' must be a Species.")
        
        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")
            
        MichaelisMenten.__init__(self=self, name=name, mechanism_type = mechanism_type)

        if filter_dict is None:
            filter_dict = {"rna":True, "notdegradable":False}

        GlobalMechanism.__init__(self, name = name, mechanism_type = mechanism_type, default_on = default_on,
                                 filter_dict = filter_dict, recursive_species_filtering = recursive_species_filtering)

    def update_species(self, s, mixture):
        species = [s, self.nuclease, self.nucleotide]
        #RNA species that are inside a ComplexSpecies are not degraded in this mechanism
        #So if the material type is simply RNA, break it up.
        complex_deg = Complex([s, self.nuclease])
        species.append(complex_deg)

        return species

    def update_reactions(self, s, mixture):
        
        rxns = []

        #If the material type is simply RNA, break it up.
        if s.material_type == "rna":
            k_tx_4b = self.get_parameter(s, "k_tx_4b", mixture)
            k_tx_4u = self.get_parameter(s, "k_tx_4u", mixture)
            k_tx_5 = self.get_parameter(s, "k_tx_5", mixture)
            b_tx_4u = self.get_parameter(s, "b_tx_4u", mixture)

            n_delta = mixture.get_parameter(param_name = "n_delta", part_id = None, mechanism = "mrna_degradation")
            c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
            
            complex_deg = Complex([s, self.nuclease])
            
            r1 = Reaction.from_massaction([s, self.nuclease], [complex_deg], k_forward = k_tx_4b)
            rxns.append(r1)

            propensity_complex = delta_with_basal(k = k_tx_4u, c = self.cell_count, c_max = c_max, x = complex_deg, 
            b = b_tx_4u, positive = True, n = n_delta)

            r2 = Reaction([complex_deg], [s, self.nuclease], propensity_type = propensity_complex)
            rxns.append(r2)

            r3 = Reaction.from_massaction([complex_deg], [self.nuclease, self.nucleotide], k_forward = k_tx_5)
            rxns.append(r3)
            
        return rxns


# Class for all the Translation related reactions in the CRN
class BacterialTranslation(Mechanism):
    """Translation in growing bacteria that consumes amino acids and ATP.

        Et + AA <--> Ct
        Ct + tRNA <--> Caa + Et
        mRNA + Ribo <--> Ctic
        Caa + Ctic --> mRNA + Ribo + tRNA + unfolded protein 
        unfolded protein --> folded protein 
        ø <--> tRNA
        ø <--> Ribosome
        ø <--> tRNA synthetase
        ø --> AA
        ø <--> Protease
        unfolded protein + Protease --> C_deg --> Pc + Protease
        folded protein + Protease --> C_deg --> Pc + Protease

    """

    def __init__(self, aminoacid: Species,
                 unfolded_protein: Species, peptide_chain: Species,
                 Ct : Species, Caa : Species, Ctic : Species, 
                 C_deg_unfolded : Species, C_deg_folded : Species, 
                 name = "bacterial_translation", degtag = True, cell_count = None, **keywords):
        """Initializes a Translation instance.
        
        :param ribosome: Species instance that is representing a ribosome
        :param name: name of the Mechanism, default: energy_translation_mm
        :param cell_count: Species representing cell count
        """
        if isinstance(aminoacid, Species):
            self.aminoacid = aminoacid
        else:
            raise ValueError("aminoacid must be a Species!")
        
        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")
            
        if isinstance(unfolded_protein, Species):
            self.unfolded_protein = unfolded_protein
        else:
            raise ValueError("'unfolded_protein' parameter must be a Species.")
            
        if isinstance(peptide_chain, Species):
            self.peptide_chain = peptide_chain
        else:
            raise ValueError("'peptide_chain' parameter must be a Species.")
        
        if isinstance(Ct, Species):
            self.Ct = Ct
        else:
            raise ValueError("'Ct' parameter must be a Species.")
        
        if isinstance(Caa, Species):
            self.Caa = Caa
        else:
            raise ValueError("'Caa' parameter must be a Species.")
        
        if isinstance(Ctic, Species):
            self.Ctic = Ctic
        else:
            raise ValueError("'Ctic' parameter must be a Species.")
        
        if isinstance(C_deg_unfolded, Species):
            self.C_deg_unfolded = C_deg_unfolded
        else:
            raise ValueError("'C_deg_unfolded' parameter must be a Species.")
        
        if isinstance(C_deg_folded, Species):
            self.C_deg_folded = C_deg_folded
        else:
            raise ValueError("'C_deg_folded' parameter must be a Species.")
        
        self.degtag = degtag
        
            
            
        Mechanism.__init__(self = self, name=name, mechanism_type = "translation")

    def update_species(self, transcript, protein, **keywords):
        species = [self.aminoacid, protein, self.unfolded_protein, 
                  self.peptide_chain, self.Ct, self.Caa, 
                  self.Ctic, self.C_deg_unfolded, self.C_deg_folded]
        
        return species

    def update_reactions(self, transcript, protein, component, part_id = None, complex=None, **keywords):
        rxns = []

        # Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        #k_charge = component.get_parameter("k_charge", part_id = part_id, mechanism = self)
        k_tl_1b = component.get_parameter("k_tl_1b", part_id = part_id, mechanism = self) # AA + Et binding
        k_tl_1u = component.get_parameter("k_tl_1u", part_id = part_id, mechanism = self) # AA + Et unbinding
        k_tl_2b = component.get_parameter("k_tl_2b", part_id = part_id, mechanism = self) # Ct + tRNA binding
        k_tl_2u = component.get_parameter("k_tl_2u", part_id = part_id, mechanism = self) # Ct + tRNA unbinding
        k_tl_3b = component.get_parameter("k_tl_3b", part_id = part_id, mechanism = self) # mRNA Ribosome Binding
        k_tl_3u = component.get_parameter("k_tl_3u", part_id = part_id, mechanism = self) # mRNA Ribosome Unbinding
        k_tl_4 = component.get_parameter("k_tl_4", part_id = part_id, mechanism = self) # TL elongation
        k_tl_5 = component.get_parameter("k_tl_5", part_id = part_id, mechanism = self) # Protein folding
        b_tl_5 = component.get_parameter("b_tl_5", part_id = part_id, mechanism = self) # Protein folding basal
        k_tl_6 = component.get_parameter("k_tl_6", part_id = part_id, mechanism = self) # C_deg degradation
        k_tl_8 = component.get_parameter("k_tl_8", part_id = part_id, mechanism = self) # Amino acid synthesis
        k_tl_9b_P = component.get_parameter("k_tl_9b_P", part_id = part_id, mechanism = self) # Protease binding rate for unfolded Protein
        k_tl_9b_Pm = component.get_parameter("k_tl_9b_Pm", part_id = part_id, mechanism = self) # Protease binding rate for folded Protein
        k_tl_9u = component.get_parameter("k_tl_9u", part_id = part_id, mechanism = self) # Protease unbinding
        k_tl_10 = component.get_parameter("k_tl_10", part_id = part_id, mechanism = self) # Peptide chain degradation

        Et_max = component.get_parameter("Et_max", part_id = part_id, mechanism = self) # Max conc. of Et
        n_gamma_Et = component.get_parameter("n_gamma_Et", part_id = part_id, mechanism = self) # Exponent of gamma for Et
        tRNA_max = component.get_parameter("tRNA_max", part_id = part_id, mechanism = self) # Max conc. of tRNA
        n_gamma_tRNA = component.get_parameter("n_gamma_tRNA", part_id = part_id, mechanism = self) # Exponent of gamma for tRNA
        Ribo_max = component.get_parameter("Ribo_max", part_id = part_id, mechanism = self) # Max conc. of Ribosome
        n_gamma_Ribo = component.get_parameter("n_gamma_Ribo", part_id = part_id, mechanism = self) # Exponent of gamma for Ribosome
        Protease_max = component.get_parameter("Protease_max", part_id = part_id, mechanism = self) # Max conc. of Protease
        n_gamma_Protease = component.get_parameter("n_gamma_Protease", part_id = part_id, mechanism = self) # Exponent of gamma for Protease
        n_gamma_folding = component.get_parameter("n_gamma_folding", part_id = part_id, mechanism = self) # Exponent of gamma for Protein folding
        n_gamma_syn = component.get_parameter("n_gamma_syn", part_id = part_id, mechanism = self) # Exponent of gamma for Amino acid synthesis

        c_max = component.get_parameter("c_max", part_id = part_id, mechanism = "logistic_cell_growth")
       

        rxns  = []
        
        # Amino acid-tRNA synthetase reaction

        # sigma and free RNAP Binding
        aa_Et_binding_propensity = resource_allocation_propensity(k_f = k_tl_1b, k_r = k_tl_1u, reactant = self.aminoacid, 
                                                                                  R_max = Et_max, n = n_gamma_Et,
                                                                                  c = self.cell_count, c_max = c_max,
                                                                                  product = self.Ct)
        
        r1 = Reaction([self.aminoacid], [self.Ct], propensity_type = aa_Et_binding_propensity)

        rxns.append(r1)

        # Ct tRNA reaction

        tRNA_Ct_binding_propensity = resource_allocation_propensity(k_f = k_tl_2b, k_r = k_tl_2u, reactant = self.Ct, 
                                                                                  R_max = tRNA_max, n = n_gamma_tRNA,
                                                                                  c = self.cell_count, c_max = c_max,
                                                                                  product = self.Caa, reversible = False)
        
        r2 = Reaction([self.Ct], [self.Caa], propensity_type = tRNA_Ct_binding_propensity)

        rxns.append(r2)
        
        # Translation complex reaction

        mRNA_Ribosome_binding_propensity = resource_allocation_propensity(k_f = k_tl_3b, k_r = k_tl_3u, reactant = transcript, 
                                                                                  R_max = Ribo_max, n = n_gamma_Ribo,
                                                                                  c = self.cell_count, c_max = c_max,
                                                                                  product = self.Ctic)
        
        r3 = Reaction([transcript], [self.Ctic], propensity_type = mRNA_Ribosome_binding_propensity)

        rxns.append(r3)

        # Translation elongation reaction

        r4 = Reaction.from_massaction([self.Caa, self.Ctic], [transcript, self.unfolded_protein], k_forward = k_tl_4)
    
        rxns.append(r4)
        
        # Protein folding

        folding_propensity = apply_gamma_propensity(k = k_tl_5, x = self.unfolded_protein, b = b_tl_5, c = self.cell_count,
                                                        c_max = c_max, n = n_gamma_folding)

        r5 = Reaction([self.unfolded_protein], [protein], propensity_type = folding_propensity)
        rxns.append(r5)
        
        
        # Amino acid synthesis

        tl_aa_synthesis_propensity = gamma(k = k_tl_8, c = self.cell_count, c_max = c_max, n = n_gamma_syn)

        r6 = Reaction([], [self.aminoacid], propensity_type = tl_aa_synthesis_propensity)
        rxns.append(r6)

        # Peptide to Amino Acid

        r7 = Reaction.from_massaction([self.peptide_chain], [self.aminoacid],
                                      k_forward = k_tl_10)
        rxns.append(r7)

        ## degradation due to degradation tag 

        if self.degtag: 

            # Unfolded peptide degradation by protease
            unfolded_P_protease_binding_propensity = resource_allocation_propensity_2_products(k_f = k_tl_9b_P, k_r = k_tl_9u, reactant = self.unfolded_protein, 
                                                                                  R_max = Protease_max, n = n_gamma_Protease,
                                                                                  c = self.cell_count, c_max = c_max,
                                                                                  product_1 = self.C_deg_unfolded, product_2 = self.C_deg_folded, 
                                                                                  product_case = self.C_deg_unfolded)
        
            r8 = Reaction([self.unfolded_protein], [self.C_deg_unfolded], propensity_type = unfolded_P_protease_binding_propensity)

            rxns.append(r8)

            r9 = Reaction.from_massaction([self.C_deg_unfolded], [self.peptide_chain], k_forward = k_tl_6)
            rxns.append(r9)

            # Folded peptide degradation by protease
            folded_P_protease_binding_propensity = resource_allocation_propensity_2_products(k_f = k_tl_9b_Pm, k_r = k_tl_9u, reactant = protein, 
                                                                                  R_max = Protease_max, n = n_gamma_Protease,
                                                                                  c = self.cell_count, c_max = c_max,
                                                                                  product_1 = self.C_deg_unfolded, product_2 = self.C_deg_folded, 
                                                                                  product_case = self.C_deg_folded)
        
            r10 = Reaction([protein], [self.C_deg_folded], propensity_type = folded_P_protease_binding_propensity)

            rxns.append(r10)

            r11 = Reaction.from_massaction([self.C_deg_folded], [self.peptide_chain], k_forward = k_tl_6)
            rxns.append(r11)
        
      
        return rxns



# This class contains all the non-deg tag protein degradation related reactions for the CRN 
class BacterialDegradation_active(GlobalMechanism):
    """Non deg-tag degradation of proetin or the endogenous protein degradation 
    """
    def __init__(self,aminoacid, peptide_chain, name="non_tag_degradation", 
                 cell_count = None, mechanism_type = "non_tag_degradation", 
                 filter_dict= None, recursive_species_filtering = False, default_on = False, **keywords):
        
        if isinstance(aminoacid, Species):
            self.aminoacid = aminoacid
        else:
            raise ValueError("'aminoacid' must be a Species.")
        
        if isinstance(peptide_chain, Species):
            self.peptide_chain = peptide_chain
        else:
            raise ValueError("'peptide_chain' must be a Species.")
            
        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")

        
            
        GlobalMechanism.__init__(self, name = name, mechanism_type = mechanism_type, default_on = default_on,
                                 filter_dict = filter_dict, recursive_species_filtering = recursive_species_filtering)

    def update_species(self, s, mixture):
        species = [s, self.aminoacid, self.peptide_chain]
        return species

    def update_reactions(self, s, mixture):
        if s.material_type == "protein": 

            k_tl_7 = mixture.get_parameter(param_name = "k_tl_7", part_id = None, mechanism = "non_tag_degradation")
            n_delta = mixture.get_parameter(param_name = "n_delta", part_id = None, mechanism = "non_tag_degradation")
            c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
            rxns = []
            
            propensity_dp = delta_propensity(k = k_tl_7, c = self.cell_count, c_max = c_max, x = s, 
        positive = True, n = n_delta)
            
            r1 = Reaction([s], [self.peptide_chain], propensity_type = propensity_dp)
            rxns.append(r1)
            
        return rxns



# Create a Dilution mechanism appropriate for the cell growth mechanism: 

class DilutionLogGrowth(GlobalMechanism):
    """A global mechanism to represent dilution based on Logistic Growth profile."""

    def __init__(self, name = "global_dilution_logistic_growth",
                 mechanism_type = "dilution", filter_dict = None, cell_count = None,
                 default_on = True, recursive_species_filtering = True):
        GlobalMechanism.__init__(self, name = name,
                                 mechanism_type = mechanism_type,
                                 default_on = default_on,
                                 filter_dict = filter_dict,
                                 recursive_species_filtering = recursive_species_filtering)
        
        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")


    def update_reactions(self, s: Species, mixture):
        
        k_gr = mixture.get_parameter(param_name = "k", part_id = None, mechanism = "logistic_cell_growth")
        c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
        rxn = []
        if isinstance(s, ComplexSpecies):
            dilution_propensity = alpha_propensity(k = k_gr, c = self.cell_count, c_max = c_max, x = s)
            r = Reaction([s], [], propensity_type = dilution_propensity)
            rxn.append(r)
        else:
            dilution_propensity = alpha_propensity(k = k_gr, c = self.cell_count, c_max = c_max, x = s)
            r = Reaction([s], [], propensity_type = dilution_propensity)
            rxn.append(r)
        return rxn

# This clas defines the growth law reactions 
class LogisticCellGrowth(Mechanism):
    """Cell growth modeled as logistic growth function using LogisticPositive propensity:
     --> cell_count, propensity = k*cell_count*(1 - cell_count/c_max)
    """

    def __init__(self, name = "logistic_cell_growth", cell_count = None, 
                 **keywords):
        """Initializes a LogisticCellGrowth instance.
        :param cell_count: Species representing cell count
        """
        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' must be a Species.")
        Mechanism.__init__(self=self, name=name, mechanism_type="logistic_cell_growth")

    def update_species(self, **keywords):
        
        species = [self.cell_count]
        return species

    def update_reactions(self, component, part_id = None, **keywords):
        
        #Get Parameters
        rxns = []
        #Get Parameters
        k = component.get_parameter(param_name = "k", part_id = part_id, mechanism = self)
        c_max = component.get_parameter(param_name = "c_max", part_id = part_id, mechanism = self)
        c0 = component.get_parameter(param_name = "c0", part_id = part_id, mechanism = self)
        cell_growth_propensity = LogisticPositive(k = k, c = self.cell_count, c_max = c_max)
        r0 = Reaction([],[self.cell_count], propensity_type = cell_growth_propensity)
        return [r0]