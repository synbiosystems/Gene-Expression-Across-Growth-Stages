from biocrnpyler import *

from GEAGS_propensities_4 import *

# Class for defining all the reactions in the transcription mechanism in the CRN
class BacterialTranscription(Mechanism):
    """Transcription with sigma factor and holoenzyme binding.
        sigma + RNAP <--> holoenzyme
        holoenzyme + DNA --> C_open (Open promoter complex)
        C_open + nucleotide --> DNA + RNAP + T (Transcription consumes nucleotides)
    """

    def __init__(self, rnap: Species, sigma: Species, nucleotide: Species, 
                 name="bacterial_transcription", cell_count = None, **keywords):
        """Initializes a Transcription_MM instance.
        :param sigma: Species instance that represents sigma factor
        :param nucleotide: Species instance that represents nucleotides
        :param rnap: Species instance that is representing an RNA polymerase
        :param name: name of the Mechanism, default: transcription_mm
        :param cell_count: Species representing cell count
        """
        if isinstance(rnap, Species):
            self.rnap = rnap
        else:
            raise ValueError("'rnap' argument must be a Species.")
        
        if isinstance(sigma, Species):
            self.sigma = sigma
        else:
            raise ValueError("'sigma' argument must be a Species.")
            
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
        species = [dna, self.rnap, self.sigma, self.nucleotide, transcript]
        holoenzyme = Complex([self.sigma, self.rnap])
        species.append(holoenzyme)
        open_complex = Complex([dna, holoenzyme], attribute = ["open"])
        species.append(open_complex)
        return species

    def update_reactions(self, dna, component, part_id = None, complex=None, transcript=None, protein=None,
                         **keywords):
        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name


        k_tx_1b = component.get_parameter("k_tx_1b", part_id = part_id, mechanism = self)
        k_tx_1u = component.get_parameter("k_tx_1u", part_id = part_id, mechanism = self)
        k_tx_2b = component.get_parameter("k_tx_2b", part_id = part_id, mechanism = self)
        k_tx_2u = component.get_parameter("k_tx_2u", part_id = part_id, mechanism = self)
        k_tx_3 = component.get_parameter("k_tx_3", part_id = part_id, mechanism = self)
        # b_tx_4u = component.get_parameter("b_tx_4u", mechanism = "mrna_degradation")
        n = component.get_parameter("n", part_id = part_id, mechanism = self)
        c_max = component.get_parameter("c_max", part_id = part_id, mechanism = "logistic_cell_growth")
        c0 = component.get_parameter("c0", part_id = part_id, mechanism = "logistic_cell_growth")
        
        holoenzyme = Complex([self.sigma, self.rnap])
        open_complex = Complex([dna, holoenzyme], attributes = ["open"])
        
        # sigma and RNAP Binding
        r1 = Reaction.from_massaction([self.sigma, self.rnap], [holoenzyme], k_forward = k_tx_1b, k_reverse = k_tx_1u)
        
        
        # holoenzyme and dna binding
        r2 = Reaction.from_massaction([holoenzyme, dna], [open_complex], k_forward = k_tx_2b)
        
        
        # reverse reaction (unbinding) of Open promoter to DNA + Holoenzyme. This reaction is favoured when cells 
        # are in stationary phase, as the Open promoter complex is destabilized

        propensity_open_complex = hill_growth(k = k_tx_2u, c = self.cell_count, c_max = c_max, 
                                              x = open_complex, n = n, positive = True)
        # propensity_open_complex = hill_growth_basal(k = k_tx_2u, c = self.cell_count, c_max = c_max, x = open_complex, 
        #     b = b_tx_4u, positive = True, n = n)
        r3 = Reaction([open_complex], [holoenzyme, dna], propensity_type = propensity_open_complex)
        
        
        # Transcription elongation
        r4 = Reaction.from_massaction([self.nucleotide, open_complex], [dna, self.rnap, transcript, self.sigma],
                                      k_forward = k_tx_3)
        
        # RNAP with growth
        k_rnap = component.get_parameter("k_rnap", part_id = part_id, mechanism = self)
        rnap_min = component.get_parameter("rnap_min", part_id = part_id, mechanism = self)
        
        propensity_rnap = three_species_source(k = k_rnap, c = self.cell_count, c_max = c_max, S1 = self.rnap, S2 = holoenzyme,
                                               S3 = open_complex, S_min = rnap_min)
    
        r5 = Reaction([], [self.rnap], propensity_type = propensity_rnap)

        # # Sigma with growth
        # k_sigma = component.get_parameter("k_sigma", part_id = part_id, mechanism = self)
        # sigma_min = component.get_parameter("sigma_min", part_id = part_id, mechanism = self)
        
        # propensity_sigma = three_species_hill_source(k = k_sigma, c = self.cell_count, c_max = c_max, S1 = self.sigma, S2 = holoenzyme,
        #                                        S3 = open_complex, S_min = sigma_min, n = n)
    
        # r6 = Reaction([self.sigma],[], propensity_type = propensity_sigma)
        
        
        return [r1, r2, r3, r4, r5]
    

# Class for defining all the reactions for mRNA degradation mechanism in the CRN    
class BacterialDegradation_mRNA(GlobalMechanism, MichaelisMenten):
    """Michaelis Menten mRNA Degredation by Endonucleases that results in nucleotides being free.
       mRNA + Endo <--> mRNA:Endo --> Endo + NT
       All species of type "rna" are degraded by this mechanisms, including those inside of a ComplexSpecies.
       ComplexSpecies are seperated by this process, including embedded ComplexSpecies. 
       OrderedPolymerSpecies are ignored.
    """
    def __init__(self, nuclease, nucleotide, 
                 name="mrna_degradation", 
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
            n = self.get_parameter(s, "n", mixture)
            c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
            
            complex_deg = Complex([s, self.nuclease])
            
            r1 = Reaction.from_massaction([s, self.nuclease], [complex_deg], k_forward = k_tx_4b)
            rxns.append(r1)

            propensity_complex = hill_growth_basal(k = k_tx_4u, c = self.cell_count, c_max = c_max, x = complex_deg, 
            b = b_tx_4u, positive = True, n = n)

            r2 = Reaction([complex_deg], [s, self.nuclease], propensity_type = propensity_complex)
            rxns.append(r2)

            r3 = Reaction.from_massaction([complex_deg], [self.nuclease, self.nucleotide], k_forward = k_tx_5)
            rxns.append(r3)
            
        return rxns


# Class for all the Translation related reactions in the CRN
class BacterialTranslation(Mechanism):
    """Translation in growing bacteria that consumes amino acids and ATP.

        Et + AA <--> Ct
        Ct --> Ct-amp
        Ct-amp + tRNA --> Cin 
        Cin --> Cac
        Cac --> Et + Caa
        Caa + mRNA + Ribo <--> Ctic
        Ctic --> mRNA + Ribo + tRNA + unfolded protein 
        unfolded protein --> folded protein 

    """

    def __init__(self, ribosome: Species, aminoacid: Species, tRNA: Species, 
                 unfolded_protein: Species, Et: Species, peptide_chain: Species,
                 protease: Species, 
                 name="bacterial_translation", degtag = True, cell_count = None, **keywords):
        """Initializes a Translation instance.
        
        :param ribosome: Species instance that is representing a ribosome
        :param name: name of the Mechanism, default: energy_translation_mm
        :param cell_count: Species representing cell count
        """
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        else:
            raise ValueError("ribosome must be a Species!")
        
        if isinstance(aminoacid, Species):
            self.aminoacid = aminoacid
        else:
            raise ValueError("aminoacid must be a Species!")
        
        if isinstance(tRNA, Species):
            self.tRNA = tRNA
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
            
        if isinstance(Et, Species):
            self.Et = Et
        else:
            raise ValueError("'Et' parameter must be a Species.")

        if isinstance(peptide_chain, Species):
            self.peptide_chain = peptide_chain
        else:
            raise ValueError("'peptide_chain' parameter must be a Species.")
        
        if isinstance(protease, Species):
            self.protease = protease
        else:
            raise ValueError("'protease' parameter must be a Species.")
        
        self.degtag = degtag
        
            
            
        Mechanism.__init__(self = self, name=name, mechanism_type = "translation")

    def update_species(self, transcript, protein, **keywords):
        species = [self.ribosome, self.aminoacid, self.tRNA, protein, self.unfolded_protein, 
                  self.Et, self.peptide_chain, self.protease]
        

        Complex_t = Complex([self.aminoacid, self.Et])
        species.append(Complex_t)

        Complex_t_amp = Complex([self.aminoacid, self.Et], attribute = "amp")
        species.append(Complex_t_amp)

        Complex_inactive = Complex([Complex_t_amp, self.tRNA], attribute = "inactive")
        species.append(Complex_inactive)

        Complex_active = Complex([Complex_t_amp, self.tRNA], attribute = "active")
        species.append(Complex_active)

        Complex_charged = Complex([self.aminoacid, self.tRNA])
        species.append(Complex_charged)

        Complex_tic = Complex([Complex_charged, transcript, self.ribosome])
        species.append(Complex_tic)

        deg_complex_unfolded = Complex([self.unfolded_protein, self.protease])
        species.append(deg_complex_unfolded)

        deg_complex_folded = Complex([protein, self.protease])
        species.append(deg_complex_folded)
        
        return species

    def update_reactions(self, transcript, protein, component, part_id = None, complex=None, **keywords):
        rxns = []

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        #k_charge = component.get_parameter("k_charge", part_id = part_id, mechanism = self)
        k_tl_1b = component.get_parameter("k_tl_1b", part_id = part_id, mechanism = self)
        k_tl_1u = component.get_parameter("k_tl_1u", part_id = part_id, mechanism = self)
        k_tl_2 = component.get_parameter("k_tl_2", part_id = part_id, mechanism = self)
        k_tl_3 = component.get_parameter("k_tl_3", part_id = part_id, mechanism = self)
        k_tl_4 = component.get_parameter("k_tl_4", part_id = part_id, mechanism = self)
        k_tl_5 = component.get_parameter("k_tl_5", part_id = part_id, mechanism = self)
        k_tl_6b = component.get_parameter("k_tl_6b", part_id = part_id, mechanism = self)
        k_tl_6u = component.get_parameter("k_tl_6u", part_id = part_id, mechanism = self)
        k_tl_7 = component.get_parameter("k_tl_7", part_id = part_id, mechanism = self)
        b_tl_7 = component.get_parameter("b_tl_7", part_id = part_id, mechanism = self)
        k_tl_8 = component.get_parameter("k_tl_8", part_id = part_id, mechanism = self)
        P_max = component.get_parameter("P_max", part_id = part_id, mechanism = self)
        k_tl_9 = component.get_parameter("k_tl_9", part_id = part_id, mechanism = self)
        b_tl_8 = component.get_parameter("b_tl_8", part_id = part_id, mechanism = self)
        k_tl_11 = component.get_parameter("k_tl_11", part_id = part_id, mechanism = self)
        k_tl_12 = component.get_parameter("k_tl_12", part_id = part_id, mechanism = self)
        k_tl_13b = component.get_parameter("k_tl_13b", part_id = part_id, mechanism = self)
        k_tl_13u = component.get_parameter("k_tl_13u", part_id = part_id, mechanism = self)
        
        c_max = component.get_parameter("c_max", part_id = part_id, mechanism = "logistic_cell_growth")
       
        k_gr = component.get_parameter("k", part_id = part_id, mechanism = "logistic_cell_growth")
        
       
        Complex_t = Complex([self.aminoacid, self.Et])

        Complex_t_amp = Complex([self.aminoacid, self.Et], attributes = ["amp"])

        Complex_inactive = Complex([Complex_t_amp, self.tRNA], attributes = ["inactive"])

        Complex_active = Complex([Complex_t_amp, self.tRNA], attributes = ["active"])

        Complex_charged = Complex([self.aminoacid, self.tRNA])

        Complex_tic = Complex([Complex_charged, transcript, self.ribosome])

        deg_complex_unfolded = Complex([self.unfolded_protein, self.protease])

        deg_complex_folded = Complex([protein, self.protease])

        rxns  = []
        
        # Amino acid-tRNA synthetase reaction

        r1 = Reaction.from_massaction([self.aminoacid, self.Et], [Complex_t], 
                                      k_forward = k_tl_1b, k_reverse = k_tl_1u)
        rxns.append(r1)

        # AA-synthetase reaction

        r2 = Reaction.from_massaction([Complex_t], [Complex_t_amp], k_forward = k_tl_2)
        rxns.append(r2)
        
        # AA-amp complex reaction

        r3 = Reaction.from_massaction([Complex_t_amp, self.tRNA], [Complex_inactive], k_forward = k_tl_3)
        rxns.append(r3)
        
        # aminoacyl-tRNA active complex reaction

        r4 = Reaction.from_massaction([Complex_inactive], [Complex_active], k_forward = k_tl_4)
        rxns.append(r4)
        
        # Charged complex reaction

        r5 = Reaction.from_massaction([Complex_active], [Complex_charged, self.Et], k_forward = k_tl_5)
        rxns.append(r5)
        
        # Translation complex reaction

        r6 = Reaction.from_massaction([Complex_charged, self.ribosome, transcript], [Complex_tic],
                                      k_forward = k_tl_6b, k_reverse = k_tl_6u)
        rxns.append(r6)

        # Translation elongation reaction

        tl_elongation_propensity = decrease_linear_growth_basal(k = k_tl_7, c = self.cell_count, c_max = c_max, x = Complex_tic, 
                                                            b = b_tl_7)
        r7 = Reaction([Complex_tic], [transcript, self.ribosome, self.unfolded_protein, self.tRNA], 
        propensity_type = tl_elongation_propensity)

        # r7 = Reaction.from_massaction([Complex_tic], [transcript, self.ribosome, self.unfolded_protein, self.tRNA],
        #                               k_forward = k_tl_7)
        rxns.append(r7)
        
        # Protein folding

        folding_propensity = protein_folding_propensity(k = k_tl_8, P_max = P_max, x = self.unfolded_protein, b = b_tl_8, c = self.cell_count,
                                                        c_max = c_max, p = protein)

        r8 = Reaction([self.unfolded_protein], [protein], propensity_type = folding_propensity)
        rxns.append(r8)
        

        
        # tRNA with growth
        k_tRNA = component.get_parameter("k_tRNA", part_id = part_id, mechanism = self)
        tRNA_min = component.get_parameter("tRNA_min", part_id = part_id, mechanism = self)

        propensity_trna = five_species_source(k = k_tRNA, c = self.cell_count, c_max = c_max, S1 = self.tRNA, S2 = Complex_inactive,
                                               S3 = Complex_active, S4 = Complex_charged, S5 = Complex_tic, S_min = tRNA_min)

        r9 = Reaction([], [self.tRNA], propensity_type = propensity_trna)
        rxns.append(r9)
        
        # Ribosome with growth
        k_ribo = component.get_parameter("k_ribo", part_id = part_id, mechanism = self)
        ribo_min = component.get_parameter("ribo_min", part_id = part_id, mechanism = self)
        
        propensity_ribo = two_species_source(k = k_ribo, c = self.cell_count, c_max = c_max, S1 = self.ribosome, S2 = Complex_tic,
                                               S_min = ribo_min)
        
        r10 = Reaction([], [self.ribosome], propensity_type = propensity_ribo)
        rxns.append(r10)
        
        
        # tRNA synthetase with growth
        k_Et = component.get_parameter("k_Et", part_id = part_id, mechanism = self)
        Et_min = component.get_parameter("Et_min", part_id = part_id, mechanism = self)

        propensity_Et = five_species_source(k = k_Et, c = self.cell_count, c_max = c_max, S1 = self.Et, S2 = Complex_t,
                                               S3 = Complex_t_amp, S4 = Complex_inactive, S5 = Complex_active, S_min = Et_min)
      
        r11 = Reaction([], [self.Et], propensity_type = propensity_Et)
        rxns.append(r11)
        

        # Peptide to Amino Acid

        r12 = Reaction.from_massaction([self.peptide_chain], [self.aminoacid],
                                      k_forward = k_tl_11)
        rxns.append(r12)
        
        # Amino acid replenishment 

        tl_aa_replenish_propensity = decrease_linear_growth_blank(k = k_tl_12, c = self.cell_count, c_max = c_max)

        r13 = Reaction([], [self.aminoacid], propensity_type = tl_aa_replenish_propensity)
        rxns.append(r13)

        # Protease with growth
        k_protease = component.get_parameter("k_protease", part_id = part_id, mechanism = self)
        protease_min = component.get_parameter("protease_min", part_id = part_id, mechanism = self)
        
        propensity_protease = three_species_source(k = k_protease, c = self.cell_count, c_max = c_max, S1 = self.protease, S2 = deg_complex_unfolded,
                                               S3 = deg_complex_folded, S_min = protease_min)
    
        r14 = Reaction([], [self.protease], propensity_type = propensity_protease)
        rxns.append(r14)

        if self.degtag: 

            r15 = Reaction.from_massaction([self.unfolded_protein, self.protease], [deg_complex_unfolded], k_forward = k_tl_13b, k_reverse = k_tl_13u)
            rxns.append(r15)

            r16 = Reaction.from_massaction([deg_complex_unfolded], [self.peptide_chain, self.protease], k_forward = k_tl_9)
            rxns.append(r16)

            r17 = Reaction.from_massaction([protein, self.protease], [deg_complex_folded], k_forward = k_tl_13b, k_reverse = k_tl_13u)
            rxns.append(r17)

            r18 = Reaction.from_massaction([deg_complex_folded], [self.peptide_chain, self.protease], k_forward = k_tl_9)
            rxns.append(r18)
        
      
        return rxns


class AminoAcidPool(GlobalMechanism):
    """A global mechanism to model the amino acid maintainence by degradation"""

    def __init__(self, aminoacid: Species, peptide_chain: Species, unfolded_protein: Species, 
                 name = "amino_acid_pool_maintainence", degtag = True,
                 mechanism_type = "amino_acid_pool_maintainence", filter_dict = None, cell_count = None,
                 default_on = True, recursive_species_filtering = True):
        GlobalMechanism.__init__(self, name = name,
                                 mechanism_type = mechanism_type,
                                 default_on = default_on,
                                 filter_dict = filter_dict,
                                 recursive_species_filtering = recursive_species_filtering)
        
        if isinstance(aminoacid, Species):
            self.aminoacid = aminoacid
        else:
            raise ValueError("aminoacid must be a Species!")
        
        if isinstance(peptide_chain, Species):
            self.peptide_chain = peptide_chain
        else:
            raise ValueError("peptide_chain must be a Species!")

        if isinstance(unfolded_protein, Species):
            self.unfolded_protein = unfolded_protein
        else:
            raise ValueError("unfolded_protein must be a Species!")

        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")

        self.degtag = degtag

    def update_reactions(self, s: Species, mixture):

        k_tl_9 = self.get_parameter(s, "k_tl_9", mixture)
        k_tl_10 = self.get_parameter(s, "k_tl_10", mixture)
        a_tl_10 = self.get_parameter(s, "a_tl_10", mixture)
        
        n = self.get_parameter(s, "n", mixture)
        
        c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
       
        rxn = []

        if self.degtag:

                deg_tag_propensity = hill_growth(k = k_tl_9, c = self.cell_count, c_max = c_max, x = s, 
                positive = False, n = n)

                r1 = Reaction([], [self.peptide_chain], propensity_type = deg_tag_propensity)
                rxn.append(r1)

        non_deg_tag_propensity = hill_growth_fraction(k = k_tl_10, c = self.cell_count, c_max = c_max, x = s, 
        a = a_tl_10, positive = True, n = n)

        r2 = Reaction([], [self.peptide_chain], propensity_type = non_deg_tag_propensity)
        rxn.append(r2)

        # dilution_propensity = decrease_linear_growth_fraction(k = k_gr, c = self.cell_count, c_max = c_max,
        # x = s, a = r_dil)

        # r3 = Reaction([], [self.aminoacid], propensity_type = dilution_propensity)
        # rxn.append(r3)

        # r4 = Reaction.from_massaction([self.peptide_chain], [self.aminoacid],
        #                               k_forward = k_tl_11)
        # rxn.append(r4)

        return rxn




# This class contains all the non-deg tag protein degradation related reactions for the CRN 
class BacterialDegradation_active(GlobalMechanism):
    """Non deg-tag degradation of proetin or the endogenous protein degradation 
    This involves not just the protease but also other forms of protein degradation like bleaching 
       All species with the attribute degtagged and material_type protein are degraded. The method is not recursive.
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

            k_tl_10 = self.get_parameter(s, "k_tl_10", mixture)
            a_tl_10 = self.get_parameter(s, "a_tl_10", mixture)
            n = self.get_parameter(s, "n", mixture)
            c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
            rxns = []
            
        #     propensity_dp = hill_growth_fraction(k = k_tl_10, c = self.cell_count, c_max = c_max, x = s, 
        # a = a_tl_10, b = b_tl_10, positive = positive_tl_10, n = n)
            propensity_dp = hill_growth(k = k_tl_10, c = self.cell_count, c_max = c_max, x = s, 
        positive = True, n = n)
            
            r1 = Reaction([s], [self.peptide_chain], propensity_type = propensity_dp)
            rxns.append(r1)
            
        return rxns


# This class contains all the deg-tag protein degradation related reactions for the CRN 
class BacterialDegradationTagged(GlobalMechanism):
    """Deg-tag degradation of proetin or the ssrA tag protein degradation 
   
       All species with the attribute degtagged and material_type protein are degraded. 
    """
    def __init__(self,aminoacid, peptide_chain, protease, name="deg_tag_degradation", 
                 cell_count = None, mechanism_type = "deg_tag_degradation", 
                 filter_dict= None, recursive_species_filtering = False, default_on = False, **keywords):
        
        if isinstance(aminoacid, Species):
            self.aminoacid = aminoacid
        else:
            raise ValueError("'aminoacid' must be a Species.")
        
        if isinstance(peptide_chain, Species):
            self.peptide_chain = peptide_chain
        else:
            raise ValueError("'peptide_chain' must be a Species.")
        
        if isinstance(protease, Species):
            self.protease = protease
        else:
            raise ValueError("'protease' must be a Species.")
        
        if isinstance(deg_complex, Species):
            self.deg_complex = deg_complex
        else:
            raise ValueError("'deg_complex' must be a Species.")
            
        if isinstance(cell_count, Species):
            self.cell_count = cell_count
        else:
            raise ValueError("'cell_count' parameter must be a Species.")

        
            
        GlobalMechanism.__init__(self, name = name, mechanism_type = mechanism_type, default_on = default_on,
                                 filter_dict = filter_dict, recursive_species_filtering = recursive_species_filtering)

    def update_species(self, s, mixture):
        species = [s, self.protease, self.aminoacid, self.deg_complex]
        return species

    def update_reactions(self, s, mixture):
        if s.material_type == "protein": 

            rxns = []

            k_tl_9 = self.get_parameter(s, "k_tl_9", mixture)
            k_tl_13b = self.get_parameter(s, "k_tl_13b", mixture)
            k_tl_13u = self.get_parameter(s, "k_tl_13u", mixture)

            
            n = self.get_parameter(s, "n", mixture)
            c_max = mixture.get_parameter(param_name = "c_max", part_id = None, mechanism = "logistic_cell_growth")
            
            
            propensity_dp = hill_growth(k = k_tl_9, c = self.cell_count, c_max = c_max, x = self.deg_complex, 
        positive = False, n = n)
            
            r1 = Reaction.from_massaction([s, self.protease], [self.deg_complex], k_forward = k_tl_13b, k_reverse = k_tl_13u)
            rxns.append(r1)

            r2 = Reaction.from_massaction([self.deg_complex], [self.peptide_chain, self.protease], k_forward = k_tl_9)
            #r2 = Reaction([deg_complex], [self.peptide_chain, self.protease], propensity_type = propensity_dp)
            rxns.append(r2)
            
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
            dilution_propensity = decrease_linear_growth(k = k_gr, c = self.cell_count, c_max = c_max, x = s)
            r = Reaction([s], [], propensity_type = dilution_propensity)
            rxn.append(r)
        else:
            dilution_propensity = decrease_linear_growth(k = k_gr, c = self.cell_count, c_max = c_max, x = s)
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