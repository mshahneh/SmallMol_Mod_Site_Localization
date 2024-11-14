from abc import ABC, abstractmethod
from modifinder.classes.Spectrum import Spectrum
from modifinder.classes.Compound import Compound
from modifinder.classes.ModiFinder import ModiFinder
from modifinder.classes.EdgeDetail import EdgeDetail
from typing import List, Tuple


# Base class for alignment engines
class AlignmentEngine(ABC):
    @abstractmethod
    def align(self, network: ModiFinder, **kwargs):
        pass
    
    @abstractmethod
    def single_align(self, SpectrumTuple1: Spectrum,
                      SpectrumTuple2: Spectrum, 
                      fragment_mz_tolerance: float = 0.02, 
                      fragment_ppm_tolerance: float = 100.0, 
                      **kwargs) -> EdgeDetail:
        """
        Aligns two spectra using cosine similarity and returns the cosine score and the matched peaks

        Parameters:
            :SpectrumTuple1 (SpectrumTuple): First spectrum
            :SpectrumTuple2 (SpectrumTuple): Second spectrum
            :fragment_mz_tolerance (float): Fragment mz tolerance
            :fragment_ppm_tolerance (float): Fragment ppm tolerance
            :kwargs: additional arguments
        
        Returns:
            :EdgeDetail: the edge detail object
        """
        pass

# Base class for annotation engines
class AnnotationEngine(ABC):
    @abstractmethod
    def annotate(self, network, **kwargs):
        pass

    
    @abstractmethod
    def annotate_single(self, compound: Compound,  modify_compound: bool = True, **kwargs) -> List[Tuple[int, List[str]]]:
        """
        provides annotation for the peaks in a single compound

        Parameters:
            :compound (Compound): the compound to be annotated
            :modify_compound (bool): whether to modify the passed compound with the annotations
            :kwargs: additional arguments
        
        Returns:
            :Mapping[int, List[str]]: a dictionary with the indices of the peaks as keys and the list of annotations as values
        """
        pass


    @abstractmethod
    def get_fragment_info(self, Compound: Compound, fragment):
        """
        converts a fragment to a SMILES string

        Parameters:
            :Compound (Compound): the compound
            :fragment (int): the fragment
        
        Returns:
            :tuple(atomlist -> List[int], edge_list -> List[Tuple[int, int]], formula -> str, smiles -> str): a tuple containing the atom list, the edge list, the formula and the SMILES string
        """
        pass

# Base class for prediction engines
class PredictionEngine(ABC):
    @abstractmethod
    def predict(self, network: ModiFinder, **kwargs):
        pass

    @abstractmethod
    def confidence(self, network:ModiFinder, prediction, **kwargs):
        pass