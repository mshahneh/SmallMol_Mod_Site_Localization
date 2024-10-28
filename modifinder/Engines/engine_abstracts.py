from abc import ABC, abstractmethod
from modifinder.utilities.gnps_types import SpectrumTuple
from typing import List, Tuple


# Base class for alignment engines
class AlignmentEngine(ABC):
    @abstractmethod
    def align(self, network):
        pass
    
    def single_align(self, SpectrumTuple1: SpectrumTuple,
                      SpectrumTuple2: SpectrumTuple, 
                      fragment_mz_tolerance: float = 0.02, 
                      fragment_ppm_tolerance: float = 100.0) -> Tuple[float, List[Tuple[int, int]]]:
        """
        Aligns two spectra using cosine similarity and returns the cosine score and the matched peaks

        Parameters:
            SpectrumTuple1 (SpectrumTuple): First spectrum
            SpectrumTuple2 (SpectrumTuple): Second spectrum
            fragment_mz_tolerance (float): Fragment mz tolerance
            fragment_ppm_tolerance (float): Fragment ppm tolerance
        
        Returns:
            Tuple[float, List[Tuple[int, int]]]: alignment score and the list of matched peaks, each value in the list is a tuple of the indices of the matched peaks in the two spectra
        """
        pass

# Base class for annotation engines
class AnnotationEngine(ABC):
    @abstractmethod
    def annotate(self, network, alignment):
        pass

# Base class for prediction engines
class PredictionEngine(ABC):
    @abstractmethod
    def predict(self, network, annotation, alignment):
        pass

    @abstractmethod
    def calculate_confidence(self, network, annotation, alignment, prediction):
        pass