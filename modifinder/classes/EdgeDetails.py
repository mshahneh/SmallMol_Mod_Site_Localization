from dataclasses import dataclass, field
from typing import List, Tuple
from modifinder.classes.Alignment import Alignment

@dataclass
class EdgeDetails:
    number_of_modifications: int = None
    # make a default alignment object for the edge
    alingment: Alignment = None
    is_smaller: bool = False