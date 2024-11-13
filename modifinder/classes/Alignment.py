from dataclasses import dataclass, field
from typing import List, Tuple

@dataclass
class Alignment:
    matches: List[Tuple[int, int]] = field(default_factory=list) # List of matched peaks
    score: float = 0.0                               # Float for score
    shifts: List[int] = field(default_factory=list)   # index in matches that are shifted
    unshifted: List[int] = field(default_factory=list) # index in matches that are not shifted