from dataclasses import dataclass

@dataclass(frozen=True)
class VariableSpec:
    name: str
    dims: tuple
    timeDependent: bool