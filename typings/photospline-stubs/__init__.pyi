from typing import Annotated, Sequence, Union, Any, overload
from os import PathLike
import numpy as np
from scipy import sparse
import numpy.typing as npt

Float1D = npt.ArrayLike[tuple[int], np.floating[Any]]
Int1D = npt.ArrayLike[tuple[int], np.integer[Any]]

class SplineTable:
    coefficients: npt.NDArray[np.float32]
    extents: tuple[tuple[float, float], ...]
    knots: tuple[npt.NDArray[np.float32], ...]
    ndim: int
    order: tuple[int, ...]
    def aux_value(self, key: str) -> str: ...
    def __init__(self, path: PathLike): ...
    @classmethod
    def stack(
        cls,
        tables: Sequence[SplineTable],
        coordinates: Float1D,
        stackOrder: int = 2,
    ) -> SplineTable:
        """ """
        ...
    def convolve(self, dim: int, knots: Float1D) -> None: ...
    def permute_dimensions(self, permutation: Sequence[int]) -> None: ...
    def write(self, path: str) -> None: ...
    @overload
    def search_centers(self, x: Sequence[float]) -> tuple[int]: ...
    @overload
    def search_centers(
        self, x: npt.NDArray[np.floating[Any]]
    ) -> npt.NDArray[np.longlong]: ...
    def search_centers(
        self, x: Union[npt.NDArray[np.floating[Any]], Sequence[float]]
    ) -> Union[npt.NDArray[np.longlong], tuple[int]]: ...
    def deriv(
        self, x: Float1D, centers: Int1D, derivatives: Int1D
    ) -> float: ...
    @overload
    def evaluate(
        self, x: Float1D, centers: Int1D, derivatives: Int1D = 0
    ) -> float: ...
    @overload
    def evaluate(
        self,
        x: npt.NDArray[np.floating[Any]],
        centers: npt.NDArray[np.longlong],
        derivatives: int = 0,
    ) -> npt.NDArray[np.floating[Any]]: ...
    def evaluate(
        self,
        x: Union[npt.NDArray[np.floating[Any]], Sequence[float]],
        centers: Union[npt.NDArray[np.longlong], Sequence[int]],
        derivatives: int = 0,
    ) -> Union[npt.NDArray[np.floating[Any]], float]: ...
    def evaluate_gradient(
        self, x: Float1D, centers: Int1D
    ) -> float: ...
    def evaluate_simple(self, x: Float1D) -> float: ...
    def __call__(self, x: Float1D) -> float: ...
    def grideval(self, coords: Sequence[Float1D]) -> npt.NDArray[np.float64]: ...

class ndsparse:
    def __init__(self, rows: int, ndim: int) -> None: ...
    @classmethod
    def from_data(
        cls, values: npt.ArrayLike, weights: None | npt.ArrayLike = None
    ) -> tuple[ndsparse, ndsparse]: ...
    def insert(self, value: float, indices: Sequence[int]) -> None: ...

def bspline(knots: npt.ArrayLike, x: float, index: int, order: int) -> float: ...
def glam_fit(
    data: ndsparse,
    weights: ndsparse,
    coordinates: Sequence[Float1D],
    knots: Sequence[Float1D],
    order: Int1D,
    smoothing: Sequence[float],
    penaltyOrder: Sequence[int],
    monodim: None | int = None,
    verbose: int = 1,
) -> SplineTable: ...

def nnls(
    A: sparse.csc_matrix,
    b: npt.NDArray[np.float64],
    tolerance: float=0.,
    min_iterations: int=0,
    max_iterations: int=np.iinfo(np.int32).max,
    npos: int=0,
    normaleq: bool=False,
    verbose: bool=False,
) -> npt.NDArray[np.float64]: ...
