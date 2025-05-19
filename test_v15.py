from typing import overload, Callable, Union
import numbers


class Signal:
    def __init__(self, fn: Callable[[float], float]) -> None:
        if not callable(fn):
            raise TypeError("fn 必須是 Callable")
        self._fn = fn

    @classmethod
    def from_scalar(cls, value: float) -> "Signal":
        return cls(lambda _t, v=value: v)
sig1 = Signal(lambda t: t**2)
sig2 = Signal.from_scalar(42)
