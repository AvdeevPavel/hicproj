import re
from typing import List, Tuple


# TODO add soft and hard clipping
# TODO create unit tests
class Cigar(object):
    _SUPPORTED_LETTERS: str = "M=IDX"  # SNHP]

    def __init__(self, inp_str: str = ''):
        self._str: str = inp_str
        self._events: List[Tuple[int, str]] = [(int(x), l) for x, l in
                                               re.findall('(\d+)([' + self._SUPPORTED_LETTERS + '])', inp_str)]

    def __str__(self):
        return self._str

    def _complement(self, c: str) -> str:
        if c == 'I':
            return 'D'
        elif c == 'D':
            return 'I'
        elif c in self._SUPPORTED_LETTERS:
            return c
        else:
            print(c)
            assert False

    def size_on_ref(self) -> int:
        sm: int = 0
        for n, l in self._events:
            if l != 'I':
                sm += n
        return sm

    def size_on_query(self) -> int:
        sm: int = 0
        for n, l in self._events:
            if l != 'D':
                sm += n
        return sm

    def complement(self) -> "Cigar":
        """
        Reverse-complement cigar: Switches I/D and reverses the order
        :return:
        """
        out: str = ''
        for n, l in self._events[::-1]:
            out += '{0}{1}'.format(n, self._complement(l))
        return Cigar(out)
