from .constants import RNA_CODON_TABLE, COMPLEMENTARY_TABLE
from typing import List
import re

# |------------------------------------------------------------------------------|#


class Translator():
    def __init__(self, rna: str, is_reversed: bool = False, is_forward: bool = True):
        self.rna = rna
        self.frames = self.get_frames(is_reversed, is_forward)
        self.translated_frames = self.translate_frames()
        self.open_reading_frames = self.find_open_reading_frames()

    # |------------------------------------------------------------------------------|#

    def get_frames(self, is_reversed: bool, is_forward: bool) -> List[str]:
        _frames = {}

        if is_forward:
            _frames["5'3'"] = [self.rna[i:]for i in range(3)]

        if is_reversed:
            reversed_frame = "".join(
                list(map(lambda x: COMPLEMENTARY_TABLE[x[::-1]], self.rna))[::-1])

            _frames["3'5'"] = [reversed_frame[i:] for i in range(3)]

        return _frames

    # |------------------------------------------------------------------------------|#

    def translate_frames(self):
        _transladed_frames = {}

        for direction, frames in self.frames.items():
            _transladed_frames[direction] = []
            for frame in frames:
                open_frame = ""
                for i in range(0, len(frame), 3):
                    codon = frame[i:i + 3]
                    if codon in RNA_CODON_TABLE:
                        open_frame += RNA_CODON_TABLE[codon]
                    else:
                        break

                _transladed_frames[direction].append(open_frame)
        return _transladed_frames

    # |------------------------------------------------------------------------------|#

    def find_open_reading_frames(self):
        _open_frames = {}

        for direction, frames in self.translated_frames.items():
            _open_frames[direction] = {}
            for i, frame in enumerate(frames):
                _open_frames[direction][i] = []
                for protein in frame.split("-"):
                    if len(protein) > 0:
                        if protein[0] != "M":
                            _open_frames[direction][i].append(protein[0])
                        _open_frames[direction][i].extend(
                            re.findall(r"M[A-Z]+", protein))

        print(_open_frames)
        return _open_frames


# |------------------------------------------------------------------------------|#
