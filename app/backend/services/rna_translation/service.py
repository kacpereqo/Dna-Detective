from .constants import RNA_CODON_TABLE, COMPLEMENTARY_TABLE
from typing import List
import re
from services.database.db import DB

# |------------------------------------------------------------------------------|#


class Translator():
    def __init__(self, rna: str, is_reversed: bool = False, is_forward: bool = True):
        self.rna = rna
        self.is_reversed = is_reversed
        self.is_forward = is_forward
        self.ids = {}
        self.db = DB()

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

    def translate_frames(self, sequence: str) -> str:

        translated_frame = ""
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            if codon in RNA_CODON_TABLE:
                translated_frame += RNA_CODON_TABLE[codon]
            else:
                break

        return translated_frame.split("-")


# |------------------------------------------------------------------------------|#

    def parse(self):

        frames = self.get_frames(self.is_reversed, self.is_forward)
        result = {}

        for i, (direction, _frames) in enumerate(frames.items()):
            for j, frame in enumerate(_frames):
                translated_frames = self.translate_frames(frame)

                result[str(i * 3 + j)] = {
                    "frame": frame,
                    "shift": j,
                    "direction": ["5'3", "3'5"] if direction == "5'3'" else ["3'5", "5'3"],
                    "translatedFrames": translated_frames,
                }
        DB().insert_translation(self.rna , result)

        return result

# |------------------------------------------------------------------------------|#
