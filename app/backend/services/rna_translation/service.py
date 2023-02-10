from .constants import RNA_CODON_TABLE, COMPLEMENTARY_TABLE
from typing import List
import re
from backend.services.database.db import DB

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

        return translated_frame

    # |------------------------------------------------------------------------------|#

    def find_open_reading_frames(self, translated_sequence: str) -> List[str]:
        seq = ""
        index = 0

        for protein in translated_sequence.split("-"):
            if len(protein) > 0:
                if protein[0] != "M":
                    if protein[0] not in self.ids:
                        self.ids[protein[0]] = self.db.post_frame(protein)['id']

                    protein = f"""<a href="#/analize/{self.ids[protein[0]]}" class = "frame">{protein[0]}</a>""" + protein[1:]
                    seq += protein
                else:

                    seq += re.sub(
                        r"M[A-Z]+", lambda x: f"""<a href="#/analize/{self.db.post_frame(x.group(0))['id'] if x.group(0) not in self.ids else self.ids[x.group(0)]}" class="frame">{x.group(0)}</a>""", protein)

            seq += "-"
            index += len(protein) + 1

        return seq[:-1]


# |------------------------------------------------------------------------------|#


    def parse(self):
        frames = self.get_frames(self.is_reversed, self.is_forward)
        result = {}

        for i, (direction, _frames) in enumerate(frames.items()):
            for j, frame in enumerate(_frames):
                translated_frame = self.translate_frames(frame)
                open_reading_frames = self.find_open_reading_frames(
                    translated_frame)

                result[i * 3 + j] = {
                    "frame": frame,
                    "shift": j,
                    "direction": ["5'3", "3'5"] if direction == "5'3'" else ["3'5", "5'3"],
                    "translatedFrame": translated_frame,
                    "openReadingFrames": open_reading_frames,
                }

        return result

# |------------------------------------------------------------------------------|#
