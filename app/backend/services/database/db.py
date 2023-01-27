import sqlite3

from .models import SEQUENCE_TABLE, FRAME_TABLE

# |-----------------------------------------------------------------------------|#


class DB:
    def __init__(self):
        self.con = sqlite3.connect("backend/database/database.sqlite")

# |-----------------------------------------------------------------------------|#

    def connect(func):
        def wrapper(self, *args):
            with self.con as cur:
                return func(self, cur, *args)
        return wrapper

# |-----------------------------------------------------------------------------|#

    @connect
    def migrate(self, cur):
        cur.execute(SEQUENCE_TABLE)
        cur.execute(FRAME_TABLE)

        cur.commit()

# |-----------------------------------------------------------------------------|#

    @connect
    def get_sequence(self, cur, _id : int):
        query = """SELECT sequence FROM sequences WHERE id = :id"""
        return cur.execute(query, {"id": _id}).fetchone()[0]

# |-----------------------------------------------------------------------------|#

    @connect
    def get_frame(self, cur, _id : int):
        query = """SELECT sequence FROM frames WHERE id = :id"""
        return cur.execute(query, {"id": _id}).fetchone()[0]

# |-----------------------------------------------------------------------------|#

    @connect
    def post_sequence(self, cur, sequence: str):
        query = """INSERT INTO sequences (sequence) VALUES (:sequence) ON CONFLICT DO NOTHING"""
        cur.execute(query, {"sequence": sequence})

        cur.commit()

        query = """SELECT id FROM sequences WHERE sequence = :sequence"""
        return {"id": cur.execute(query, {"sequence": sequence}).fetchone()[0]}

# |-----------------------------------------------------------------------------|#

    @connect
    def post_frame(self, cur, frame: str):
        query = """INSERT INTO frames (sequence) VALUES (:frame) ON CONFLICT DO NOTHING"""
        cur.execute(query, {"frame": frame})

        cur.commit()

        query = """SELECT id FROM frames WHERE sequence = :frame"""
        return {"id": cur.execute(query, {"frame": frame}).fetchone()[0]}
