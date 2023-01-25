import sqlite3

from .models import SEQUENCE_TABLE

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


# |-----------------------------------------------------------------------------|#


    async def get_sequence(self, _id : int):
        query = """SELECT id FROM sequences WHERE id = :id"""
        return await self.db.fetch_one(query=query, values={"id": _id})

# |-----------------------------------------------------------------------------|#

    async def post_sequence(self, sequence: str):

        query = """INSERT INTO sequences (sequence) VALUES (:sequence) ON CONFLICT DO NOTHING"""
        await self.db.execute(query=query, values={"sequence": sequence})

        query = """SELECT id FROM sequences WHERE sequence = :sequence"""
        return await self.db.fetch_one(query=query, values={"sequence": sequence})
