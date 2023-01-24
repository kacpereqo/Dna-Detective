from databases import Database

# |-----------------------------------------------------------------------------|#


class DB:
    def __init__(self):
        self.db = Database("sqlite://app/backend/database/database.sqlite")

# |-----------------------------------------------------------------------------|#

    async def connect(self):
        await self.db.connect()

# |-----------------------------------------------------------------------------|#

    async def disconnect(self):
        await self.db.disconnect()

# |-----------------------------------------------------------------------------|#

    async def migrate(self):
        await self.db.execute("""CREATE TABLE IF NOT EXISTS sequences (id INTEGER PRIMARY KEY AUTOINCREMENT, sequence TEXT NOT NULL UNIQUE)""")

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
