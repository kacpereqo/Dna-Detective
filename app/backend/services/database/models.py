SEQUENCE_TABLE = """CREATE TABLE IF NOT EXISTS sequences (id INTEGER PRIMARY KEY AUTOINCREMENT, sequence TEXT NOT NULL UNIQUE)"""

FRAME_TABLE = """CREATE TABLE IF NOT EXISTS frames (id INTEGER PRIMARY KEY AUTOINCREMENT, sequence TEXT NOT NULL UNIQUE)"""
