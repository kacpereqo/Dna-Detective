import sqlite3
from pymongo import MongoClient
from datetime import datetime
from bson.objectid import ObjectId
from .models import SEQUENCE_TABLE, FRAME_TABLE
from fastapi import HTTPException, status

# |-----------------------------------------------------------------------------|#


class DB:
    def __init__(self):
        DB = "mongodb+srv://kacperek:NymgvVKCliOlC6T9@cluster0.omn7jln.mongodb.net/test"
        self.client = MongoClient(DB)
        self.db = self.client["dna-detective"]

# |-----------------------------------------------------------------------------|#

    def migrate(self):
        try:
            self.db.frames.create_index(
                [('frame', 1)], unique=True)
        except Exception as e:
            pass

        try:
            self.db.sequences.create_index(
                [('sequence', 1)], unique=True)
        except Exception as e:
            pass

        try:
            self.db.users.create_index(
                [('email', 1)], unique=True)
        except Exception as e:
            pass

        # |-----------------------------------------------------------------------------|#

    def get_sequence(self, _id : str):
        """ Get sequence from database by ID"""
        return self.db.sequences.find_one({"_id": ObjectId(_id)})

        # |-----------------------------------------------------------------------------|#

    def get_frame(self, _id : str):
        """ Get frame from database by ID"""
        return self.db.frames.find_one({"_id": ObjectId(_id)})["frame"]

# |-----------------------------------------------------------------------------|#

    def post_sequence(self, sequence: str) -> str:
        """ Post sequence to database and return ID"""
        try:
            _id = self.db.sequences.insert_one(
                {"sequence": sequence}).inserted_id
            return str(_id)
        except Exception as e:
            _id = self.db.sequences.find_one({"sequence": sequence})["_id"]
            return str(_id)

# |-----------------------------------------------------------------------------|#

    def post_frame(self, frame):
        """ Post frame to database and return ID"""
        try:
            _id = self.db.frames.insert_one(
                {"frame": frame}).inserted_id
            return str(_id)
        except Exception as e:
            _id = self.db.frames.find_one({"frame": frame})["_id"]
            return str(_id)

# |-----------------------------------------------------------------------------|#

    def insert_translation(self, sequence, translation):
        """ Add translation to frame in database"""
        try:
            self.db.sequences.update_one(
                {"sequence": sequence}, {"$set": {"translation": translation}})
        except Exception as e:
            pass

# |-----------------------------------------------------------------------------|#

    def get_user(self, email : str) -> str:
        """ Get user from database by e-mail"""
        return self.db.users.find_one({"email": email})

# |-----------------------------------------------------------------------------|#

    def register_user(self, email : str, password : str) -> str:
        """ Register user in database"""
        try:
            self.db.users.insert_one(
                {"email": email, "hashed_password": password})
            return status.HTTP_200_OK
        except Exception as e:
            raise HTTPException(status_code=400, detail="E-mail already used")

# |-----------------------------------------------------------------------------|#

    def add_sequence_to_user(self, _id : str, user_id):
        """ Add sequence to user in database"""
        try:
            self.db.users.update_one(
                {"_id": ObjectId(user_id)}, {"$addToSet": {"sequences": _id}})
        except Exception as e:
            pass

# |-----------------------------------------------------------------------------|#

    def get_user_sequences(self, user_id):
        """ Get user sequences from database"""
        try:
            sequences = self.db.users.find_one(
                {"_id": ObjectId(user_id)})["sequences"]
            return sequences
        except Exception as e:
            return []
