from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError


def has_local_mongo():
    try:
        max_server_selection_delay = 1
        client = MongoClient("localhost", serverSelectionTimeoutMS=max_server_selection_delay)
        client.server_info()  # force connection on a request as the
        return True
    except ServerSelectionTimeoutError as err:
        print(err)
        return False

