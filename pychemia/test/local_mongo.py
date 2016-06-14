import pymongo


def has_local_mongo():
    try:
        maxSevSelDelay = 1
        client = pymongo.MongoClient("localhost", serverSelectionTimeoutMS=maxSevSelDelay)
        client.server_info()  # force connection on a request as the
        return True
    except pymongo.errors.ServerSelectionTimeoutError as err:
        print(err)
        return False
