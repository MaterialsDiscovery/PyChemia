class StructureDistances:

    def __init__(self, pcdb):
        self.pcdb = pcdb

    def get_distance(self, ids_pair):
        return self.pcdb.db.distances.find_one({'pair': ids_pair})

    def set_distance(self, ids_pair, distance):
        self.pcdb.db.distances.insert_one({'pair': ids_pair, 'distance': distance})


class FingerPrints:

    def __init__(self, pcdb):
        self.pcdb = pcdb

    def set_fingerprint(self, fingerprint):
        return self.pcdb.db.fingerprints.insert_one(fingerprint)

    def get_fingerprint(self, entry_id):
        return self.pcdb.db.fingerprints.find_one({'_id': entry_id})

    def update(self, entry_id, fingerprint):
        self.pcdb.db.fingerprints.update({'_id': entry_id}, fingerprint)
