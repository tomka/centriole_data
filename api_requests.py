import requests
from requests.auth import HTTPBasicAuth

import config


class CatmaidApiTokenAuth(HTTPBasicAuth):
    """Attaches HTTP X-Authorization Token headers to the given Request.
    Optionally, Basic HTTP Authentication can be used in parallel.
    """

    def __init__(self, token, username=None, password=None):
        super(CatmaidApiTokenAuth, self).__init__(username, password)
        self.token = token

    def __call__(self, r):
        r.headers["X-Authorization"] = "Token {}".format(self.token)
        if self.username and self.password:
            super(CatmaidApiTokenAuth, self).__call__(r)
        return r


AUTH_TOKEN = CatmaidApiTokenAuth(config.token, config.username, config.password)
