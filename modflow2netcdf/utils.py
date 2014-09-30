import time

class LoggingTimer(object):
    def __init__(self, message, logger_func=None):
        self.message     = message
        self.logger_func = logger_func

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.logger_func:
            self.logger_func("[%s seconds] - %s" % (round(self.secs, 4), self.message))
