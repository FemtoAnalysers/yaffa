'''
Init script
'''

import logging

print("Welcome to yaffa!")

# Create logging formatter for beautiful log messages
class CustomFormatter(logging.Formatter):
    '''
    DEBUG: for very detailed output
    INFO: for tracking the progress of the program
    WARNING: for something unexpected that doesn't cause problems
    ERROR: for something unexpected that can cause problems
    CRITICAL: for a failure that cause the program to terminate
    '''
    grey = "\x1b[90;20m"
    blue = "\x1b[34;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    purple = "\x1b[35;1m"
    reset = "\x1b[0m"
    style = "%(module)s::%(funcName)s %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    criticalStyle = "%(module)s::%(funcName)s %(levelname)s - %(message)s ---> EXIT! (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + style + reset,
        logging.INFO: blue + style + reset,
        logging.WARNING: yellow + style + reset,
        logging.ERROR: red + style + reset,
        logging.CRITICAL: purple + criticalStyle + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

# Create handler to exit the program for CRITICAL errors
class ExitOnExceptionHandler(logging.StreamHandler):
    '''Exit on CRITICAL'''

    def emit(self, record):
        super().emit(record)
        if record.levelno is logging.CRITICAL:
            raise SystemExit(-1)

# create logger
logger = logging.getLogger("testprj")
logger.setLevel(logging.INFO)
ch = ExitOnExceptionHandler()
ch.setFormatter(CustomFormatter())
logger.addHandler(ch)
