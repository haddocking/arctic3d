"""arctic3d"""
import logging
import sys

log = logging.getLogger(__name__)
log.handlers.clear()
log.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(
    logging.Formatter("[%(asctime)s %(module)s %(levelname)s] %(message)s")
)
log.addHandler(handler)
