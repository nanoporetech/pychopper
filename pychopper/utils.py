
from collections import OrderedDict

def parse_config_string(s):
    res = OrderedDict()
    s = s.strip()
    for token in s.split("|"):
        token = token.strip().split(":")
        if len(token) != 2:
            raise Exception("Invalid token: " + token)
        strand = token[0].strip()
        if strand not in ('+', '-'):
            raise Exception("Invalid strand: " + strand)
        token = token[1].strip().split(",")
        res[(token[0].strip(), token[1].strip())] = strand
    return res

