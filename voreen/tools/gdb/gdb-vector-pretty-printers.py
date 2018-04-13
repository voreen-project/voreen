import gdb
import re

class TgtVec2Printer(object):
    "Print a tgt::Vector2"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return "{}({}, {})".format(self.val.type.name, self.val['x'], self.val['y'])

class TgtVec3Printer(object):
    "Print a tgt::Vector3"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return "{}({}, {}, {})".format(self.val.type.name, self.val['x'], self.val['y'], self.val['z'])

class TgtVec4Printer(object):
    "Print a tgt::Vector4"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return "{}({}, {}, {}, {})".format(self.val.type.name, self.val['x'], self.val['y'], self.val['z'], self.val['w'])

def tgtvec_lookup_function(val):
    lookup_tag = val.type.name
    if lookup_tag is None:
        return None

    regex_vec2 = re.compile("^tgt::.?vec2$")
    regex_vec2_long = re.compile("^tgt::Vector2<.*>$")
    regex_vec3 = re.compile("^tgt::.?vec3$")
    regex_vec3_long = re.compile("^tgt::Vector3<.*>$")
    regex_vec4 = re.compile("^tgt::.?vec4$")
    regex_vec4_long = re.compile("^tgt::Vector4<.*>$")
    if regex_vec2.match(lookup_tag) or regex_vec2_long.match(lookup_tag):
        return TgtVec2Printer(val)
    if regex_vec3.match(lookup_tag) or regex_vec3_long.match(lookup_tag):
        return TgtVec3Printer(val)
    if regex_vec4.match(lookup_tag) or regex_vec4_long.match(lookup_tag):
        return TgtVec4Printer(val)

    return None

gdb.pretty_printers.append(tgtvec_lookup_function)
