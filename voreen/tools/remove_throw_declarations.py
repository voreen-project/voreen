#! /usr/bin/python3

# This script is intended to remove all throw (<exception>) function declarations
# since they are currently deprecated and whose support will be dropped soon.
# Can be executed by the following command (example for .cpp files):
#   find <voreen_dir> -type f -name "*.cpp" -exec python tools/remove_throw_declarations.py {} \;

# Usage: All parameters are interpreted as paths to files that are
# supposed to be processed by this script.

import sys
import re

warningMsg = "TODO: THROW SHOULD NOT USE DYNAMIC OBJECTS"

def re_sub(pattern, replacement, string):
    def _r(m):
        # Now this is ugly.
        # Python has a "feature" where unmatched groups return None
        # then re.sub chokes on this.
        # see http://bugs.python.org/issue1519638
        
        # this works around and hooks into the internal of the re module...

        # the match object is replaced with a wrapper that
        # returns "" instead of None for unmatched groups

        class _m():
            def __init__(self, m):
                self.m=m
                self.string=m.string
            def group(self, n):
                return m.group(n) or ""

        return re._expand(pattern, _m(m), replacement)
    return re.sub(pattern, _r, string)

# Remove throw declaration
replacements = [
        [r"\)\s*?(\s?const)?\s*throw\s*\([^\)].*\)\s*?(\{|(\n\s*)?:)", r")\1 \2"], # definition
        [r"\)\s*?(\s?const)?\s*throw\s*\([^\)].*\)\s*?(\s?(=\s*0|override)?;)", r")\1\2"], # declaration        
        [r"(throw\s*?\(?\s*new)", r"/*" + warningMsg + r"*/ \1"], # dynamic objects
        #[r"(throw\s*?\(.*?\))", r"/*\1*/"] # uncomment
    ]

def needsProcessing(content):
    for replacement in replacements:
        if re.search(replacement[0], content):
            return True

    return False

def handleReplacement(content):
    for replacement in replacements:
        content = re_sub(replacement[0], replacement[1], content)
    return content

def processFile(filename):
    data = ""
    with open(filename, 'rb') as f:
        try:
            data = f.read().decode('utf8')
        except:
            print("WARNING: Cannot decode {} as utf8, trying iso-8859-1".format(filename))
            data = f.read().decode('iso-8859-1')

    if needsProcessing(data):
        print("Remove throw declaration in {} ".format(filename))
        updatedDate = handleReplacement(data)
        with open(filename, 'w') as f:
            f.write(updatedDate.encode('utf8'))
    else:
        print("{} does not contain throw declarations".format(filename))

def main():
    for filename in sys.argv[1:]:
        print("Processing file {}".format(filename))
        processFile(filename)
        print("Done: {}\n".format(filename))

if __name__ == "__main__":
    main()
