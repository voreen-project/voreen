import os
import string

dirsToIgnore = ['.svn', 'debug', 'release', 'Debug', 'Release', 'ext']

def removeTabsFromFile(file):
    f = open(file)
    lines = f.readlines()
    f.close()

    f = open(file, 'w')
    for line in lines:
        f.write(line.expandtabs(4))
    f.close()

def removeTrailingWhiteSpacesFromFile(file):
    f = open(file)
    lines = f.readlines()
    f.close()

    f = open(file, 'w')
    for line in lines:
        a = line.rstrip()
        f.write(a + '\n')
    f.close()

def pantera(directory):
    for root, dirs, files in os.walk(directory):
        for dir in dirsToIgnore:
            try:
                dirs.remove(dir)
            except ValueError:
                pass
        for file in files:
            if (file.endswith('.cpp') or file.endswith('.h') or file.endswith('.vert') or file.endswith('.frag') or file.endswith('.pro')):
                removeTabsFromFile(root + '\\' + file)
                removeTrailingWhiteSpacesFromFile(root + '\\' + file)

pantera('..\\..\\src')
pantera('..\\..\\include\\voreen')
