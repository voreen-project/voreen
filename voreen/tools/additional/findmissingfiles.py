import os

dirsToIgnore = ['.svn', 'debug', 'release', 'Debug', 'Release', 'ext']

# tests, if the file 'file' occurs in any line of the project file 'proFile'
# If so, the line is returned
# otherwise, None is returned
def isInProFile(proFile, file):
    f = open(proFile, 'r')
    lines = f.readlines()

    file = file.replace('\\', '/')
    file = file.lstrip()

    for line in lines:
        line = line.lstrip()
        line = line[:len(file)]

        if (line == file):
            return file

# Walks the 'baseCPP' and 'baseH' directories and tests all files if they
# are present in the project file 'proFile'
# prints all files to the stdout which have no equivalent in the project file
def pantera(baseCPP, baseH, proFile):
    print proFile + ' ::: cpp'
    print '----------'
    for root, dirs, files in os.walk(baseCPP):
        for dir in dirsToIgnore:
            try:
                dirs.remove(dir)
            except ValueError:
                pass
        
        for file in files:
            if (file.find('.cpp') != -1):
                complete = os.path.join(root, file)[len(baseCPP):]
                if (isInProFile(baseCPP + proFile, complete) == None):
                   print os.path.join(root, file)

    print '\n'
    print proFile + ' ::: h'
    print '----------'

    for root, dirs, files in os.walk(baseH):
        for dir in dirsToIgnore:
            try:
                dirs.remove(dir)
            except ValueError:
                pass
        
        for file in files:
            if (file.find('.h') != -1):
                complete = os.path.join(root, file)[len(baseH):]
                if (isInProFile(baseCPP + proFile, '../' + baseH + complete) == None):
                    print os.path.join(root, file)

    
pantera('..\\src\\core\\', '..\\include\\voreen\\core\\', 'voreenlib_core.pro')
print '\n'
pantera('..\\src\\qt\\', '..\\include\\voreen\qt\\', 'voreenlib_qt.pro')
