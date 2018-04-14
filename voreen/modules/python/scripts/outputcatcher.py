# internal helper script used to redirect sys.stdout/sys.stderr
# to the C++ PythonModule.
import sys
import voreen_internal

class OutputCatcherStdout:
    def write(self, string):
        voreen_internal.vrnPrint(string, 0)

class OutputCatcherStderr:
    def write(self, string):
        voreen_internal.vrnPrint(string, 1)
        
sys.stdout = OutputCatcherStdout()  
sys.stderr = OutputCatcherStderr()  
