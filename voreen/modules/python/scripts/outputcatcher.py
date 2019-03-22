# internal helper script used to redirect sys.stdout/sys.stderr
# to the C++ PythonModule.
import sys
import voreen_internal

class OutputCatcherStdout:
    def write(self, string):
        voreen_internal.vrnPrint(string, 0, __voreen_script_id__)
    def flush(self):
        pass

class OutputCatcherStderr:
    def write(self, string):
        voreen_internal.vrnPrint(string, 1, __voreen_script_id__)
    def flush(self):
        pass
        
sys.stdout = OutputCatcherStdout()  
sys.stderr = OutputCatcherStderr()  
