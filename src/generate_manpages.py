'''
(re)Generate Phyx manpages.

Requires:
1) help2man
2) all programs to be compiled (i.e., present in the src dir)
manpages will be written to man/px*.1.in

Make sure to run this script from the src dir.

The idea here is that developers will make/regenerate these pages.
When a user compiles a program, the manpage px*.1.in is copied to px*.1, and
it is this copy that is installed with make install. (Likewise, when a user
does make clean, man/px*.1 files are remove, and man/px*.1.in remain)

The rationale behind this is 1) users need not worry about having help2man,
and 2) manpages are installed only for those programs that are actually compiled.
This makes even more sense if we want to idiosyncratically edit the px*.1.in
files in some way that cannot be automated.
'''

import sys
import os
import subprocess

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

# feel free to edit the command below to add extras (e.g., citations for specific methods)
def make_manpage(pname):
    cmd = "help2man --output=man/" + pname + ".1.in --include=man/citation.h2m --no-info ./" + pname
    print ("Writing manpage for " + pname)
    os.system(cmd)
    return None

if __name__ == "__main__":
    if len(sys.argv) != 1:
        print ("python run_tests.py")
        sys.exit(0)
    
    # check that help2man is installed
    if which("help2man") is None:
        print ("Error: help2man must be installed to generate manpages. Exiting.")
        sys.exit()
       
    n = 0
    print ("=================")
    for i in sorted(os.listdir(".")):
        if i[:2] == "px":
            make_manpage(i)
            print ("=================")
            n +=1
    print ("\nWrote " + str(n) + " manpages to /man")
    print ("Fin.")
