import os, sys

def check_if_directory_exists(mydir):
    if os.path.isdir(mydir):
        sys.exit("Directory: \"" + mydir + "\" already exists. You must remove/rename it to run. Exiting.")
    else:
        print("Making \"" + mydir + "\". This is where all outputs will be stored.")
        os.makedirs(mydir)
        return
