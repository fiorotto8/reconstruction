echo "Entering the CYGNO virtual environment software..."
echo "Setting up python3.8 from central installation on /nfs..."
export PYTHONPATH=$PYTHONPATH:/nfs/cygno/users/$USER/local/lib/python3.8/site-packages/
alias python="python3.8"
alias pip="pip3.8"
echo "Setting up ROOT..."
source /nfs/cygno/software/root-v6-22-02-py38-install/bin/thisroot.sh
echo "DONE."
