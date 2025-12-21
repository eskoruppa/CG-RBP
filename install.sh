# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# the CG-DNA package cannot be used without
# the MOLECULE and ASPHERE packages installed.

if (test $1 = 1) then
  if (test ! -e ../bond_rbp.cpp) then
    echo "Must install MOLECULE package with CG-RBP"
    exit 1
  fi
  if (test ! -e ../fix_nve_asphere.cpp) then
    echo "Must install ASPHERE package with CG-RBP"
    exit 1
  fi
fi


# list of files with dependcies

action  bond_rbp.cpp
action  bond_rbp.h
action  bond_rbp_fene.cpp
action  bond_rbp_fene.h
action  so3.h
action  lamath.h
action  angle_rbp.cpp
action  angle_rbp.h
action  so3.h
action  so3.h
