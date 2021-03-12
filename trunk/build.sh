#!/bin/bash
help(){
  echo "Enter comiler: "
  echo "  Options: ifort, gfortran"
}

if [[ $# -lt 1 ]]; then
  help
  exit 1
fi

for var in "$@"
do
  if [ "$var" == "--help" -o "$var" == "-help" -o "$var" == "-h" ]; then
    help
    exit 0

  elif [ "$var" == "ifort" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    cd serial
    make ifort
    cd ..

  elif [ "$var" == "gfortran" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    cd serial
    make gfortran
    cd ..
  fi 
done

echo
echo 
echo ">>>> Executable available in ./executables"

