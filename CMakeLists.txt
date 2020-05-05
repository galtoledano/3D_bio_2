cmake_minimum_required(VERSION 3.14)
project(3D_bio_2)

set(CMAKE_CXX_STANDARD 14)

add_executable(3D_bio_2 structalign.cc structalign.h Vector3.h Atom.h RigidTrans3.h Matrix3.h Molecule.h PDB.h Match.h GeomHash.h Triangle.h Vector3.cc Atom.cc RigidTrans3.cc Matrix3.cc PDB.cc Match.cc Triangle.cc)