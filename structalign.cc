//
// Created by galtoledano and yahel.yed on 05/05/2020.
//

#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include "Triangle.h"
#include <chrono>
#include <iostream>

int main(int argc , char* argv[]){

  if(argc !=4) {
    std::cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << std::endl;
    exit(1);
  }

  //********Parameters********************
  float epsilon = atof(argv[1]); // distance threshold on atoms in correspondence

  // read the two files into Molecule
  Molecule<Atom> molModel, molTarget, allModelAtoms;

  std::ifstream fileModel(argv[3]);
  std::ifstream fileTarget(argv[2]);

  if(!fileModel) {
    std::cout<< "File " << argv[3] << "does not exist." << std::endl;
    return 0;
  }
  if(!fileTarget) {
    std::cout << "File " << argv[2] << "does not exist." << std::endl;
    return 0;
  }

  molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
  molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
  if (molModel.size() == 0 || molTarget.size() == 0){
      std::ifstream fileModel(argv[3]);
      std::ifstream fileTarget(argv[2]);
      molModel.readPDBfile(fileModel, PDB::PSelector());
      molTarget.readPDBfile(fileTarget, PDB::PSelector());
  }

  // next we insert the target molecule into hash
  // this will help us to find atoms that are close faster
  GeomHash <Vector3,int> gHash(3, epsilon); // 3 is a dimension and epsilon is the size of the hash cube
  for(unsigned int i=0; i<molTarget.size(); i++) {
    gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
  }

  // now we try random rotations and choose the best alignment from random rotations
  unsigned int iMaxSize=0;
  RigidTrans3 rtransBest;
  float rmsd = 0;
  for (int i_model = 0; i_model < molModel.size() - 2; ++i_model) {
    Triangle* modelTriangel = new Triangle(molModel[i_model].position(), molModel[i_model + 1].position(),
            molModel[i_model + 2].position());
    for (int j_target = 0; j_target < molTarget.size() - 2; ++j_target) {
        Triangle* targetTriangle = new Triangle(molTarget[j_target].position(), molTarget[j_target + 1].position(),
                                                  molTarget[j_target + 2].position());

        RigidTrans3 rotation = *targetTriangle | *modelTriangel;
      // match is a class that stores the correspondence list, eg.
      // pairs of atoms, one from each molecule, that are matching
        Match match;

          // apply rotation on each atom in the model molecule and
          // add the pairs of atoms (one from target and one from model)
          // that are close enough to the match list
          for(unsigned int i=0; i< molModel.size(); i++) {
              Vector3 mol_atom = rotation*molModel[i].position(); // rotate

              // find close target molecule atoms using the hash
              HashResult<int> result;
              gHash.query(mol_atom, epsilon, result); // key is mol atom coordinate

              // check if the atoms in the result are inside the distance threshold
              // the hash is a cube shape, there can be atoms further that the threshold
              for(auto x = result.begin(); x != result.end(); x++) {
                  float dist = mol_atom.dist(molTarget[*x].position());
                  if(dist <= epsilon) {
                      float score = (1 / (1 + dist));
                      match.add( *x , i, score, score );
                  }
              }
              result.clear();
          }

          //calculates transformation that is a little better than "rotation"
          match.calculateBestFit(molTarget, molModel);

          if(iMaxSize < match.size() ){
              iMaxSize = match.size();
              rtransBest=match.rigidTrans();
              rmsd=match.rmsd();
          }

          }
      }

    // printing the values
    std::cout <<  iMaxSize <<" "<< rmsd<< " " << rtransBest << std::endl;

    std::ifstream allModel(argv[3]);
    allModelAtoms.readPDBfile(allModel);
    allModelAtoms *= rtransBest;

    // writing to transformed file
    std::ofstream myfile;
    myfile.open ("transformed.pdb");
    myfile << allModelAtoms;
    myfile.close();

}
