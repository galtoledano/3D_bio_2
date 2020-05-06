////
//// Created by galtoledano on 05/05/2020.
////
//
//#include "structalign.h"
//#include "Vector3.h"
//#include "Atom.h"
//#include "RigidTrans3.h"
//#include "Matrix3.h"
//#include "Molecule.h"
//#include "PDB.h"
//#include "Match.h"
//#include "GeomHash.h"
//#include "Triangle.h"
//#include <chrono>
//#include <iostream>
//
//int main(int argc , char* argv[]){
//
//    // measure the run time
//    auto start = std::chrono::system_clock::now();
//
//
//    if(argc != 4) {
//        std::cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << std::endl;
//        exit(1);
//    }
//
//    //********Parameters********************
//    float epsilon = atof(argv[1]); // distance threshold on atoms in correspondence
//
//
//    std::cout << "Distance threshold: "<< epsilon  << std::endl;
//
//
//
//    std::ifstream fileModel(argv[3]);
//    std::ifstream fileTarget(argv[2]);
//
//    if(!fileModel) {
//        std::cout<< "File " << argv[3] << " does not exist." << std::endl;
//        return 0;
//    }
//    if(!fileTarget) {
//        std::cout << "File " << argv[2] << " does not exist." << std::endl;
//        return 0;
//    }
//
//    // read the two files into Molecule
//    Molecule<Atom> molModel, molTarget;
//
//    molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
//    molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
//
//    if(molModel.size() == 0 || molTarget.size() == 0){
//        std::ifstream fileModel(argv[3]);
//        std::ifstream fileTarget(argv[2]);
//        molModel.readPDBfile(fileModel, PDB::PSelector());
//        molTarget.readPDBfile(fileTarget, PDB::PSelector());
//    }
//
//    int c = 0;
//    // next we insert the target molecule into hash
//    // this will help us to find atoms that are close faster
//    GeomHash <Vector3,int> gHash(3, epsilon); // 3 is a dimension and epsilon is the size of the hash cube
//    for(unsigned int i=0; i<molTarget.size(); i++) {
//        gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
//        c++;
//    }
//    printf("hash counter : %d \n", c);
//
//    unsigned int iMaxSize=0;
//    RigidTrans3 bestRigid;
//    float rmsd = 0;
//    int counter =0;
//
//    for (unsigned int i = 0; i < molModel.size() - 2; i++) {
//        Triangle triangleModel(molModel[i].position(), molModel[i + 1].position(),
//                               molModel[i + 2].position());
//        for (unsigned int j = 0; j < molTarget.size() - 2; j++) {
//            Triangle triangleTarget(molTarget[j].position(), molTarget[j + 1].position(),
//                                   molTarget[j + 2].position());
//            RigidTrans3 currentRigid = triangleModel | triangleTarget;
//
//            Match match;
//            for(unsigned int k=0; k < molModel.size(); k++) {
//                Vector3 mol_atom = currentRigid * molModel[k].position();
//
//                // find close target molecule atoms using the hash
//                HashResult<int> result;
//                gHash.query(mol_atom, epsilon, result); // key is mol atom coordinate
//
////                printf("result : %d\n", result.size());
//                // check if the atoms in the result are inside the distance threshold
//                // the hash is a cube shape, there can be atoms further that the threshold
//                for (auto x = result.begin(); x != result.end(); x++) {
//                    float dist = mol_atom.dist(molTarget[*x].position());
//                    if (dist <= epsilon) {
//                        float score = (1 / (1 + dist));
//                        match.add(*x, i, score, score);
//                        counter ++;
////                        std::cout<< "match ? "  <<std::endl;
//                    }
//                }
////                printf("counter : %d\n", counter);
//                counter = 0;
//                result.clear();
//            }
//                //calculates transformation that is a little better than "rotation"
//                match.calculateBestFit(molTarget, molModel);
//                if(iMaxSize < match.size() ){
//                printf("match size: %d \n", match.size());
//                    iMaxSize = match.size();
//                    bestRigid = match.rigidTrans();
//                    rmsd = match.rmsd();
//                }
//          }
//        }
//
//    std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
//    std::cout << "rmsd: " << rmsd << std::endl;
//    std::cout << "counter: " << counter << std::endl;
////    std::cout << "Rigid Trans: " <<
////              RigidTrans3(Vector3(0,0,0),vectTargetMass)*
////              rtransBest*
////              RigidTrans3(Vector3(0,0,0),(-vectModelMass)) << std::endl;
//
//    auto end = std::chrono::system_clock::now();
//
//    std::chrono::duration<double> elapsed_seconds = end-start;
//    std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
//}
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
  // measure the run time
  auto start = std::chrono::system_clock::now();

  if(argc !=4) {
    std::cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << std::endl;
    exit(1);
  }

  //********Parameters********************
  float m_fDistThr = atof(argv[1]); // distance threshold on atoms in correspondence

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
  GeomHash <Vector3,int> gHash(3,m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
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
              gHash.query(mol_atom, m_fDistThr, result); // key is mol atom coordinate

              // check if the atoms in the result are inside the distance threshold
              // the hash is a cube shape, there can be atoms further that the threshold
              for(auto x = result.begin(); x != result.end(); x++) {
                  float dist = mol_atom.dist(molTarget[*x].position());
                  if(dist <= m_fDistThr) {
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

    std::cout <<  iMaxSize <<" "<< rmsd<< " " << rtransBest << std::endl;
    std::ifstream allModel(argv[3]);
    allModelAtoms.readPDBfile(allModel);
    allModelAtoms *= rtransBest;
    std::ofstream myfile;
    myfile.open ("transformed.pdb");
    myfile << allModelAtoms;
    myfile.close();

}
