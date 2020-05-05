//
// Created by galtoledano on 05/05/2020.
//

#include "structalign.h"
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


    if(argc != 4) {
        std::cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << std::endl;
        exit(1);
    }

    //********Parameters********************
    float epsilon = atof(argv[1]); // distance threshold on atoms in correspondence


    std::cout << "Distance threshold: "<< epsilon  << std::endl;



    std::ifstream fileModel(argv[3]);
    std::ifstream fileTarget(argv[2]);

    if(!fileModel) {
        std::cout<< "File " << argv[3] << " does not exist." << std::endl;
        return 0;
    }
    if(!fileTarget) {
        std::cout << "File " << argv[2] << " does not exist." << std::endl;
        return 0;
    }

    // read the two files into Molecule
    Molecule<Atom> molModel, molTarget;

    molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
    molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());

    if(molModel.size() == 0 || molTarget.size() == 0){
//        std::ifstream fileModel(argv[3]);
        molModel.readPDBfile(fileModel, PDB::PSelector());
        molTarget.readPDBfile(fileTarget, PDB::PSelector());
    }

    // next we insert the target molecule into hash
    // this will help us to find atoms that are close faster
    GeomHash <Vector3,int> gHash(3,epsilon); // 3 is a dimension and m_fDistThr is the size of the hash cube
    for(unsigned int i=0; i<molTarget.size(); i++) {
        gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
    }

    RigidTrans3 bestRigid;
    unsigned int iMaxSize=0;
    float rmsd = 0;
    int counter =0;

    for (int i = 0; i < molModel.size()-2; ++i) {
        auto triangleModel = new Triangle(molModel[i].position(), molModel[i + 1].position(),
                                          molModel[i + 2].position());
        for (int j = 0; j < molTarget.size() - 2; ++j) {
            auto triangleTarget = new Triangle(molTarget[j].position(), molTarget[j + 1].position(),
                                               molTarget[j + 2].position());
            auto currentRigid = *triangleModel | *triangleTarget;

            Match match;
            for(unsigned int k=0; k < molModel.size(); k++) {
                Vector3 mol_atom = currentRigid * molModel[k].position();

                // find close target molecule atoms using the hash
                HashResult<int> result;
                gHash.query(mol_atom, epsilon, result); // key is mol atom coordinate

                // check if the atoms in the result are inside the distance threshold
                // the hash is a cube shape, there can be atoms further that the threshold
                for (auto x = result.begin(); x != result.end(); x++) {
                    float dist = mol_atom.dist(molTarget[*x].position());
                    if (dist <= epsilon) {
                        float score = (1 / (1 + dist));
                        match.add(*x, i, score, score);
                        counter ++;
//                        std::cout<< "match ? " << a <<std::endl;
                    }
                }
                printf("counter : %d\n", counter);
                counter = 0;
                result.clear();
            }
                //calculates transformation that is a little better than "rotation"
                match.calculateBestFit(molTarget, molModel);
                printf("match size: %d \n ", match.size());
                if(iMaxSize < match.size() ){
                    iMaxSize = match.size();
                    bestRigid = match.rigidTrans();
                    rmsd = match.rmsd();
                }
          }
        }

    std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
    std::cout << "rmsd: " << rmsd << std::endl;
    std::cout << "counter: " << counter << std::endl;
//    std::cout << "Rigid Trans: " <<
//              RigidTrans3(Vector3(0,0,0),vectTargetMass)*
//              rtransBest*
//              RigidTrans3(Vector3(0,0,0),(-vectModelMass)) << std::endl;

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
}

